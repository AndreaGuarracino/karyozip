//! --db: the labels are featureIDs from a database, so one merged track
//! replaces the per-set BEDs and the mapping stays in features.tsv.

use crate::kst::{self, FNV0};
use std::collections::HashMap;
use std::io::Write;

pub struct Feats {
    pub sets: Vec<Vec<u8>>,
    /// featureID to its label per set.
    pub id2lab: HashMap<u32, Vec<Vec<u8>>>,
    /// the label tuple, tab-joined, to its featureID
    tup: HashMap<Vec<u8>, u32>,
    pub hash: u64,
}

impl Feats {
    pub fn load(path: &str) -> Feats {
        let data = std::fs::read(path)
            .unwrap_or_else(|_| crate::bail("cannot open the features table"));
        let mut f = Feats { sets: Vec::new(), id2lab: HashMap::new(),
                            tup: HashMap::new(), hash: FNV0 };
        let mut first = true;
        for raw in data.split_inclusive(|&c| c == b'\n') {
            f.hash = kst::fnv(f.hash, raw);
            let mut line = raw.strip_suffix(b"\n").unwrap_or(raw);
            if line.last() == Some(&b'\r') { line = &line[..line.len() - 1]; }
            let cols: Vec<&[u8]> = line.split(|&c| c == b'\t').collect();
            if first {
                if cols.len() < 2 {
                    crate::bail("features.tsv: the first line must be the tab-separated header");
                }
                for s in &cols[1..] {
                    if s.is_empty() || s.contains(&b'/') {
                        crate::bail("features.tsv: a set name is empty or has a slash");
                    }
                    if f.sets.iter().any(|x| x == s) {
                        crate::bail("features.tsv names a set twice");
                    }
                    f.sets.push(s.to_vec());
                }
                first = false;
                continue;
            }
            if cols.len() < 2 { if !line.is_empty() { crate::bail("features.tsv: a row without a tab"); } continue; }
            if cols[0].first() == Some(&b'0') && cols[0].len() > 1 {
                crate::bail("features.tsv: a featureID with a leading zero");
            }
            let id = num(cols[0], (1u32 << 26) - 1,
                         "features.tsv: a featureID is not a number below 2^26");
            if id == 0 { crate::bail("features.tsv defines featureID 0, which is reserved for novel"); }
            if cols.len() - 1 < f.sets.len() { crate::bail("a features.tsv row has fewer labels than sets"); }
            if cols.len() - 1 > f.sets.len() { crate::bail("a features.tsv row has more labels than sets"); }
            if f.id2lab.contains_key(&id) { crate::bail("features.tsv defines a featureID twice"); }
            f.tup.insert(cols[1..].join(&b'\t'), id);
            f.id2lab.insert(id, cols[1..].iter().map(|c| c.to_vec()).collect());
        }
        if f.sets.is_empty() { crate::bail("features table has no header"); }
        f
    }
}

fn num(b: &[u8], max: u32, what: &str) -> u32 {
    if b.is_empty() { crate::bail(what); }
    let mut v: u64 = 0;
    for &c in b {
        if !c.is_ascii_digit() { crate::bail(what); }
        v = v * 10 + (c - b'0') as u64;
        if v > max as u64 { crate::bail(what); }
    }
    v as u32
}

/// One input BED, read a row at a time.
struct Brd { lines: Vec<Vec<u8>>, i: usize, seq: Vec<u8>, st: u64, en: u64, lab: Vec<u8>, eof: bool }

impl Brd {
    fn open(path: &str) -> Brd {
        let data = if path == "-" {
            use std::io::Read;
            let mut v = Vec::new();
            std::io::stdin().read_to_end(&mut v).unwrap_or_else(|_| crate::bail("cannot open an input BED"));
            v
        } else {
            std::fs::read(path).unwrap_or_else(|_| {
                eprintln!("cannot open {path}");
                crate::bail("cannot open an input BED")
            })
        };
        let lines = data.split(|&c| c == b'\n').filter(|l| !l.is_empty())
            .map(|l| l.to_vec()).collect();
        let mut b = Brd { lines, i: 0, seq: Vec::new(), st: 0, en: 0, lab: Vec::new(), eof: false };
        b.next();
        if b.eof { crate::bail("an input BED is empty"); }
        b
    }

    fn next(&mut self) {
        if self.i >= self.lines.len() { self.eof = true; return; }
        let line = std::mem::take(&mut self.lines[self.i]);
        self.i += 1;
        if line.contains(&0) { crate::bail("a NUL byte inside a row"); }
        let c: Vec<&[u8]> = line.splitn(4, |&x| x == b'\t').collect();
        if c.len() < 4 { crate::bail("bad line"); }
        if c[0].len() >= 512 { crate::bail("name too long"); }
        self.seq = c[0].to_vec();
        self.st = crate::bedcoord(c[1], "start is not a number");
        self.en = crate::bedcoord(c[2], "end is not a number");
        self.lab = c[3].to_vec();
    }
}

/// One BED per feature set in header order, or one merged featureID track.
pub fn encode(files: &[String], feats: &Feats, w: &mut crate::writer::Writer) {
    let dbids = files.len() != feats.sets.len();
    if dbids && files.len() != 1 {
        crate::bail("--db: give one BED per feature set in features.tsv column order, or one merged featureID track");
    }
    let mut b: Vec<Brd> = files.iter().map(|f| Brd::open(f)).collect();
    let (mut rows, mut missing) = (0u64, 0u64);
    let (mut prevfid, mut lastend) = (u32::MAX, 0u64);
    let mut lastseq: Vec<u8> = Vec::new();
    let mut haveseq = false;
    while !b[0].eof {
        for i in 1..b.len() {
            if b[i].eof || b[i].seq != b[0].seq {
                crate::bail("the input BEDs disagree on sequence order");
            }
        }
        let lo = b[0].st;
        let mut hi = b[0].en;
        for i in 1..b.len() {
            if b[i].st != lo { crate::bail("the input BEDs disagree on a run start"); }
            if b[i].en < hi { hi = b[i].en; }        // cut at the first change
        }
        if !haveseq || lastseq != b[0].seq {
            haveseq = true;
            lastseq = b[0].seq.clone();
            prevfid = u32::MAX;                      // a new sequence starts fresh
        }
        let id = if dbids {
            if b[0].lab.first() == Some(&b'0') && b[0].lab.len() > 1 {
                crate::bail("--db: a featureID with a leading zero");
            }
            let v = num(&b[0].lab, (1u32 << 26) - 1,
                        "--db: a label is not a featureID; give one BED per feature set instead");
            if v != 0 && !feats.id2lab.contains_key(&v) {
                crate::bail("--db: featureID not in the features table");
            }
            v
        } else {
            let key = b.iter().map(|x| x.lab.clone()).collect::<Vec<_>>().join(&b'\t');
            match feats.tup.get(&key) {
                Some(&v) => v,
                None => {
                    if !b.iter().all(|x| x.lab == b"novel") {
                        crate::bail("a label tuple is not in the features table");
                    }
                    missing += 1;
                    0
                }
            }
        };
        if prevfid != u32::MAX && id == prevfid && lastend == lo {
            crate::bail("two touching runs carry the same featureID: merge them first");
        }
        prevfid = id;
        lastend = hi;
        let mut lab = Vec::with_capacity(12);
        kst::putu64(&mut lab, id as u64);
        w.row(&b[0].seq, hi - lo, &lab, lo)
            .unwrap_or_else(|m| crate::bail(&m));
        rows += 1;
        for x in b.iter_mut() {
            x.st = hi;
            if x.st == x.en { x.next(); }
        }
    }
    for i in 1..b.len() {
        if !b[i].eof { crate::bail("the input BEDs have different lengths"); }
    }
    eprintln!("db: {rows} rows from {} tracks, {missing} all-novel", b.len());
}

/// decode --db --out DIR: one BED per feature set, the featureIDs expanded.
pub fn decode(r: &mut crate::reader::Reader, feats: &Feats, dir: &str) {
    let n = feats.sets.len();
    std::fs::create_dir_all(dir).unwrap_or_else(|_| crate::bail("cannot create the output directory"));
    let mut out: Vec<std::io::BufWriter<std::fs::File>> = feats.sets.iter()
        .map(|s| {
            let p = format!("{dir}/{}.bed", String::from_utf8_lossy(s));
            std::io::BufWriter::new(std::fs::File::create(&p)
                .unwrap_or_else(|_| crate::bail("cannot write a derived BED")))
        }).collect();
    // One open run per set, extended while its label does not change.
    let mut cur: Vec<Option<(u64, u64, Vec<u8>)>> = vec![None; n];
    let mut cseq: Vec<u8> = Vec::new();
    let mut rows = 0u64;
    let novel = b"novel".to_vec();
    r.for_each_row(&mut |seq: &[u8], st: u64, en: u64, lab: &[u8]| {
        if seq != cseq {
            for j in 0..n {
                if let Some((a, b2, l)) = cur[j].take() {
                    wr(&mut out[j], &cseq, a, b2, &l);
                }
            }
            cseq = seq.to_vec();
        }
        let id = num(lab, (1u32 << 26) - 1, "corrupt --db track: a label is not a featureID");
        if id != 0 && !feats.id2lab.contains_key(&id) {
            crate::bail("a featureID in the track is not in the features table");
        }
        for j in 0..n {
            let l: &[u8] = if id == 0 { &novel } else { &feats.id2lab[&id][j] };
            match &mut cur[j] {
                Some((_, e, cl)) if *e == st && cl == l => *e = en,
                slot => {
                    if let Some((a, b2, cl)) = slot.take() { wr(&mut out[j], &cseq, a, b2, &cl); }
                    *slot = Some((st, en, l.to_vec()));
                }
            }
        }
        rows += 1;
    });
    for j in 0..n {
        if let Some((a, b2, l)) = cur[j].take() { wr(&mut out[j], &cseq, a, b2, &l); }
        out[j].flush().unwrap_or_else(|_| crate::bail("write failed"));
    }
    eprintln!("db: {rows} rows into {n} tracks under {dir}");
}

fn wr(f: &mut impl Write, seq: &[u8], st: u64, en: u64, lab: &[u8]) {
    let mut v = Vec::with_capacity(seq.len() + lab.len() + 24);
    v.extend_from_slice(seq);
    v.push(b'\t');
    crate::kst::putu64(&mut v, st);
    v.push(b'\t');
    crate::kst::putu64(&mut v, en);
    v.push(b'\t');
    v.extend_from_slice(lab);
    v.push(b'\n');
    f.write_all(&v).unwrap_or_else(|_| crate::bail("write failed"));
}
