//! --db: the labels are featureIDs from a database, so one merged track
//! replaces the per-set BEDs and the mapping stays in features.tsv.

use crate::kst::{self, FNV0, FnvMap as HashMap};
use std::io::Write;

pub struct Feats {
    pub sets: Vec<Vec<u8>>,
    /// featureID to its label per set: a dense position table, as the C's
    /// fid2lab, 4 bytes a slot so a sparse id costs less than the C's 8.
    idpos: Vec<u32>,
    labs: Vec<Vec<Vec<u8>>>,
    pub nids: usize,
    /// the label tuple, tab-joined, to its featureID
    tup: HashMap<Vec<u8>, u32>,
    /// the C's integer-tuple index: labels interned per set, tuples of ids
    intern: Vec<HashMap<Vec<u8>, u32>>,
    tupint: HashMap<kst::Tup, u32>,
    pub hash: u64,
    /// --hier: per set, features.tsv column plus the hierarchy nodes plus
    /// "novel", sorted; a tuple outside features.tsv is 2^26 + its rank
    voc: Vec<Vec<Vec<u8>>>,
    stride: Vec<u64>,
    pub prod: u64,
    pub cart: bool,
    pub hierhash: u64,
}
pub const CARTBASE: u32 = 1 << 26;

impl Feats {
    pub fn load(path: &str) -> Feats {
        let data = crate::slurp(path, "cannot open the features table");
        let mut f = Feats { sets: Vec::new(), idpos: Vec::new(), labs: Vec::new(), nids: 0,
                            tup: HashMap::default(), intern: Vec::new(), tupint: HashMap::default(), hash: FNV0, voc: Vec::new(), stride: Vec::new(),
                            prod: 0, cart: false, hierhash: 0 };
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
            if f.has(id) { crate::bail("features.tsv defines a featureID twice"); }
            f.tup.insert(cols[1..].join(&b'\t'), id);
            if f.idpos.len() <= id as usize {
                let n = id as usize + 1024 - f.idpos.len();
                f.idpos.try_reserve(n).unwrap_or_else(|_| crate::bail("oom"));
                f.idpos.resize(id as usize + 1024, u32::MAX);
            }
            f.idpos[id as usize] = f.labs.len() as u32;
            f.labs.push(cols[1..].iter().map(|c| c.to_vec()).collect());
            f.nids += 1;
            if f.intern.is_empty() { f.intern = vec![HashMap::default(); f.sets.len()]; }
            let key = kst::Tup(cols[1..].iter().enumerate().map(|(k, c)| {
                let n = f.intern[k].len() as u32; *f.intern[k].entry(c.to_vec()).or_insert(n) }).collect());
            f.tupint.entry(key).or_insert(id);
        }
        if f.sets.is_empty() { crate::bail("features table has no header"); }
        f
    }
    fn has(&self, id: u32) -> bool { self.idpos.get(id as usize).is_some_and(|&p| p != u32::MAX) }
    /// A label's interned id in set k, or MAX when no table tuple has it.
    fn labid(&self, k: usize, lab: &[u8]) -> u32 { self.intern[k].get(lab).copied().unwrap_or(u32::MAX) }

    /// A track label as a featureID, checked as decode checks it.
    pub fn featid(&self, lab: &[u8]) -> u32 {
        let id = num(lab, self.idmax(), "corrupt --db track: a label is not a featureID");
        if id != 0 && id < CARTBASE && !self.has(id) {
            crate::bail("a featureID in the track is not in the features table");
        }
        id
    }
    /// The largest featureID the track may hold.
    pub fn idmax(&self) -> u32 { if self.cart { CARTBASE - 1 + self.prod as u32 } else { CARTBASE - 1 } }

    /// The label of one set for a featureID.
    pub fn lab(&self, id: u32, k: usize) -> &[u8] {
        if id == 0 { b"novel" }
        else if id < CARTBASE { &self.labs[self.idpos[id as usize] as usize][k] }
        else { &self.voc[k][((id - CARTBASE) as u64 / self.stride[k] % self.voc[k].len() as u64) as usize] }
    }

    /// The vocabulary of every set, and the hash continued over the file.
    pub fn load_hier(&mut self, path: &str) {
        let data = crate::slurp(path, "cannot open the hierarchy file");
        let n = self.sets.len();
        let mut voc: Vec<Vec<Vec<u8>>> = (0..n).map(|k| {
            let mut v = vec![b"novel".to_vec()];
            v.extend(self.labs.iter().map(|t| t[k].clone()));
            v
        }).collect();
        self.hierhash = self.hash;
        for raw in data.split_inclusive(|&c| c == b'\n') {
            self.hierhash = kst::fnv(self.hierhash, raw);
            let mut line = raw.strip_suffix(b"\n").unwrap_or(raw);
            if line.last() == Some(&b'\r') { line = &line[..line.len() - 1]; }
            if line.is_empty() { continue; }
            let cols: Vec<&[u8]> = line.split(|&c| c == b'\t').collect();
            if cols.len() < 3 { continue; }
            if cols.len() > 3 { crate::bail("hierarchy: a row has more than three columns"); }
            if cols[0] == b"feature_set" { continue; }
            if let Some(k) = self.sets.iter().position(|s| s == cols[0]) {
                voc[k].push(cols[1].to_vec()); voc[k].push(cols[2].to_vec());
            }
        }
        self.stride = vec![0; n];
        self.prod = 1;
        for k in (0..n).rev() {
            voc[k].sort(); voc[k].dedup();
            self.stride[k] = self.prod;
            let u = voc[k].len() as u64;
            if self.prod > (u32::MAX - CARTBASE) as u64 / u { crate::bail("the label vocabulary is too large for featureIDs"); }
            self.prod *= u;
        }
        self.voc = voc;
        self.cart = true;
        eprintln!("hierarchy: {} label tuples possible", self.prod);
    }

    /// 2^26 + the mixed-radix rank of a tuple outside features.tsv.
    fn cart_id(&self, labs: &[&[u8]]) -> u32 {
        let mut rank = 0u64;
        for (k, l) in labs.iter().enumerate() {
            let x = self.voc[k].binary_search_by(|v| v.as_slice().cmp(l))
                .unwrap_or_else(|_| crate::bail("a label is neither in features.tsv nor in the hierarchy of its set"));
            rank += x as u64 * self.stride[k];
        }
        CARTBASE + rank as u32
    }

    /// --set A,B as columns of features.tsv; None or "all" is every set.
    pub fn setcols(&self, spec: Option<&str>) -> Vec<usize> {
        let Some(spec) = spec.filter(|s| *s != "all") else { return (0..self.sets.len()).collect() };
        let mut cols = Vec::new();
        for t in spec.split(',').filter(|t| !t.is_empty()) {
            if t == "all" { crate::bail("--set all takes no other set"); }
            let k = self.sets.iter().position(|s| s == t.as_bytes())
                .unwrap_or_else(|| crate::bail(&format!("--set: {t} is not a set in features.tsv")));
            if cols.contains(&k) { crate::bail("--set names a set twice"); }
            cols.push(k);
        }
        if cols.is_empty() { crate::bail("--set needs at least one set"); }
        cols
    }

    /// The labels of the picked sets for one featureID, tab-joined.
    pub fn joined(&self, id: u32, cols: &[usize], out: &mut Vec<u8>) {
        out.clear();
        for (i, &k) in cols.iter().enumerate() {
            if i > 0 { out.push(b'\t'); }
            out.extend_from_slice(self.lab(id, k));
        }
    }

    pub fn check_hash(&self, r: &crate::reader::Reader) {
        if r.cart() && !self.cart { crate::bail("this track holds label tuples outside features.tsv: give --hier HIERARCHY"); }
        if r.feathash() != if r.cart() { self.hierhash } else { self.hash } {
            crate::bail(if r.cart() { "this features.tsv and hierarchy are not the ones the track was encoded with" }
                        else { "this features.tsv is not the one the track was encoded with" });
        }
    }
}

/// Digits, in 64 bits, at most `max`: the C's num_arg, same words in the same order.
fn num(b: &[u8], max: u32, what: &str) -> u32 {
    if b.is_empty() { crate::bail(what); }
    let mut v: u64 = 0;
    for &c in b {
        if !c.is_ascii_digit() { crate::bail(what); }
        if v > (u64::MAX - (c - b'0') as u64) / 10 { crate::bail("a coordinate does not fit in 64 bits"); }
        v = v * 10 + (c - b'0') as u64;
    }
    if v > max as u64 { crate::bail(what); }
    v as u32
}

/// One input BED, streamed a row at a time.
struct Brd { lines: crate::Lines, seq: Vec<u8>, st: u64, en: u64, lab: Vec<u8>, labid: u32,
             pseq: Vec<u8>, pen: u64, had: bool, eof: bool }

impl Brd {
    fn open(path: &str) -> Brd {
        let lines = if path == "-" { crate::Lines::new(std::io::stdin()) } else {
            if std::fs::metadata(path).is_err() { eprintln!("cannot open {path}"); }
            crate::Lines::open(path, "cannot open an input BED")
        };
        let mut lines = lines; lines.err = "read error on an input BED";
        let mut b = Brd { lines, seq: Vec::new(), st: 0, en: 0, lab: Vec::new(), labid: u32::MAX,
                          pseq: Vec::new(), pen: 0, had: false, eof: false };
        b.next();
        if b.eof { crate::bail("an input BED is empty"); }
        b
    }

    /// The row parsed straight from the reader's buffer, no copy of the line,
    /// and the previous label compared, not copied: as the C's bnext.
    fn next(&mut self) {
        let Some((line, nl)) = self.lines.next_line_nl() else { self.eof = true; return; };
        if line.is_empty() { crate::bail("bad line"); }   // as the C: a blank row is not skipped
        if !nl { crate::bail("input does not end with a newline"); }
        if crate::ffi::find(0, line).is_some() { crate::bail("a NUL byte inside a row"); }
        let mut it = line.splitn(4, |&x| x == b'\t');
        let mut c = [&line[..0]; 4];
        for (i, s) in c.iter_mut().enumerate() { *s = it.next().unwrap_or_else(|| { if i < 4 { crate::bail("bad line") } else { &[] } }); }
        if c[0].len() >= 512 { crate::bail("name too long"); }
        if c[0].is_empty() { crate::bail("empty sequence name"); }
        if self.had {
            std::mem::swap(&mut self.pseq, &mut self.seq);
            self.pen = self.en;
        }
        self.seq.clear(); self.seq.extend_from_slice(c[0]);
        self.st = crate::bedcoord(c[1], "start is not a number");
        self.en = crate::bedcoord(c[2], "end is not a number");
        if c[1..3].iter().any(|v| v.len() > 1 && v[0] == b'0') {
            crate::bail("a coordinate has a leading zero: it would not decode to the same bytes");
        }
        if self.en < self.st { crate::bail("end is before start"); }
        if self.en == self.st {
            crate::bail("end equals start: a run must cover at least one base");
        }
        let samelab = c[3] == self.lab.as_slice();
        if !samelab { self.lab.clear(); self.lab.extend_from_slice(c[3]); self.labid = u32::MAX - 1; }
        if self.had && self.st == self.pen && self.seq == self.pseq && samelab {
            crate::bail("two touching runs carry the same label: merge them first");
        }
        self.had = true;
    }
}

/// One BED per feature set in header order, or one merged featureID track.
/// True when some featureID came from the vocabulary.
pub fn encode(files: &[String], feats: &Feats, w: &mut crate::writer::Writer) -> bool {
    let mut cart_used = false;
    let dbids = files.len() != feats.sets.len();
    if dbids && files.len() != 1 {
        crate::bail("--db: give one BED per feature set in features.tsv column order, or one merged featureID track");
    }
    let mut b: Vec<Brd> = files.iter().map(|f| Brd::open(f)).collect();
    let (mut rows, mut missing) = (0u64, 0u64);
    let (mut prevfid, mut lastend) = (u32::MAX, 0u64);
    let mut lastseq: Vec<u8> = Vec::new();
    let mut haveseq = false;
    let mut keyi = kst::Tup(Vec::with_capacity(b.len()));
    let mut lab = Vec::with_capacity(12);
    while !b[0].eof {
        for i in 1..b.len() {
            if b[i].eof || b[i].seq != b[0].seq {
                crate::bail("the input BEDs disagree on sequence order");
            }
        }
        let lo = b[0].st;
        let mut hi = b[0].en;
        for x in &b[1..] {
            if x.st != lo { crate::bail("the input BEDs disagree on a run start"); }
            if x.en < hi { hi = x.en; }              // cut at the first change
        }
        if !haveseq || lastseq != b[0].seq {
            haveseq = true;
            lastseq.clone_from(&b[0].seq);   // reuses the buffer, one a sequence
            prevfid = u32::MAX;                      // a new sequence starts fresh
        }
        for (k, x) in b.iter_mut().enumerate() { if x.labid == u32::MAX - 1 { x.labid = feats.labid(k, &x.lab); } }
        let id = if dbids {
            if b[0].lab.first() == Some(&b'0') && b[0].lab.len() > 1 {
                crate::bail("--db: a featureID with a leading zero");
            }
            let v = num(&b[0].lab, feats.idmax(),
                        "--db: a label is not a featureID; give one BED per feature set instead");
            if v != 0 && v < CARTBASE && !feats.has(v) {
                crate::bail("--db: featureID not in the features table");
            }
            if v >= CARTBASE { cart_used = true; }
            v
        } else {
            keyi.0.clear(); keyi.0.extend(b.iter().map(|x| x.labid));
            let hit = if b.iter().all(|x| x.labid != u32::MAX) { feats.tupint.get(&keyi).copied() } else { None };
            match hit.or_else(|| feats.tup.get(&b.iter().map(|x| x.lab.clone()).collect::<Vec<_>>().join(&b'\t')).copied()) {
                Some(v) => v,
                None if b.iter().all(|x| x.lab == b"novel") => { missing += 1; 0 }
                None if !feats.cart => crate::bail("a label tuple is not in the features table"),
                None => { cart_used = true; feats.cart_id(&b.iter().map(|x| x.lab.as_slice()).collect::<Vec<_>>()) }
            }
        };
        if prevfid != u32::MAX && id == prevfid && lastend == lo {
            crate::bail("two touching runs carry the same featureID: merge them first");
        }
        prevfid = id;
        lastend = hi;
        lab.clear();
        kst::putu64(&mut lab, id as u64);
        w.push(&b[0].seq, hi - lo, &lab, lo, false)
            .unwrap_or_else(|m| crate::bail(&m));
        rows += 1;
        for x in b.iter_mut() {
            x.st = hi;
            if x.st == x.en { x.next(); }
        }
    }
    if b[1..].iter().any(|x| !x.eof) { crate::bail("the input BEDs have different lengths"); }
    eprintln!("db: {rows} rows from {} tracks, {missing} all-novel", b.len());
    cart_used
}

/// decode --db --out-dir DIR: one BED per feature set, the featureIDs expanded.
/// The staged BEDs of a decode in flight, removed by `die` and on the way
/// out, so a failed run never leaves half an output behind.
static TMP: std::sync::Mutex<Vec<std::path::PathBuf>> = std::sync::Mutex::new(Vec::new());
/// DIR/<set><suffix>, from the set's bytes: a set name need not be UTF-8,
/// and two that differ only past a lossy conversion must not share a file.
fn setpath(dir: &str, s: &[u8], suffix: &str) -> std::path::PathBuf {
    use std::os::unix::ffi::OsStrExt;
    let mut f = s.to_vec(); f.extend_from_slice(suffix.as_bytes());
    std::path::Path::new(dir).join(std::ffi::OsStr::from_bytes(&f))
}

pub fn drop_tmp() {
    crate::ffi::sig_stage(&[]);
    if let Ok(mut v) = TMP.try_lock() {
        for p in v.drain(..) { let _ = std::fs::remove_file(&p); }
    }
}

/// decode --db. With DIR: one BED per picked set, staged then swapped in
/// together. Without: the picked columns, tab-joined, merged where they all
/// agree, on stdout; one set gives that set's own track.
pub fn decode(r: &mut crate::reader::Reader, feats: &Feats, dir: Option<&str>, cols: &[usize]) {
    let sets: Vec<&[u8]> = cols.iter().map(|&k| feats.sets[k].as_slice()).collect();
    let mut tmps: Vec<std::path::PathBuf> = Vec::new();
    // each output its own 1 MB row buffer, written whole, as the C's ow:
    // a row through a BufWriter was a third copy and a virtual call a row
    let mut out: Vec<(Box<dyn Write>, Vec<u8>)> = Vec::new();
    if let Some(dir) = dir {
        tmps = stage_open(dir, &sets);
        out = tmps.iter().map(|p| (Box::new(std::fs::File::create(p)
                .unwrap_or_else(|_| crate::bail("cannot write a derived BED"))) as Box<dyn Write>, Vec::with_capacity(OCAP))).collect();
    } else {
        out.push((Box::new(std::io::stdout()), Vec::with_capacity(OCAP)));
    }
    // One open run per output, extended while its label does not change.
    let n = out.len();
    let mut joined = Vec::new();
    // (open, start, end, label) per output: the label is the table's own
    // bytes, as the C's pointer, so a run start copies nothing; on stdout
    // the joined label of the open run is kept in cj
    let mut cur: Vec<(bool, u64, u64, &[u8])> = vec![(false, 0, 0, &[][..]); n];
    let mut cj: Vec<u8> = Vec::new();
    let mut cseq: Vec<u8> = Vec::new();
    let mut rows = 0u64;
    // a label's featureID, parsed the first time its table entry is seen, as the C's fidof
    let mut fidof: Vec<u32> = Vec::new();
    r.for_each_row(&mut |seq: &[u8], st: u64, en: u64, lab: &[u8], tix: Option<usize>| {
        if seq != cseq {
            for j in 0..n {
                let c = &mut cur[j];
                if c.0 { wr(&mut out[j], &cseq, c.1, c.2, if dir.is_none() { &cj } else { c.3 }); c.0 = false; }
            }
            cseq.clear(); cseq.extend_from_slice(seq);
        }
        let id = match tix {
            Some(t) => {
                if t >= fidof.len() { fidof.resize(t + 1, u32::MAX); }
                if fidof[t] == u32::MAX { fidof[t] = feats.featid(lab); }
                fidof[t]
            }
            None => feats.featid(lab),
        };
        if dir.is_none() { feats.joined(id, cols, &mut joined); }
        for j in 0..n {
            let l: &[u8] = if dir.is_none() { &[] } else { feats.lab(id, cols[j]) };
            let c = &mut cur[j];
            if c.0 && c.2 == st && if dir.is_none() { cj == joined } else { c.3 == l } { c.2 = en; continue; }
            if c.0 { wr(&mut out[j], &cseq, c.1, c.2, if dir.is_none() { &cj } else { c.3 }); }
            *c = (true, st, en, l);
            if dir.is_none() { cj.clear(); cj.extend_from_slice(&joined); }
        }
        rows += 1;
    });
    for j in 0..n {
        let c = &mut cur[j];
        if c.0 { wr(&mut out[j], &cseq, c.1, c.2, if dir.is_none() { &cj } else { c.3 }); }
        let (f, b) = &mut out[j];
        f.write_all(b).and_then(|_| f.flush()).unwrap_or_else(|_| crate::bail("write failed"));
    }
    drop(out);
    match dir {
        None => eprintln!("db: {rows} rows, {} sets on stdout", cols.len()),
        Some(dir) => {
            stage_done(dir, &sets, &tmps);
            eprintln!("db: {rows} rows into {n} tracks under {dir}");
        }
    }
}

const OCAP: usize = 1 << 20;
/// One row into its output's buffer; the buffer goes out when the row would not fit, as the C's owput.
fn wr(o: &mut (Box<dyn Write>, Vec<u8>), seq: &[u8], st: u64, en: u64, lab: &[u8]) {
    let (f, v) = o;
    if v.len() + seq.len() + lab.len() + 48 > OCAP.max(v.capacity()) {
        f.write_all(v).unwrap_or_else(|_| crate::bail("write failed"));
        v.clear();
    }
    crate::kst::put_row(v, seq, st, en, lab);
}

/// --out-dir DIR: one staged BED a set, beside the old ones, as the C's stage_open.
pub fn stage_open(dir: &str, sets: &[&[u8]]) -> Vec<std::path::PathBuf> {
    std::fs::create_dir_all(dir).unwrap_or_else(|_| crate::bail("cannot create the output directory"));
    for s in sets {
        let p = setpath(dir, s, ".bed");
        if p.as_os_str().len() + 16 >= 4096 { crate::bail("--out-dir: the output path is too long"); }
        if let Ok(m) = std::fs::symlink_metadata(&p) {
            if !m.is_file() { crate::bail("--out-dir: a set's BED exists and is not a regular file"); }
        }
    }
    let pid = std::process::id();
    let tmps: Vec<std::path::PathBuf> = sets.iter().map(|s| setpath(dir, s, &format!(".bed.tmp.{pid}"))).collect();
    *TMP.lock().unwrap() = tmps.clone();
    crate::ffi::sig_stage(&tmps);
    tmps
}

/// The staged BEDs swapped in together, as the C's stage_swap.
pub fn stage_done(dir: &str, sets: &[&[u8]], tmps: &[std::path::PathBuf]) {
    let owned: Vec<Vec<u8>> = sets.iter().map(|s| s.to_vec()).collect();
    swap_in(dir, &owned, tmps);
    crate::ffi::sig_stage(&[]);
    TMP.lock().unwrap().clear();
}

/// Move every staged BED into place, or put back what was there.
fn swap_in(dir: &str, sets: &[Vec<u8>], tmps: &[std::path::PathBuf]) {
    let pid = std::process::id();
    let live = |s: &[u8]| setpath(dir, s, ".bed");
    let mut old: Vec<Option<std::path::PathBuf>> = vec![None; sets.len()];
    let mut bad = false;
    for (j, s) in sets.iter().enumerate() {
        if bad { break; }
        let o = setpath(dir, s, &format!(".bed.old.{pid}"));
        match std::fs::rename(live(s), &o) {
            Ok(()) => old[j] = Some(o),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(_) => bad = true,
        }
    }
    let mut moved = 0usize;
    for (j, s) in sets.iter().enumerate() {
        if bad { break; }
        if std::fs::rename(&tmps[j], live(s)).is_err() { bad = true; } else { moved = j + 1; }
    }
    if bad {
        for (j, s) in sets.iter().enumerate().take(moved) {
            let _ = std::fs::rename(live(s), &tmps[j]);
        }
        for (j, s) in sets.iter().enumerate() {
            if let Some(o) = &old[j] { let _ = std::fs::rename(o, live(s)); }
        }
        crate::bail("cannot rename a derived BED into place");
    }
    for o in old.into_iter().flatten() { let _ = std::fs::remove_file(o); }
}
