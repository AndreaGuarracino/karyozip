//! Reading a KST track. Plain tracks only for now: no --aln, no --db, no
//! sidecar, no permutation. Every refusal message is the C's, because the
//! two are meant to be interchangeable.

use std::collections::HashMap;
use crate::ffi;
use crate::kst::{self, FNV0};
use std::io::{Read, Seek, SeekFrom};

const TOKFLAG: u64 = 1 << 63;
const INTFLAG: u64 = 1 << 62;
const DROPPED: u64 = 1 << 60;
const ALNFLAG: u64 = 1 << 59;
const DBFLAG: u64 = 1 << 58;
const SIDEFLAG: u64 = 1 << 57;
const LENFLAG: u64 = 1 << 56;
const TAILFLAG: u64 = 1 << 55;
const PERMFLAG: u64 = 1 << 54;
const PERMTOK: u64 = 1 << 53;
const FLAGS: u64 = TOKFLAG | INTFLAG | DROPPED | ALNFLAG | DBFLAG | SIDEFLAG
                 | LENFLAG | TAILFLAG | PERMFLAG | PERMTOK;
const IDXMAGIC: u64 = 0x3130_5844_4954_534b;   // "KSTIDX01"
const ICHUNK: u64 = 65536;

fn ordinals(n: u64) -> Vec<u8> {
    let mut v = Vec::with_capacity(crate::ckalloc(n as usize * 8));
    for i in 0..n { crate::kst::putu64(&mut v, i); v.push(b'\n'); }
    v
}

/// --names FILE: a FASTA, a FASTQ, or one name per line. Rank N is record N.
fn read_name_file(path: &str, nreads: u64) -> Vec<u8> {
    let data = crate::slurp(path, "cannot open the file given to --names");
    // getline gives no line after a final newline, and the split does, which
    // made an @-line index land on the empty tail.
    let mut lines: Vec<&[u8]> = data.split(|&c| c == b'\n').collect();
    if data.last() == Some(&b'\n') { lines.pop(); }
    // 1 fastq, 0 fasta, 2 a plain list.
    let fq = if lines.len() >= 3 && lines[0].first() == Some(&b'@') && lines[2].first() == Some(&b'+') { 1 }
             else if lines.len() >= 2 && lines[0].first() == Some(&b'>') && lines[1].first() != Some(&b'>') { 0 }
             else { 2 };
    let mut out = Vec::new();
    let mut got = 0u64;
    for (i, ln) in lines.iter().enumerate() {
        let mut q: &[u8] = ln;
        if fq == 1 { if i % 4 != 0 || q.is_empty() { continue; } q = &q[1..]; }
        else if fq == 0 { if q.first() != Some(&b'>') { continue; } q = &q[1..]; }
        let stop: &[u8] = if fq == 2 { b"\r" } else { b" \t\r" };
        let l = q.iter().position(|c| stop.contains(c)).unwrap_or(q.len());
        if l == 0 { continue; }
        if l >= 512 { crate::bail("--names: a name is longer than any the track could hold"); }
        out.extend_from_slice(&q[..l]);
        out.push(b'\n');
        got += 1;
    }
    if got != nreads {
        eprintln!("--names: the file has {got} records, the track has {nreads} sequences. Rank N is the Nth record only when every record produced rows, so this file cannot name it.");
        crate::bail("name file does not match the track");
    }
    out
}

/// One block: where it is, both frame sizes, both decoded sizes, its rows,
/// the bases before it, and its own hash.
struct Blk { off: u64, size: u64, l: u64, raw: u64, lraw: u64, rows: u64, hash: u64 }

/// seq, start, end, label, as decode hands each row out.
pub type RowFn<'a> = dyn FnMut(&[u8], u64, u64, &[u8]) + 'a;

#[derive(Default)]
pub struct ReadOpts {
    pub raw: bool,
    pub db: bool,
    pub aln: Option<std::sync::Arc<crate::aln::Aln>>,
    pub ordnames: bool,
    pub names_file: Option<String>,
}

pub struct Reader {
    f: std::fs::File,
    pub nreads: u64,
    pub rows: u64,
    lab: Vec<Vec<u8>>,
    block: u64,
    b: Vec<Blk>,
    data0: u64,
    idcodec: i32,
    lencodec: i32,
    taildec: bool,
    tailpred: u64,
    tailskip: u64,
    /// every name, newline-free, back to back, and noff[i]..noff[i+1] is
    /// name i: decoded from nz/pz on first use
    names: std::cell::OnceCell<std::sync::Arc<(Vec<u8>, Vec<u32>)>>,
    nz: Vec<u8>,
    pz: Vec<u8>,
    names_file: Option<String>,
    ordnames: bool,
    /// side, perm_ext, dropped, int, tok, perm, perm_tok
    nf: u8,
    /// the appended index: offset of its key table, and its trailer's binding word
    idx: Option<(u64, u64)>,
    /// name -> rank, built on the first lookup, which is also where the C
    /// notices that two sequences share a name.
    nameix: std::cell::OnceCell<HashMap<Vec<u8>, u64>>,
    ibx: Vec<u8>,
    hasidx: bool,
    mhash: u64,
    cur_ri: u64,
    cur_q: usize,
    cur_pend: u64,
    cached: Option<(usize, std::rc::Rc<Vec<u8>>, usize, u64)>,
    /// where the last seq_rows stopped: (block, rows into it, cursor), so a
    /// walk in order continues instead of re-skipping the block's rows
    sc: Option<(usize, u64, Cursor)>,
    feathash: u64,
    aln: Option<std::sync::Arc<crate::aln::Aln>>,
    lenfile: bool,
    /// --raw shows what is stored, so a zero length is the answer, not a fault.
    raw: bool,
    alnq: Option<usize>,
    blk_nos: Vec<u8>,
}

fn rd<T>(r: Result<T, std::io::Error>, what: &str) -> T {
    r.unwrap_or_else(|_| crate::bail(what))
}

impl Reader {

    pub fn feathash(&self) -> u64 { self.feathash }

    pub fn open_with(path: &str, o: &ReadOpts) -> Reader {
        let mut f = rd(std::fs::File::open(path), "open");
        let mut m = [0u8; 8];
        rd(f.read_exact(&mut m), "not a KST file: shorter than the magic");
        if &m != crate::writer::MAGIC {
            let sh: String = m.iter()
                .map(|&c| if (32..127).contains(&c) { c as char } else { '.' }).collect();
            if m[..3] == crate::writer::MAGIC[..3] {
                eprintln!("karyozip: this is KST version {}, this build writes {}",
                          String::from_utf8_lossy(&m[3..]),
                          String::from_utf8_lossy(&crate::writer::MAGIC[3..]));
            } else {
                eprintln!("karyozip: file starts with \"{sh}\", not \"KST00001\"");
            }
            crate::bail("not a KST file this build can read");
        }
        let mut hb = [0u8; 128];
        rd(f.read_exact(&mut hb), "short header");
        let mut h = [0u64; 16];
        for (i, v) in h.iter_mut().enumerate() {
            *v = u64::from_le_bytes(hb[i * 8..i * 8 + 8].try_into().unwrap());
        }
        let h1_raw = h[1];
        let mut fsz = rd(f.seek(SeekFrom::End(0)), "seek");
        let mut hasidx = false;
        let mut idx = None;
        // An index trailer sits past the data: take it off, so every size
        // check below still bounds the data alone.
        if fsz >= 32 + 8 + 128 {
            let mut t = [0u8; 32];
            rd(f.seek(SeekFrom::Start(fsz - 32)), "seek");
            rd(f.read_exact(&mut t), "seek");
            let w: Vec<u64> = (0..4)
                .map(|i| u64::from_le_bytes(t[i * 8..i * 8 + 8].try_into().unwrap())).collect();
            if w[3] == IDXMAGIC && w[0] == h[0] && w[0] <= (fsz - 32) / 20
                && fsz - 32 - w[0] * 20 >= 8 + 128 {
                fsz -= 32 + w[0] * 20; hasidx = true; idx = Some((fsz, w[2]));
            }
        }
        let nml = h[1] & !FLAGS;
        if nml > fsz || h[2] > fsz || h[3] > fsz || h[6] > fsz || h[11] > fsz
            || h[14] > fsz || h[11] + h[14] > fsz
            || 8 + 128 + nml + h[2] + h[3] + h[6] + h[11] + h[14] > fsz {
            crate::bail("corrupt header: the stream sizes do not fit the file");
        }
        let was_side = h[1] & SIDEFLAG != 0;
        let dropped = h[1] & DROPPED != 0;
        if !o.raw {
            if h[1] & ALNFLAG != 0 && o.aln.is_none() {
                crate::bail("this track was written with --aln: decode it with --aln PAF --ref R, or --raw to see what is stored");
            }
            if h[1] & ALNFLAG == 0 && o.aln.is_some() {
                crate::bail("this track was not written with --aln");
            }
            if h[1] & DBFLAG != 0 && !o.db {
                crate::bail("this track holds featureIDs: decode it with --db features.tsv --out DIR, or --raw to see what is stored");
            }
            if h[1] & DBFLAG == 0 && o.db {
                crate::bail("this track does not hold featureIDs");
            }
            if h[1] & LENFLAG != 0 && o.aln.is_none() {
                crate::bail("this track leaves out every length the reference implies: decode it with --aln PAF --ref R");
            }
        }
        if let Some(a) = &o.aln {
            if h[1] & ALNFLAG != 0 && (h[13] >> 8) != (a.hash & ((1u64 << 56) - 1)) {
                crate::bail("this PAF, reference or name map is not the one the track was encoded with");
            }
        }
        let was_perm = h[1] & PERMFLAG != 0;
        let perm_tok = h[1] & PERMTOK != 0;
        rd(f.seek(SeekFrom::Start(8 + 128)), "seek");

        let lzl = h[6] as usize;
        let mut lzraw = vec![0u8; crate::ckalloc(lzl)];
        if lzl > 0 { rd(f.read_exact(&mut lzraw), "short labels"); }
        let mut lab: Vec<Vec<u8>> = Vec::new();
        if lzl > 0 {
            let lbuf = zdec(&lzraw, "label table");
            if h[7] > lbuf.len() as u64 { crate::bail("label count exceeds the label table"); }
            if h[7] == 0 { crate::bail("corrupt file: a label table with no labels"); }
            let mut it = lbuf.split(|&c| c == b'\n');
            for _ in 0..h[7] {
                lab.push(it.next().unwrap_or_else(|| crate::bail("bad label table")).to_vec());
            }
            if it.next() != Some(&[][..]) || it.next().is_some() {
                crate::bail("corrupt file: the label table has more lines than its count");
            }
        } else if h[7] != 0 {
            crate::bail("corrupt file: a label count with no label table");
        }

        let names_tok = h[1] & TOKFLAG != 0;
        let names_int = h[1] & INTFLAG != 0;
        if names_tok && names_int {
            crate::bail("corrupt file: names flagged both tokenised and integer");
        }
        let taildec = h[1] & TAILFLAG != 0;
        if taildec && h[13] & 0xff != 64 {
            crate::bail("corrupt file: tail skip is not the format constant");
        }
        let nzl = nml as usize;
        let nz_off = 8 + 128 + lzl as u64;
        let mut nz = vec![0u8; crate::ckalloc(nzl)];
        if nzl > 0 { rd(f.read_exact(&mut nz), "short names"); }
        let mut iz = vec![0u8; crate::ckalloc(h[2] as usize)];
        rd(f.read_exact(&mut iz), "short index");
        let mut sz = vec![0u8; crate::ckalloc(h[3] as usize)];
        rd(f.read_exact(&mut sz), "short sizes");
        let mut itz = vec![0u8; crate::ckalloc(h[11] as usize)];
        rd(f.read_exact(&mut itz), "short index table");

        // A sidecar holds the names --names drop:FILE wrote.
        let mut nzl = nzl;
        if nzl == 0 && was_side {
            let side = o.names_file.clone()
                .unwrap_or_else(|| format!("{path}.names"));
            match std::fs::read(&side) {
                Ok(v) => { nzl = v.len(); nz = v; }
                Err(_) if o.raw => {}
                Err(_) => {
                    eprintln!("cannot open the names sidecar {side}");
                    crate::bail("this track keeps its names in a sidecar and it is missing");
                }
            }
        }
        let mut pz = vec![0u8; crate::ckalloc(if was_perm { h[14] as usize } else { 0 })];
        if !pz.is_empty() { rd(f.read_exact(&mut pz), "short input-order names"); }

        // The metadata hash covers the header with h[15] zeroed, then every
        // stream, and the sidecar's bytes stand in for the names when they
        // are not in the file. A foreign sidecar is caught here.
        if !(was_side && nzl == 0 && o.raw) {
            let mut hz = h;
            hz[15] = 0;
            hz[1] = h1_raw;
            let mut b = Vec::with_capacity(128);
            for v in hz { b.extend_from_slice(&v.to_le_bytes()); }
            let mut mh = kst::fnv(FNV0, &b);
            mh = kst::fnv(mh, &lzraw);
            mh = kst::fnv(mh, &nz);
            mh = kst::fnv(mh, &iz);
            mh = kst::fnv(mh, &sz);
            mh = kst::fnv(mh, &itz);
            mh = kst::fnv(mh, &pz);
            if mh != h[15] {
                if was_side {
                    crate::bail("this names file does not belong to this track: give the file --names drop:FILE wrote, not a FASTA or a list");
                }
                crate::bail("corrupt file: the metadata does not match its hash");
            }
        }
        if was_side && nzl == 0 && o.raw {
            eprintln!("sidecar missing: names are ordinals and the metadata hash is unchecked");
        }
        let _ = nz_off;

        // The block table.
        let was_len = h[1] & LENFLAG != 0;
        let mut blk_nos: Vec<u8> = Vec::new();
        let sbuf = zdec(&sz, "block sizes");
        let mut b = Vec::new();
        let (mut q, mut acc) = (0usize, 0u64);
        while q < sbuf.len() {
            let (a, n1) = vget(&sbuf, q);
            let (b2, n2) = vget(&sbuf, n1);
            let (c2, n3) = vget(&sbuf, n2);
            let (mut d2, n4) = vget(&sbuf, n3);
            let nos = if was_len { let x = (d2 & 3) as u8; d2 >>= 2; x } else { 0 };
            blk_nos.push(nos);
            let (e2, n5) = vget(&sbuf, n4);
            let (_bp, n6) = vget(&sbuf, n5);
            let (g2, n7) = vget(&sbuf, n6);
            q = n7;
            if a > fsz || b2 > fsz || a + b2 > fsz || acc + a + b2 > fsz {
                crate::bail("corrupt index: a block runs past the end of the file");
            }
            if e2 > h[4] { crate::bail("corrupt index: a block claims more rows than the file has"); }
            if e2 != 0 && (c2 == 0 || d2 == 0) {
                crate::bail("corrupt index: a block with rows but an empty stream");
            }
            if c2 > h[4] * 16 + 64 || d2 > h[4] * 16 + 64 {
                crate::bail("corrupt index: a block claims an impossible decoded size");
            }
            b.push(Blk { off: acc, size: a + b2, l: a, raw: c2, lraw: d2, rows: e2, hash: g2 });
            acc += a + b2;
        }
        let tot: u64 = b.iter().map(|x| x.rows).sum();
        let blk = if h[8] != 0 { h[8] } else { 512 };
        for (i, x) in b.iter().enumerate() {
            let last = i + 1 == b.len();
            if (!last && x.rows != blk) || (last && (x.rows == 0 || x.rows > blk)) {
                crate::bail("corrupt index: a block does not hold the block size of rows");
            }
        }
        if tot != h[4] || (b.is_empty() && h[4] != 0) {
            crate::bail("corrupt index: the blocks do not add up to the row count");
        }

        // The index: chunk table, then the concatenated chunk frames.
        let tbuf = zdec(&itz, "index table");
        let mut coff = vec![0u64];
        let mut cpend = vec![];
        let (mut q, mut acc2) = (0usize, 0u64);
        for _ in 0..h[10] {
            let (cl, n1) = vget(&tbuf, q);
            let (sp, n2) = vget(&tbuf, n1);
            q = n2;
            acc2 += cl;
            coff.push(acc2);
            cpend.push(sp);
        }
        cpend.push(h[4]);
        if acc2 != h[2] { crate::bail("index chunk table does not match the index"); }
        if h[10] != h[0].div_ceil(ICHUNK) {
            crate::bail("corrupt file: index chunks do not match the sequence count");
        }
        // Each chunk must hold exactly its share of the sequences and add up
        // to the run count its table entry promised.
        let mut ibx = Vec::new();
        for c in 0..h[10] as usize {
            let d = zdec(&iz[coff[c] as usize..coff[c + 1] as usize], "index chunk");
            let want = if (c as u64) + 1 < h[10] { ICHUNK } else { h[0] - c as u64 * ICHUNK };
            let (mut i, mut sum, mut k) = (0usize, 0u64, 0u64);
            while k < want && i < d.len() {
                let (v, n) = vget(&d, i);
                if v > u32::MAX as u64 { crate::bail("corrupt index: a run count past 32 bits"); }
                sum += v;
                i = n;
                k += 1;
            }
            if k != want || i != d.len() || cpend[c] + sum != cpend[c + 1] {
                crate::bail("corrupt index: an index chunk does not match its table entry");
            }
            ibx.extend_from_slice(&d);
        }

        // A permutation track with no stored blob keeps its names outside,
        // so --names FILE is the source there, not a replacement.
        let perm_ext = was_perm && pz.is_empty();

        // The names are decoded on first use: a query through the index
        // never needs them, and a short-read track has 20 million.
        let mut nf = 0u8;
        for (bit, on) in [(1, was_side), (2, perm_ext), (4, dropped), (8, names_int), (16, names_tok), (32, was_perm), (64, perm_tok)] {
            if on { nf |= bit; }
        }
        let (lenc, idc) = (h[5] as u32 as i32, h[9] as u32 as i32);
        if crate::ffi::nosz(lenc) || crate::ffi::nosz(idc) {
            crate::bail("corrupt file: a NOSZ codec cannot be size-checked");
        }
        if h[5] as u32 >> 31 != 0 || h[9] as u32 >> 31 != 0 {
            crate::bail("corrupt file: codec id out of range");
        }
        let data0 = 8 + 128 + lzl as u64 + nml + h[2] + h[3] + h[11] + pz.len() as u64;
        Reader {
            f, nreads: h[0], rows: h[4], lab, block: blk, b, data0,
            nameix: std::cell::OnceCell::new(),
            idcodec: idc, lencodec: lenc,
            taildec, tailpred: h[12], tailskip: h[13] & 0xff, ibx,
            names: std::cell::OnceCell::new(), nz, pz, names_file: o.names_file.clone(),
            ordnames: o.ordnames, nf, idx,
            hasidx, mhash: h[15], cur_ri: 0, cur_q: 0, cur_pend: 0, cached: None, sc: None,
            feathash: (h[5] >> 32) | ((h[9] >> 32) << 32),
            aln: o.aln.clone(), lenfile: h[1] & LENFLAG != 0 && !o.raw, raw: o.raw,
            alnq: None, blk_nos,
        }
    }

    /// Every row, in order, to a callback. decode and --db share it.
    pub fn for_each_row(&mut self, f: &mut RowFn) {
        let mut c = Cursor { blk: usize::MAX, data: std::rc::Rc::new(Vec::new()),
                             lp: 0, ip: 0, left: 0 };
        let (mut ipos, mut emitted) = (0usize, 0u64);
        let tskip = self.tailskip;
        let (mut name, mut lab) = (Vec::new(), Vec::new());
        for i in 0..self.nreads {
            let (nruns, next) = vget(&self.ibx, ipos);
            ipos = next;
            name.clear(); name.extend_from_slice(self.name_at(i));
            self.alnq = self.aln.as_ref().and_then(|a| a.find(&name));
            let mut st = 0u64;
            for k in 0..nruns {
                let (len, id) = self.step(&mut c, st, k + 1 == nruns && i >= tskip);
                lab.clear();
                self.label(id, &mut lab, st, len);
                f(&name, st, st + len, &lab);
                st += len;
                emitted += 1;
            }
        }
        if emitted != self.rows {
            eprintln!("karyozip: emitted {emitted} rows, the header says {}", self.rows);
            crate::bail("file is truncated");
        }
    }

    /// Every sequence, in order. One forward pass over the blocks.
    pub fn decode(&mut self, out: &mut impl std::io::Write) {
        let mut buf: Vec<u8> = Vec::with_capacity(1 << 22);
        let mut c = Cursor { blk: usize::MAX, data: std::rc::Rc::new(Vec::new()), lp: 0, ip: 0, left: 0 };
        let (mut ipos, mut emitted, mut restored) = (0usize, 0u64, 0u64);
        let tskip = self.tailskip;
        let mut name = Vec::new();
        for i in 0..self.nreads {
            let (nruns, next) = vget(&self.ibx, ipos);
            ipos = next;
            name.clear(); name.extend_from_slice(self.name_at(i));
            self.alnq = self.aln.as_ref().and_then(|a| a.find(&name));
            let mut st = 0u64;
            for k in 0..nruns {
                let (len, f) = self.step(&mut c, st, k + 1 == nruns && i >= tskip);
                let en = st + len;
                buf.extend_from_slice(&name);
                buf.push(b'\t');
                kst::putu64(&mut buf, st);
                buf.push(b'\t');
                kst::putu64(&mut buf, en);
                buf.push(b'\t');
                if self.aln.is_some() && self.lab.get(f as usize).is_some_and(|l| l == b"=") { restored += 1; }
                self.label(f, &mut buf, st, len);
                buf.push(b'\n');
                if buf.len() > (1 << 22) - 4096 {
                    rd(out.write_all(&buf), "write failed");
                    buf.clear();
                }
                st = en;
                emitted += 1;
            }
        }
        rd(out.write_all(&buf), "write failed");
        if self.aln.is_some() { eprintln!("aln: {restored} rows restored from the reference"); }
        if emitted != self.rows {
            eprintln!("karyozip: emitted {emitted} rows, the header says {}", self.rows);
            crate::bail("file is truncated");
        }
    }

    /// One block's two frames back to back, the split, its rows. The last
    /// block stays decoded: without that, a track whose sequences share a
    /// block is decompressed once per sequence.
    fn load(&mut self, b: usize) -> (std::rc::Rc<Vec<u8>>, usize, u64) {
        if let Some((cb, d, sp, r)) = &self.cached {
            if *cb == b { return (d.clone(), *sp, *r); }
        }
        let (d, sp, r) = self.load_raw(b);
        let d = std::rc::Rc::new(d);
        self.cached = Some((b, d.clone(), sp, r));
        (d, sp, r)
    }

    fn load_raw(&mut self, b: usize) -> (Vec<u8>, usize, u64) {
        if b >= self.b.len() { crate::bail("corrupt index: a block number is past the end"); }
        let m = &self.b[b];
        let mut cb = vec![0u8; crate::ckalloc(m.size as usize)];
        rd(self.f.seek(SeekFrom::Start(self.data0 + m.off)), "seek");
        rd(self.f.read_exact(&mut cb), "file is truncated: a block is missing or short");
        if kst::fnv(FNV0, &cb) != m.hash {
            crate::bail("corrupt block: its bytes do not match the stored hash");
        }
        let (la, lb) = (m.l as usize, (m.size - m.l) as usize);
        let lens = frame(&cb[..la], self.lencodec, m.lraw as usize);
        let ids = frame(&cb[la..la + lb], self.idcodec, m.raw as usize);
        let split = lens.len();
        let rows = m.rows;
        let mut data = lens;
        data.extend_from_slice(&ids);
        (data, split, rows)
    }
}

fn frame(src: &[u8], codec: i32, raw: usize) -> Vec<u8> {
    if codec != 0 {
        ffi::paydecomp(src, codec, raw).unwrap_or_else(|| crate::bail("frame decode failed"))
    } else {
        ffi::zstd_decompress(src, raw).unwrap_or_else(|| crate::bail("frame decode failed"))
    }
}

fn zdec(src: &[u8], what: &str) -> Vec<u8> {
    let n = ffi::zstd_frame_size(src)
        .unwrap_or_else(|| crate::bail(&format!("{what} has no usable frame size")));
    ffi::zstd_decompress(src, n)
        .unwrap_or_else(|| crate::bail(&format!("{what} decompress failed")))
}

fn int_names_d(b: &[u8], cnt: u64) -> Vec<u8> {
    let v = zdec(b, "int names");
    let mut o = Vec::with_capacity(crate::ckalloc(cnt as usize * 8));
    let (mut q, mut cur) = (0usize, 0u64);
    for _ in 0..cnt {
        let (u, n) = vget(&v, q);
        q = n;
        cur = cur.wrapping_add(unzz(u));
        if cur > 1 << 50 { crate::bail("corrupt names: an integer name past 2^50"); }
        kst::putu64(&mut o, cur);
        o.push(b'\n');
    }
    if q != v.len() { crate::bail("corrupt names: the integer stream has more than the count"); }
    o
}

pub fn tok_d(b: &[u8]) -> Vec<u8> {
    if b.len() < 8 { crate::bail("truncated tokenised names"); }
    let nchunk = u32::from_le_bytes(b[4..8].try_into().unwrap()) as usize;
    let hdr = 8 + 8 * nchunk;
    if nchunk > (b.len() - 8) / 8 || hdr > b.len() {
        crate::bail("corrupt tokenised names header");
    }
    let mut coff = vec![0usize];
    for (i, _) in (0..nchunk).enumerate() {
        let o = 8 + 4 * (nchunk + i);
        let cl = u32::from_le_bytes(b[o..o + 4].try_into().unwrap());
        coff.push(coff[i] + cl as usize);
    }
    if coff[nchunk] > b.len() - hdr {
        crate::bail("corrupt tokenised names: chunks overrun the stream");
    }
    // One chunk per thread, as the C's tokd does. The tokeniser separates
    // with NUL; the format stores newlines.
    use rayon::prelude::*;
    let parts: Vec<Vec<u8>> = (0..nchunk).into_par_iter().map(|i| {
        let mut d = ffi::tok_decode(&b[hdr + coff[i]..hdr + coff[i + 1]])
            .unwrap_or_else(|| crate::bail("name tokeniser decode failed"));
        for c in d.iter_mut() { if *c == 0 { *c = b'\n'; } }
        d
    }).collect();
    parts.concat()
}

pub use crate::kst::{unzz, vget};

/* ------------------------------------------------ one sequence at a time */

/// A position in the row stream: the block, its bytes, a cursor per frame,
/// and how many of its rows are left.
struct Cursor { blk: usize, data: std::rc::Rc<Vec<u8>>, lp: usize, ip: usize, left: u64 }

impl Reader {
    /// (block, offset in it, run count) for sequence `ri`. A forward cursor,
    /// so walking every sequence is linear and not quadratic.
    fn locate(&mut self, ri: u64) -> (usize, u64, u64) {
        if ri < self.cur_ri { self.cur_ri = 0; self.cur_q = 0; self.cur_pend = 0; }
        let mut nruns = 0u64;
        while self.cur_ri <= ri {
            let (n, next) = vget(&self.ibx, self.cur_q);
            self.cur_q = next;
            if self.cur_pend / self.block > u32::MAX as u64 {
                crate::bail("corrupt index: a block number past 32 bits");
            }
            nruns = n;
            self.cur_pend += n;
            self.cur_ri += 1;
        }
        let start = self.cur_pend - nruns;
        ((start / self.block) as usize, start % self.block, nruns)
    }

    fn at(&mut self, blk: usize, off: u64) -> Cursor {
        let (data, split, rows) = self.load(blk);
        let mut c = Cursor { blk, data, lp: 0, ip: split, left: rows };
        for _ in 0..off { c.lp = skip(&c.data, c.lp); c.ip = skip(&c.data, c.ip); }
        c.left -= off;
        c
    }

    /// Next row: (length, label id). Steps to the next block when this one
    /// ends. Under --aln a zero length means "ask the reference".
    fn step(&mut self, c: &mut Cursor, st: u64, last: bool) -> (u64, u64) {
        if c.left == 0 {
            let b = if c.blk == usize::MAX { 0 } else { c.blk + 1 };
            let (data, split, rows) = self.load(b);
            *c = Cursor { blk: b, data, lp: 0, ip: split, left: rows };
        }
        let (v, l2) = vget(&c.data, c.lp);
        c.lp = l2;
        let len = if last && self.taildec {
            let t = self.tailpred.wrapping_add(unzz(v));
            if t <= st { crate::bail("corrupt track: sequence total is not past its last run"); }
            t - st
        } else if v != 0 {
            v
        } else if self.raw {
            0
        } else {
            let lendec = self.lenfile && self.blk_nos.get(c.blk).copied().unwrap_or(0) != 1;
            if !lendec { crate::bail("corrupt track: a zero length outside a sentinel block"); }
            let q = self.alnq.unwrap_or_else(|| crate::bail("a length sentinel has no reference answer: wrong PAF or --ref"));
            let (e, _) = self.aln.as_ref().unwrap().end(q, st);
            if e <= st { crate::bail("a length sentinel has no reference answer: wrong PAF or --ref"); }
            e - st
        };
        let (f, i2) = vget(&c.data, c.ip);
        c.ip = i2;
        c.left -= 1;
        (len, f)
    }

    fn label(&self, f: u64, out: &mut Vec<u8>, st: u64, len: u64) {
        let l: &[u8] = match self.lab.get(f as usize) {
            Some(l) => l,
            None if self.lab.is_empty() => { kst::putu64(out, f); return; }
            None => crate::bail("corrupt block: a label id past the table"),
        };
        if let (b"=", Some(a)) = (l, self.aln.as_ref()) {
            let q = self.alnq.unwrap_or_else(|| crate::bail("a delta row has no reference answer: wrong PAF or reference track"));
            let rl = a.row(q, st, len);
            if rl == u32::MAX { crate::bail("a delta row has no reference answer: wrong PAF or reference track"); }
            out.extend_from_slice(&self.lab[rl as usize]);
            return;
        }
        out.extend_from_slice(l);
    }

    /// What the C's ensurenames does: the names, by whichever codec stored
    /// them, or from --names FILE, or the ordinals.
    fn decode_names(&self) -> (Vec<u8>, Vec<u32>) {
        let (side, perm_ext, dropped) = (self.nf & 1 != 0, self.nf & 2 != 0, self.nf & 4 != 0);
        let (names_int, names_tok, was_perm, perm_tok) = (self.nf & 8 != 0, self.nf & 16 != 0, self.nf & 32 != 0, self.nf & 64 != 0);
        let (nz, pz, n) = (&self.nz, &self.pz, self.nreads);
        let mut names = if let (Some(f), false) = (&self.names_file, side || perm_ext) {
            if !dropped && !nz.is_empty() {
                crate::bail("this track stores its own names: --names FILE would replace them, so it is refused. Unset it to read the names in the file");
            }
            read_name_file(f, n)
        } else if dropped && nz.is_empty() {
            if !self.raw && !self.ordnames {
                crate::bail("the names of this track were dropped at encode time: give them back with --names FILE, or pass --names ordinal to accept the ordinals that are stored instead");
            }
            ordinals(n)
        } else if nz.is_empty() {
            ordinals(n)
        } else if names_int {
            int_names_d(nz, n)
        } else if names_tok {
            tok_d(nz)
        } else {
            zdec(nz, "names")
        };
        if names.last() == Some(&b'\n') { names.pop(); }
        let mut noff: Vec<u32> = vec![0];
        if !names.is_empty() {
            for (i, &c) in names.iter().enumerate() { if c == b'\n' { noff.push(i as u32 + 1); } }
            noff.push(names.len() as u32 + 1);
        }
        if noff.len() as u64 - 1 != n {
            crate::bail("names stream does not match the track's sequence count");
        }
        // --sort stores a rank per record; the input-order names sit beside
        // them, or come from --names FILE when the blob is empty.
        if (was_perm && !pz.is_empty()) || perm_ext {
            let resolve = if perm_ext {
                match &self.names_file {
                    Some(f) => Some(read_name_file(f, n)),
                    None if self.raw || self.ordnames => None,
                    None => crate::bail("--sort --no-names keeps the rank of every record and no name text: give the original FASTA/FASTQ in --names FILE, or --ordinal-names to see the ranks"),
                }
            } else if perm_tok {
                Some(tok_d(pz))
            } else {
                Some(zdec(pz, "input-order names"))
            };
            if let Some(nat) = resolve {
                let tab: Vec<&[u8]> = nat.split(|&c| c == b'\n').filter(|l| !l.is_empty()).collect();
                if tab.len() as u64 != n {
                    crate::bail("input-order names do not match the track's sequence count");
                }
                let (mut nn, mut no) = (Vec::new(), vec![0u32]);
                for i in 0..noff.len() - 1 {
                    let nm = &names[noff[i] as usize..noff[i + 1] as usize - 1];
                    let v: u64 = std::str::from_utf8(nm).ok().and_then(|t| t.parse().ok())
                        .unwrap_or_else(|| crate::bail("a stored rank is not a record number"));
                    if v >= n { crate::bail("a stored rank is not a record number"); }
                    nn.extend_from_slice(tab[v as usize]); nn.push(b'\n'); no.push(nn.len() as u32);
                }
                names = nn; noff = no;
            }
        }
        (names, noff)
    }

    /// The C's namesrc: an index is bound to the --names FILE it was built
    /// under, by path and stat; a pipe binds nothing.
    fn namesrc(&self) -> u64 {
        use std::os::linux::fs::MetadataExt;
        let Some(f) = self.names_file.as_ref().filter(|_| self.nf & 1 == 0) else { return 0 };
        let m = rd(std::fs::metadata(f), "cannot stat the file given to --names");
        if !m.is_file() { return 0; }
        let mut h = kst::fnv(FNV0, f.as_bytes());
        for v in [m.st_size() as u64, m.st_ino(), m.st_dev(), m.st_mtime() as u64, m.st_mtime_nsec() as u64,
                  m.st_ctime() as u64, m.st_ctime_nsec() as u64] {
            h = kst::fnv(h, &v.to_le_bytes());
        }
        h
    }

    /// The appended index, as the C reads it: a binary search over the
    /// sorted hash column, then the second hash and the ordinal confirm.
    fn idxlookup(&self, nm: &[u8]) -> i64 {
        use std::os::unix::fs::FileExt;
        let Some((base, bind)) = self.idx else { return -1 };
        if bind != self.mhash ^ self.namesrc() { return -1; }
        let n = self.nreads;
        let word = |off: u64| -> u64 { let mut b = [0u8; 8]; rd(self.f.read_exact_at(&mut b, off), "index read"); u64::from_le_bytes(b) };
        let ord = |i: u64| -> u64 { let mut b = [0u8; 4]; rd(self.f.read_exact_at(&mut b, base + n * 16 + i * 4), "index read"); u32::from_le_bytes(b) as u64 };
        let (want, want2) = (kst::fnv(FNV0, nm), h64b(nm));
        let (mut lo, mut hi) = (0u64, n);
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let k = word(base + mid * 8);
            if k == want {
                let (mut a2, mut b2) = (mid, mid + 1);
                while a2 > 0 && word(base + (a2 - 1) * 8) == want { a2 -= 1; }
                while b2 < n && word(base + b2 * 8) == want { b2 += 1; }
                let fits = |j: u64| { let o = ord(j); o < n && word(base + n * 8 + j * 8) == want2 ^ o.wrapping_mul(H64B_K) };
                if b2 - a2 == 1 {
                    let o = ord(a2);
                    let named = self.names.get().is_some_and(|_| self.name_at(o) != nm);
                    return if fits(a2) && !named { o as i64 } else { -1 };
                }
                for j in a2..b2 {
                    if fits(j) && self.name_at(ord(j)) == nm { return ord(j) as i64; }
                }
                return -1;
            }
            if k < want { lo = mid + 1; } else { hi = mid; }
        }
        -1
    }

    fn find(&self, seq: &[u8]) -> i64 {
        let hit = self.idxlookup(seq);
        if hit >= 0 { return hit; }
        let ix = self.nameix.get_or_init(|| {
            let mut m = HashMap::with_capacity(self.nreads as usize);
            for i in 0..self.nreads {
                if m.insert(self.name_at(i).to_vec(), i).is_some() {
                    crate::bail("a sequence name repeats: the names are not usable");
                }
            }
            m
        });
        ix.get(seq).map_or(-1, |&i| i as i64)
    }

    /// region F SEQ | SEQ POS | SEQ START END [HALO].
    /// --db: the sequence's rows over the picked sets, touching equal tuples
    /// merged, so one set gives what its own track gives. Needs the whole sequence.
    fn region_set(&mut self, ri: u64, seq: &[u8], lo: u64, hi: u64, halo: u64,
                  fe: &crate::db::Feats, cols: &[usize], out: &mut impl std::io::Write) {
        let (blk, off, nruns) = self.locate(ri);
        let tskip = self.tailskip;
        let mut c = self.at(blk, off);
        let (mut st, mut en, mut lab): (Vec<u64>, Vec<u64>, Vec<Vec<u8>>) = (Vec::new(), Vec::new(), Vec::new());
        let (mut s, mut l, mut t) = (0u64, Vec::new(), Vec::new());
        for i in 0..nruns {
            let (len, f) = self.step(&mut c, s, i + 1 == nruns && ri >= tskip);
            l.clear(); self.label(f, &mut l, s, len);
            fe.joined(fe.featid(&l), cols, &mut t);
            if en.last() == Some(&s) && lab.last() == Some(&t) { *en.last_mut().unwrap() = s + len; }
            else { st.push(s); en.push(s + len); lab.push(t.clone()); }
            s += len;
        }
        let n = st.len();
        let (mut a, mut b) = (0, n);
        if hi > 0 {
            let (mut f, mut l) = (n, 0usize);
            for i in 0..n { if en[i] > lo && st[i] < hi { if i < f { f = i; } l = i; } }
            if f > l { return; }
            a = f.saturating_sub(halo as usize);
            b = if n - l - 1 <= halo as usize { n } else { l + halo as usize + 1 };
        }
        let mut buf = Vec::with_capacity(1 << 20);
        for i in a..b {
            buf.extend_from_slice(seq); buf.push(b'\t');
            kst::putu64(&mut buf, st[i]); buf.push(b'\t');
            kst::putu64(&mut buf, en[i]); buf.push(b'\t');
            buf.extend_from_slice(&lab[i]); buf.push(b'\n');
            if buf.len() > (1 << 20) - 4096 { rd(out.write_all(&buf), "write failed"); buf.clear(); }
        }
        rd(out.write_all(&buf), "write failed");
    }

    pub fn region(&mut self, seq: &[u8], lo: u64, hi: u64, halo: u64,
                  feats: Option<&crate::db::Feats>, cols: &[usize], out: &mut impl std::io::Write) {
        let ri = self.find(seq);
        if ri < 0 { crate::bail("sequence not found"); }
        let ri = ri as u64;
        self.alnq = self.aln.as_ref().and_then(|a| a.find(seq));
        if let Some(fe) = feats { return self.region_set(ri, seq, lo, hi, halo, fe, cols, out); }
        let (blk, off, nruns) = self.locate(ri);
        let tskip = self.tailskip;
        let last_of = |k: u64| k + 1 == nruns && ri >= tskip;

        // With a halo the runs to print are found first, then reprinted with
        // the neighbours the halo asks for.
        let (mut kfirst, mut klast, mut lo, mut hi) = (0u64, nruns, lo, hi);
        if halo > 0 && (lo > 0 || hi > 0) {
            let mut c = self.at(blk, off);
            let (mut s2, mut first, mut last) = (0u64, nruns, 0u64);
            for k in 0..nruns {
                if hi > 0 && s2 >= hi { break; }
                let (len, _) = self.step(&mut c, s2, last_of(k));
                let e2 = s2 + len;
                if e2 > lo && (hi == 0 || s2 < hi) {
                    if k < first { first = k; }
                    last = k;
                }
                s2 = e2;
            }
            if first > last { return; }
            kfirst = first.saturating_sub(halo);
            klast = if nruns - last - 1 <= halo { nruns } else { last + halo + 1 };
            lo = 0; hi = 0;
        }

        let mut c = self.at(blk, off);
        let mut buf = Vec::with_capacity(1 << 20);
        let mut st = 0u64;
        for k in 0..klast {
            let (len, f) = self.step(&mut c, st, last_of(k));
            let en = st + len;
            if k >= kfirst && en > lo && (hi == 0 || st < hi) {
                buf.extend_from_slice(seq);
                buf.push(b'\t');
                kst::putu64(&mut buf, st);
                buf.push(b'\t');
                kst::putu64(&mut buf, en);
                buf.push(b'\t');
                self.label(f, &mut buf, st, len);
                buf.push(b'\n');
                if buf.len() > (1 << 20) - 4096 {
                    rd(out.write_all(&buf), "write failed");
                    buf.clear();
                }
            }
            st = en;
            if hi > 0 && st >= hi { break; }
        }
        rd(out.write_all(&buf), "write failed");
    }
}

fn skip(b: &[u8], mut i: usize) -> usize {
    while i < b.len() && b[i] >= 0x80 { i += 1; }
    i + 1
}

/* ------------------------------------------------------------ the index */

/// "KSTIDX01" reversed into the trailer's last word.
const H64B_K: u64 = 0x9e37_79b9_7f4a_7c15;

fn h64b(s: &[u8]) -> u64 {
    let mut h = H64B_K;
    for &c in s {
        h ^= c as u64;
        h = h.wrapping_mul(0xff51_afd7_ed55_8ccd);
        h ^= h >> 32;
    }
    h
}

impl Reader {
    /// Appends a name index to the track, magic last so a torn write is
    /// simply not seen. One file, not two.
    pub fn build_index(&mut self, path: &str) {
        if self.hasidx { eprintln!("karyozip: {path} already carries an index"); return; }
        if self.nreads > u32::MAX as u64 {
            crate::bail("more than 2^32 sequences: the index stores the ordinal in 32 bits");
        }
        let n = self.nreads as usize;
        let mut pair: Vec<(u64, u32)> = (0..n)
            .map(|i| (kst::fnv(FNV0, self.name_at(i as u64)), i as u32)).collect();
        pair.sort_unstable();
        let dup = pair.windows(2).filter(|w| w[0].0 == w[1].0).count() as u64;
        let mut t = Vec::with_capacity(crate::ckalloc(n * 20 + 32));
        for (k, _) in &pair { t.extend_from_slice(&k.to_le_bytes()); }
        for (_, o) in &pair {
            t.extend_from_slice(&(h64b(self.name_at(*o as u64))
                                  ^ (*o as u64).wrapping_mul(H64B_K)).to_le_bytes());
        }
        for (_, o) in &pair { t.extend_from_slice(&o.to_le_bytes()); }
        for v in [self.nreads, dup, self.mhash ^ self.namesrc(), IDXMAGIC] {
            t.extend_from_slice(&v.to_le_bytes());
        }
        use std::io::Write;
        let mut f = rd(std::fs::OpenOptions::new().append(true).open(path),
                       "cannot append the index to the track");
        rd(f.write_all(&t), "index write failed");
        rd(f.flush(), "index write failed");
        eprintln!("index: {} sequences, {} B appended, {dup} hash collisions",
                  self.nreads, 32 + n * 20);
    }
}

impl Reader {
    pub fn nlabels(&self) -> usize { self.lab.len() }
    pub fn has_aln(&self) -> bool { self.aln.is_some() }
    pub fn label_at(&self, i: usize) -> &[u8] { &self.lab[i] }
    pub fn name_at(&self, i: u64) -> &[u8] {
        let (names, noff) = &**self.names.get_or_init(|| std::sync::Arc::new(self.decode_names()));
        let i = i as usize; &names[noff[i] as usize..noff[i + 1] as usize - 1]
    }
    /// The decoded names, to hand to a worker's reader of the same track so
    /// the table is decoded once, as the C's smooth_mt borrows it.
    pub fn names_shared(&self) -> std::sync::Arc<(Vec<u8>, Vec<u32>)> {
        self.names.get_or_init(|| std::sync::Arc::new(self.decode_names())).clone()
    }
    pub fn adopt_names(&mut self, n: std::sync::Arc<(Vec<u8>, Vec<u32>)>) { let _ = self.names.set(n); }
    pub fn find_seq(&self, s: &[u8]) -> i64 { self.find(s) }

    /// One sequence as arrays: starts, ends, hierarchy node per run, and
    /// whether the label was "novel". What smoothing works on.
    pub fn seq_rows(&mut self, ri: u64, labmap: &[i32], isnovel: &[u8], root: i32)
        -> (Vec<u64>, Vec<u64>, Vec<i32>, Vec<u8>) {
        let (blk, off, nruns) = self.locate(ri);
        self.alnq = self.aln.as_ref().and_then(|a| a.find(self.name_at(ri)));
        let mut c = match self.sc.take() {
            Some((b, o, mut c)) if b == blk && o <= off => {
                for _ in o..off { c.lp = skip(&c.data, c.lp); c.ip = skip(&c.data, c.ip); }
                c.left -= off - o; c
            }
            _ => self.at(blk, off),
        };
        let n = nruns as usize;
        let (mut st, mut en) = (Vec::with_capacity(n), Vec::with_capacity(n));
        let (mut feat, mut nov) = (Vec::with_capacity(n), Vec::with_capacity(n));
        let tskip = self.tailskip;
        let mut s = 0u64;
        for k in 0..nruns {
            let (len, f) = self.step(&mut c, s, k + 1 == nruns && ri >= tskip);
            let mut f = f as usize;
            if f >= self.lab.len() && !self.lab.is_empty() {
                crate::bail("corrupt block: a label id past the table");
            }
            let mut m = labmap.get(f).copied().unwrap_or(root);
            if m == i32::MIN {              // stored as =: the reference row's label
                let miss = || crate::bail("a delta row has no reference answer: wrong PAF or reference track");
                let rl = self.aln.as_ref().unwrap().row(self.alnq.unwrap_or_else(miss), s, len);
                if rl == u32::MAX || rl as usize >= self.lab.len() { miss(); }
                f = rl as usize; m = labmap[f];
            }
            if m == i32::MIN + 1 { crate::bail("corrupt --db track: a label is not a featureID"); }
            let isnov = isnovel.get(f).copied().unwrap_or(0);
            st.push(s);
            en.push(s + len);
            nov.push(isnov);
            feat.push(if isnov != 0 { root } else { m });
            s += len;
        }
        let rows = self.b[c.blk].rows;
        self.sc = Some((c.blk, rows - c.left, c));
        (st, en, feat, nov)
    }
}
