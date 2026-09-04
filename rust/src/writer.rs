//! The KST writer. Rows arrive grouped by sequence, gapless, each sequence
//! starting at 0; blocks fill and are staged; the file is assembled at close.
//! Every layout decision here is c/karyozip.c's, because the two must produce
//! the same bytes.

use crate::ffi;
use crate::kst::{self, FNV0};
use rayon::prelude::*;
use std::collections::HashMap;
use std::collections::HashSet;
use std::io::{Seek, SeekFrom, Write};

pub const MAGIC: &[u8; 8] = b"KST00001";
const TOKFLAG: u64 = 1 << 63;
const INTFLAG: u64 = 1 << 62;
const TAILFLAG: u64 = 1 << 55;
const ICHUNK: u64 = 65536;
const TAILSKIP: usize = 64;

pub struct Options {
    pub block: u64,
    pub level: i32,
    pub keep_names: bool,
    pub id_codec: Option<i32>,
    pub len_codec: Option<i32>,
}

impl Default for Options {
    fn default() -> Self {
        Options { block: 262144, level: 17, keep_names: false,
                  id_codec: None, len_codec: None }
    }
}

/// One staged block, while its payload waits in the staging file.
struct BlkMeta { off: u64, len: usize, lens: usize, rows: u64, bp: u64 }

pub struct Writer {
    opt: Options,
    lens: Vec<u8>,
    ids: Vec<u8>,
    names: Vec<u8>,
    ibx: Vec<u8>,
    nidx: u64,
    ordok: bool,
    labtab: Vec<Vec<u8>>,
    labmap: HashMap<Vec<u8>, u32>,
    lastid: Option<u32>,
    cur: Vec<u8>,
    currun: u64,
    curbp: u64,
    blkbp0: u64,
    pend: u64,
    nblk: u64,
    rows: u64,
    plen: u64,
    phas: bool,
    tailsamp: Vec<u64>,
    tailpred: u64,
    tailon: Option<bool>,
    raw: std::fs::File,
    tmpdir: std::path::PathBuf,
    rawlen: u64,
    bm: Vec<BlkMeta>,
}

fn tmpfile(dir: &std::path::Path) -> std::io::Result<std::fs::File> {
    // An unlinked file: it disappears when the handle closes, however we exit.
    let p = dir.join(format!("karyozip.{}.{:x}", std::process::id(),
                             std::time::SystemTime::now()
                                 .duration_since(std::time::UNIX_EPOCH).unwrap().subsec_nanos()));
    let f = std::fs::OpenOptions::new().read(true).write(true).create_new(true).open(&p)?;
    std::fs::remove_file(&p)?;
    Ok(f)
}

impl Writer {
    pub fn new(opt: Options, tmpdir: &std::path::Path) -> std::io::Result<Self> {
        Ok(Writer {
            opt, lens: Vec::new(), ids: Vec::new(), names: Vec::new(), ibx: Vec::new(),
            nidx: 0, ordok: true, labtab: Vec::new(), labmap: HashMap::new(), lastid: None,
            cur: Vec::new(), currun: 0, curbp: 0, blkbp0: 0, pend: 0, nblk: 0, rows: 0,
            plen: 0, phas: false, tailsamp: Vec::new(), tailpred: 0, tailon: None,
            raw: tmpfile(tmpdir)?, tmpdir: tmpdir.to_path_buf(), rawlen: 0, bm: Vec::new(),
        })
    }

    fn intern(&mut self, label: &[u8]) -> u32 {
        if let Some(id) = self.lastid {
            if self.labtab[id as usize] == label { return id; }
        }
        if let Some(&id) = self.labmap.get(label) {
            self.lastid = Some(id);
            return id;
        }
        let id = self.labtab.len() as u32;
        self.labtab.push(label.to_vec());
        self.labmap.insert(label.to_vec(), id);
        self.lastid = Some(id);
        id
    }

    /// The held row was the last of its sequence: its length is either stored
    /// plainly or, once the predictor is on, as a delta from the median.
    fn tail(&mut self) {
        if !self.phas { return; }
        if self.tailsamp.len() < TAILSKIP {
            self.tailsamp.push(self.curbp);
            if self.tailsamp.len() == TAILSKIP {
                let mut s = self.tailsamp.clone();
                s.sort_unstable();
                self.tailpred = s[TAILSKIP / 2];
                let dev: u64 = self.tailsamp.iter()
                    .map(|&v| if v > self.tailpred { v - self.tailpred } else { self.tailpred - v })
                    .sum();
                self.tailon = Some(dev * 12 <= self.tailpred * TAILSKIP as u64);
            }
            kst::put(&mut self.lens, self.plen);
            self.phas = false;
            return;
        }
        if self.tailon != Some(true) {
            kst::put(&mut self.lens, self.plen);
            self.phas = false;
            return;
        }
        kst::put(&mut self.lens, kst::zz(self.curbp as i64 - self.tailpred as i64));
        self.phas = false;
    }

    fn pushread(&mut self) {
        if self.cur.is_empty() { return; }
        kst::put(&mut self.ibx, self.currun);
        if self.ordok {
            let want = self.nidx.to_string();
            if self.cur == want.as_bytes() {
                self.nidx += 1;
                return;
            }
            // Not 0,1,2,... after all: write out the ordinals so far.
            for j in 0..self.nidx {
                self.names.extend_from_slice(j.to_string().as_bytes());
                self.names.push(b'\n');
            }
            self.ordok = false;
        }
        self.nidx += 1;
        self.names.extend_from_slice(&self.cur);
        self.names.push(b'\n');
    }

    fn flush(&mut self) -> std::io::Result<()> {
        if self.pend == 0 { return Ok(()); }
        let off = self.rawlen;
        self.raw.seek(SeekFrom::Start(off))?;
        self.raw.write_all(&self.lens)?;
        self.raw.write_all(&self.ids)?;
        let len = self.lens.len() + self.ids.len();
        self.bm.push(BlkMeta { off, len, lens: self.lens.len(), rows: self.pend, bp: self.blkbp0 });
        self.rawlen += len as u64;
        self.nblk += 1;
        self.lens.clear();
        self.ids.clear();
        self.pend = 0;
        Ok(())
    }

    pub fn row(&mut self, name: &[u8], run_len: u64, label: &[u8], start: u64)
        -> Result<(), String> {
        if name.is_empty() { return Err("empty sequence name".into()); }
        if run_len == 0 {
            return Err("end equals start: a run must cover at least one base".into());
        }
        if start > u64::MAX - run_len {
            return Err("a coordinate does not fit in 64 bits".into());
        }
        if name != self.cur.as_slice() {
            if start != 0 { return Err("a sequence does not start at 0".into()); }
            if name.len() >= 512 { return Err("name too long".into()); }
            self.tail();
            self.pushread();
            self.cur.clear();
            self.cur.extend_from_slice(name);
            self.currun = 0;
            self.curbp = 0;
            if self.pend > 0 && self.pend >= self.opt.block { self.flush().map_err(io)?; }
        }
        if start != self.curbp {
            return Err("gap or overlap: start is not the previous end".into());
        }
        let id = self.intern(label) as u64;
        if self.phas {
            kst::put(&mut self.lens, self.plen);
            self.phas = false;
        }
        if self.pend >= self.opt.block { self.flush().map_err(io)?; }
        if self.pend == 0 { self.blkbp0 = self.curbp; }
        self.plen = run_len;
        self.phas = true;
        kst::put(&mut self.ids, id);
        self.pend += 1;
        self.currun += 1;
        self.rows += 1;
        self.curbp += run_len;
        Ok(())
    }

    pub fn close(mut self, out: &mut impl Write) -> Result<Stats, String> {
        self.tail();
        self.pushread();
        self.flush().map_err(io)?;
        let level = self.opt.level;

        // The codec search samples the first, middle and last block, which is
        // why it cannot run until every row has arrived.
        let mut picks: Vec<usize> = Vec::new();
        if !self.bm.is_empty() {
            picks.push(0);
            if self.bm.len() > 2 { picks.push(self.bm.len() / 2); }
            if self.bm.len() > 1 { picks.push(self.bm.len() - 1); }
        }
        let rawtot: usize = self.bm.iter().map(|b| b.len).sum();
        let seekable = self.bm.len() > 1 || rawtot as u64 > self.opt.block;
        let sampled: Vec<(Vec<u8>, Vec<u8>)> = picks.iter()
            .map(|&i| self.read_block(i))
            .collect::<std::io::Result<_>>().map_err(io)?;
        let sref: Vec<(&[u8], &[u8])> = sampled.iter()
            .map(|(l, i)| (l.as_slice(), i.as_slice())).collect();
        let p = pick(&sref, seekable, level, self.opt.id_codec, self.opt.len_codec);

        // Compress apart, write in order, a batch at a time. Holding every
        // compressed block would cost the whole file in memory; the C stages
        // them in a temp file for the same reason.
        let mut sizes = Vec::new();
        let mut body = tmpfile(&self.tmpdir).map_err(io)?;
        let (mut byt_lens, mut byt_labs) = (0u64, 0u64);
        let batch = rayon::current_num_threads() * 4;
        for chunk in (0..self.bm.len()).collect::<Vec<_>>().chunks(batch) {
            let done: Vec<(Vec<u8>, Vec<u8>, usize, usize, u64, u64)> = chunk
                .par_iter()
                .map(|&i| {
                    let (lens, ids) = self.read_block(i).expect("staged block");
                    let a = comp_lens(&lens, p.len, level).expect("lens frame");
                    let b = comp_ids(&ids, p.id, level).expect("ids frame");
                    (a, b, lens.len(), ids.len(), self.bm[i].rows, self.bm[i].bp)
                })
                .collect();
            for (a, b, nl, ni, rows, bp) in &done {
                body.write_all(a).map_err(io)?;
                body.write_all(b).map_err(io)?;
                byt_lens += a.len() as u64;
                byt_labs += b.len() as u64;
                kst::put(&mut sizes, a.len() as u64);
                kst::put(&mut sizes, b.len() as u64);
                kst::put(&mut sizes, *ni as u64);
                kst::put(&mut sizes, *nl as u64);
                kst::put(&mut sizes, *rows);
                kst::put(&mut sizes, *bp);
                kst::put(&mut sizes, kst::fnv(kst::fnv(FNV0, a), b));
            }
        }

        let mut labs = Vec::new();
        for l in &self.labtab { labs.extend_from_slice(l); labs.push(b'\n'); }

        let mut ordinal = self.ordok && self.nidx > 0;
        if self.opt.keep_names { ordinal = false; }
        if self.nidx > 1 { check_unique(&self.names)?; }
        if !ordinal && self.ordok && self.names.is_empty() && self.nidx > 0 {
            for i in 0..self.nidx {
                self.names.extend_from_slice(i.to_string().as_bytes());
                self.names.push(b'\n');
            }
        }
        if ordinal { self.names.clear(); }

        let (nz, tokused, intused) = if ordinal {
            (Vec::new(), false, false)
        } else {
            let b = best_names(&self.names, self.nidx, level);
            (b.bytes, b.kind == 2, b.kind == 1)
        };

        // The index is chunked so a lookup decodes 65536 entries, not all.
        let meta = if level < 19 { level } else { 19 };
        let (mut ichunks, mut itab, mut nichunk) = (Vec::new(), Vec::new(), 0u64);
        {
            let mut q = 0usize;
            let (mut cpend, mut seen, mut spend) = (0u64, 0u64, 0u64);
            let mut cstart = 0usize;
            while q < self.ibx.len() {
                let (nruns, next) = vget(&self.ibx, q);
                q = next;
                cpend += nruns;
                seen += 1;
                if seen % ICHUNK == 0 || q >= self.ibx.len() {
                    let cz = ffi::zstd(&self.ibx[cstart..q], meta);
                    ichunks.extend_from_slice(&cz);
                    kst::put(&mut itab, cz.len() as u64);
                    kst::put(&mut itab, spend);
                    nichunk += 1;
                    cstart = q;
                    spend = cpend;
                }
            }
        }
        let itz = ffi::zstd(&itab, meta);
        let sz = ffi::zstd(&sizes, meta);
        let lz = if labs.is_empty() { Vec::new() } else { ffi::zstd(&labs, level) };

        let mut h = [0u64; 16];
        h[0] = self.nidx;
        h[1] = nz.len() as u64
             | if tokused { TOKFLAG } else { 0 }
             | if intused { INTFLAG } else { 0 }
             | if self.tailon == Some(true) { TAILFLAG } else { 0 };
        h[2] = ichunks.len() as u64;
        h[3] = sz.len() as u64;
        h[4] = self.rows;
        h[5] = p.len as u32 as u64;
        h[6] = lz.len() as u64;
        h[7] = self.labtab.len() as u64;
        h[8] = self.opt.block;
        h[9] = p.id as u32 as u64;
        h[10] = nichunk;
        h[11] = itz.len() as u64;
        h[12] = self.tailpred;
        h[13] = TAILSKIP as u64;
        h[14] = 0;
        let mut hb = Vec::with_capacity(128);
        for v in h { hb.extend_from_slice(&v.to_le_bytes()); }
        let mut mh = kst::fnv(FNV0, &hb);
        mh = kst::fnv(mh, &lz);
        mh = kst::fnv(mh, &nz);
        mh = kst::fnv(mh, &ichunks);
        mh = kst::fnv(mh, &sz);
        mh = kst::fnv(mh, &itz);
        h[15] = mh;
        hb.truncate(120);
        hb.extend_from_slice(&h[15].to_le_bytes());

        out.write_all(MAGIC).map_err(io)?;
        out.write_all(&hb).map_err(io)?;
        out.write_all(&lz).map_err(io)?;
        out.write_all(&nz).map_err(io)?;
        out.write_all(&ichunks).map_err(io)?;
        out.write_all(&sz).map_err(io)?;
        out.write_all(&itz).map_err(io)?;
        body.seek(SeekFrom::Start(0)).map_err(io)?;
        std::io::copy(&mut body, out).map_err(io)?;
        out.flush().map_err(io)?;

        Ok(Stats {
            rows: self.rows, reads: self.nidx, blocks: self.nblk,
            names: nz.len(), index: ichunks.len(), offs: sz.len(),
            lens: byt_lens, labs: byt_labs, raw: rawtot,
            tokused, intused, ordinal,
            id_codec: p.id, len_codec: p.len, seekable,
        })
    }

    fn read_block(&self, i: usize) -> std::io::Result<(Vec<u8>, Vec<u8>)> {
        use std::os::unix::fs::FileExt;
        let m = &self.bm[i];
        let mut b = vec![0u8; m.len];
        self.raw.read_exact_at(&mut b, m.off)?;
        let ids = b.split_off(m.lens);
        Ok((b, ids))
    }
}

pub struct Stats {
    pub rows: u64, pub reads: u64, pub blocks: u64,
    pub names: usize, pub index: usize, pub offs: usize,
    pub lens: u64, pub labs: u64, pub raw: usize,
    pub tokused: bool, pub intused: bool, pub ordinal: bool,
    pub id_codec: i32, pub len_codec: i32, pub seekable: bool,
}

fn vget(b: &[u8], mut i: usize) -> (u64, usize) {
    let (mut v, mut s) = (0u64, 0u32);
    loop {
        let c = b[i];
        i += 1;
        v |= ((c & 0x7f) as u64) << s;
        s += 7;
        if c < 0x80 { return (v, i); }
    }
}

fn io(e: std::io::Error) -> String { e.to_string() }


/* ------------------------------------------------- choosing a codec */
/// The eight candidates the default search tries, beside zstd. Codec 0 means
/// zstd, so bare rANS order 0 is spelled 0x100.
pub const CAND_FAST: [i32; 8] = [
    0x100, 1, 65, 521,
    ffi::CODEC_ARITH, ffi::CODEC_ARITH | 1,
    ffi::CODEC_BZIP3 | 16, ffi::CODEC_BZIP3 | 64,
];

/// bzip3 decodes at about 8 MB/s against zstd's ~1000, so it is the wrong
/// pick for a file a query will seek into.
fn fast_enough(codec: i32) -> bool {
    codec & ffi::CODEC_BZIP3 == 0
}

pub struct Pick {
    pub id: i32,
    pub len: i32,
}

/// `samples` are whole blocks: (lens frame, ids frame).
pub fn pick(samples: &[(&[u8], &[u8])], seekable: bool, level: i32,
            forced_id: Option<i32>, forced_len: Option<i32>) -> Pick {
    if let (Some(i), Some(l)) = (forced_id, forced_len) {
        return Pick { id: i, len: l };
    }
    // 4 MB per block when seekable: a smaller sample cannot show what a BWT
    // codec does to a 100 MB stream.
    let scap = if seekable { 4 << 20 } else { 96 << 20 };
    let mut maxstream = 0usize;
    for (l, i) in samples {
        maxstream = maxstream.max(l.len()).max(i.len());
    }
    let big = maxstream >= 0x7fff_ffff;

    let mut zi = 0usize;
    let mut zl = 0usize;
    let mut ci = [0usize; CAND_FAST.len()];
    let mut cl = [0usize; CAND_FAST.len()];
    const BAD: usize = usize::MAX / 2;

    for (lens, ids) in samples {
        let nl = lens.len().min(scap);
        let ni = ids.len().min(scap);
        if ni > 0 { zi += ffi::zstd(&ids[..ni], level).len(); }
        if nl > 0 { zl += ffi::zstd(&lens[..nl], level).len(); }
        for (k, &c) in CAND_FAST.iter().enumerate() {
            // Size alone does not qualify a codec: it must decode back.
            ci[k] = accumulate(ci[k], &ids[..ni], c, BAD);
            cl[k] = accumulate(cl[k], &lens[..nl], c, BAD);
        }
    }

    let (mut wi, mut wl) = (0i32, 0i32);
    let (mut bi, mut bl) = (zi, zl);
    // The 0.5% margin, in the C's exact arithmetic.
    let need_i = (0.995_f64 * zi as f64) as usize;
    let need_l = (0.995_f64 * zl as f64) as usize;
    for (k, &c) in CAND_FAST.iter().enumerate() {
        let ok = (!seekable || fast_enough(c)) && (!big || c & ffi::CODEC_BZIP3 != 0);
        if zi > 0 && ci[k] < bi && ci[k] < need_i && ok { bi = ci[k]; wi = c; }
        if zl > 0 && cl[k] < bl && cl[k] < need_l && ok { bl = cl[k]; wl = c; }
    }
    Pick { id: forced_id.unwrap_or(wi), len: forced_len.unwrap_or(wl) }
}

fn accumulate(acc: usize, src: &[u8], codec: i32, bad: usize) -> usize {
    if src.is_empty() { return acc; }
    if acc >= bad { return bad; }
    match ffi::paycomp(src, codec) {
        Some(t) => match ffi::paydecomp(&t, codec, src.len()) {
            Some(v) if v == src => acc + t.len(),
            _ => bad,
        },
        None => bad,
    }
}

/// The lengths frame: a forced codec, else zstd.
pub fn comp_lens(src: &[u8], codec: i32, level: i32) -> Option<Vec<u8>> {
    if codec != 0 { ffi::paycomp(src, codec) } else { Some(ffi::zstd(src, level)) }
}

/// The ids frame, same rule.
pub fn comp_ids(src: &[u8], codec: i32, level: i32) -> Option<Vec<u8>> {
    if codec != 0 { ffi::paycomp(src, codec) } else { Some(ffi::zstd(src, level)) }
}

/* ------------------------------------------------- the name stream */


/// 0 zstd, 1 integer deltas, 2 the CRAM tokeniser.
pub struct BestNames {
    pub bytes: Vec<u8>,
    pub kind: u8,
}

/// Names that are integers but not exact row indices become zigzag deltas.
/// Refuses on an empty name, a leading zero, a non-digit, or a short count,
/// which is what makes the choice safe to reverse.
pub fn int_names(names: &[u8], n: u64) -> Option<Vec<u8>> {
    let mut v = Vec::new();
    let mut prev: i64 = 0;
    let mut p = 0usize;
    for _ in 0..n {
        if p >= names.len() || names[p] == b'\n' { return None; }
        if names[p] == b'0' && p + 1 < names.len() && names[p + 1] != b'\n' { return None; }
        let mut x: i64 = 0;
        while p < names.len() && names[p] != b'\n' {
            let c = names[p];
            if !c.is_ascii_digit() || x > 1i64 << 50 { return None; }
            x = x * 10 + (c - b'0') as i64;
            p += 1;
        }
        if p >= names.len() { return None; }
        p += 1;
        let d = x - prev;
        prev = x;
        kst::put(&mut v, kst::zz(d));
    }
    if p != names.len() { return None; }
    Some(v)
}

/// The tokeniser, over the chunk sizes the C tries. In fast mode it tries
/// only the ends of the list, which is what the default encode does.
pub fn tok_best(names: &[u8]) -> Option<Vec<u8>> {
    const TCAND: [usize; 4] = [250_000, 1_000_000, 3_000_000, 9_000_000];
    let mut best: Option<Vec<u8>> = None;
    for (i, &chunk) in TCAND.iter().enumerate() {
        if i > 0 && i + 1 < TCAND.len() { continue; }   // fast mode: ends only
        if let Some(b) = tok_chunked(names, chunk) {
            if best.as_ref().is_none_or(|x| b.len() < x.len()) { best = Some(b); }
        }
    }
    best
}

/// One tokeniser blob: a chunk table, then each chunk's output.
fn tok_chunked(names: &[u8], chunk: usize) -> Option<Vec<u8>> {
    if names.is_empty() { return None; }
    let mut offs = vec![0usize];
    let mut p = 0usize;
    while p < names.len() {
        let mut c = 0usize;
        let mut e = p;
        while e < names.len() && c < chunk {
            if names[e] == b'\n' { c += 1; }
            e += 1;
        }
        offs.push(e);
        p = e;
    }
    let nchunk = offs.len() - 1;
    let mut blobs = Vec::with_capacity(nchunk);
    for i in 0..nchunk {
        blobs.push(ffi::tok_encode(&names[offs[i]..offs[i + 1]])?);
    }
    // chunk size, chunk count, then the raw length of each chunk and the
    // coded length of each, then the blobs. c/karyozip.c tokc().
    let mut out = Vec::new();
    out.extend_from_slice(&(chunk as u32).to_le_bytes());
    out.extend_from_slice(&(nchunk as u32).to_le_bytes());
    for i in 0..nchunk {
        out.extend_from_slice(&((offs[i + 1] - offs[i]) as u32).to_le_bytes());
    }
    for b in &blobs {
        out.extend_from_slice(&(b.len() as u32).to_le_bytes());
    }
    for b in &blobs { out.extend_from_slice(b); }
    Some(out)
}

/// Smallest of zstd, integer deltas and the tokeniser. The losers are
/// dropped. zstd and the tokeniser are independent and neither gains from
/// threads on a buffer this size, so they run side by side.
pub fn best_names(names: &[u8], n: u64, level: i32) -> BestNames {
    let (z, t) = rayon::join(
        || ffi::zstd(names, level),
        || tok_best(names),
    );
    let i = int_names(names, n).map(|v| ffi::zstd(&v, level));
    let mut out = BestNames { bytes: z, kind: 0 };
    if let Some(v) = i {
        if !v.is_empty() && v.len() < out.bytes.len() { out = BestNames { bytes: v, kind: 1 }; }
    }
    if let Some(v) = t {
        if !v.is_empty() && v.len() < out.bytes.len() { out = BestNames { bytes: v, kind: 2 }; }
    }
    out
}

/// A repeated sequence name would make the index ambiguous, so it is refused.
pub fn check_unique(names: &[u8]) -> Result<(), String> {
    let mut seen: HashSet<u64> = HashSet::new();
    let mut exact: HashSet<&[u8]> = HashSet::new();
    for line in names.split(|&c| c == b'\n') {
        if line.is_empty() { continue; }
        let h = kst::fnv(FNV0, line);
        if !seen.insert(h) && !exact.insert(line) {
            return Err(format!("sequence {} appears twice",
                               String::from_utf8_lossy(line)));
        }
        exact.insert(line);
    }
    Ok(())
}