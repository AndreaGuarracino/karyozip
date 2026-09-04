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
const DROPPED: u64 = 1 << 60;
const SIDEFLAG: u64 = 1 << 57;
const ALNFLAG: u64 = 1 << 59;
const DBFLAG: u64 = 1 << 58;
const LENFLAG: u64 = 1 << 56;
const PERMFLAG: u64 = 1 << 54;
const PERMTOK: u64 = 1 << 53;
const TAILFLAG: u64 = 1 << 55;
const ICHUNK: u64 = 65536;
const TAILSKIP: usize = 64;

pub struct Options {
    pub block: u64,
    pub level: i32,
    pub levelset: bool,
    pub slow: bool,
    pub keep_names: bool,
    pub drop_names: bool,
    pub names_out: Option<String>,
    pub tail: Option<bool>,
    pub tok_chunk: Option<usize>,
    pub sort: bool,
    pub id_codec: Option<i32>,
    pub len_codec: Option<i32>,
}

impl Default for Options {
    fn default() -> Self {
        Options { block: 262144, level: 17, levelset: false, slow: false,
                  keep_names: false, drop_names: false, names_out: None,
                  tail: None, tok_chunk: None, sort: false,
                  id_codec: None, len_codec: None }
    }
}

/// A compressed block: lens frame, ids frame, both raw sizes, rows, bases
/// before it, and which length spelling won.
type Done = (Vec<u8>, Vec<u8>, usize, usize, u64, u64, u64);

/// One staged block, while its payload waits in the staging file.
struct BlkMeta { off: u64, len: usize, lens: usize, rows: u64, bp: u64,
                 aoff: [usize; 2], alen: [usize; 2] }

pub struct Writer {
    opt: Options,
    lens: Vec<u8>,
    ids: Vec<u8>,
    names: Vec<u8>,
    ibx: Vec<u8>,
    nidx: u64,
    ordok: bool,
    labtab: Vec<Vec<u8>>,
    labmap: HashMap<Vec<u8>, u32, kst::Fnv>,
    lastid: Option<u32>,
    cur: Vec<u8>,
    currun: u64,
    curbp: u64,
    blkbp0: u64,
    pend: u64,
    nblk: u64,
    rows: u64,
    plen: u64,
    pval: u64,
    phas: bool,
    tailsamp: Vec<u64>,
    tailpred: u64,
    tailon: Option<bool>,
    /// --sort holds every record until close: (name, runs as (len, label id)).
    srec: Vec<SRec>,
    sortbp: u64,
    natbuf: Vec<u8>,
    perm_on: bool,
    natord: bool,
    feathash: Option<u64>,
    /// --aln: the alignment, the query it is on, and the two alternate
    /// length spellings each block also keeps.
    aln: Option<crate::aln::Aln>,
    alnq: Option<usize>,
    alnhit: u64,
    lalt: [Vec<u8>; 2],
    rawalt: [Vec<u8>; 2],
    prev_rev: bool,
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
    /// --db binds the track to the features.tsv it came from.
    pub fn set_feathash(&mut self, h: u64) { self.feathash = Some(h); }

    /// The reference's labels share the track's table: loadref interns them
    /// first, so the ids the reference answers with are the track's ids.
    pub fn set_aln(&mut self, a: crate::aln::Aln) {
        for l in a.labs.clone() { self.intern(&l); }
        self.aln = Some(a);
    }
    pub fn aln_hash(&self) -> u64 { self.aln.as_ref().map_or(0, |a| a.hash) }
    pub fn aln_hit(&self) -> u64 { self.alnhit }

    pub fn new(opt: Options, tmpdir: &std::path::Path) -> std::io::Result<Self> {
        let perm_on = opt.sort && opt.names_out.is_none();
        Ok(Writer {
            perm_on,
            opt, lens: Vec::new(), ids: Vec::new(), names: Vec::new(), ibx: Vec::new(),
            nidx: 0, ordok: true, labtab: Vec::new(), labmap: HashMap::default(), lastid: None,
            cur: Vec::new(), currun: 0, curbp: 0, blkbp0: 0, pend: 0, nblk: 0, rows: 0,
            plen: 0, pval: 0, phas: false, tailsamp: Vec::new(), tailpred: 0, tailon: None,
            srec: Vec::new(), sortbp: 0, natbuf: Vec::new(), natord: true, feathash: None,
            aln: None, alnq: None, alnhit: 0,
            lalt: [Vec::new(), Vec::new()], rawalt: [Vec::new(), Vec::new()],
            prev_rev: false,
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

    /// The lengths frame, plus the two alternate spellings --aln keeps.
    fn lensput(&mut self, v: u64) {
        kst::put(&mut self.lens, v);
        if self.aln.is_some() {
            kst::put(&mut self.lalt[0], self.plen);
            let x = if self.prev_rev { self.plen } else { v };
            kst::put(&mut self.lalt[1], x);
        }
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
                    .map(|&v| v.abs_diff(self.tailpred))
                    .sum();
                self.tailon = Some(self.opt.tail
                    .unwrap_or(dev * 12 <= self.tailpred * TAILSKIP as u64));
            }
            self.lensput(self.pval);
            self.phas = false;
            return;
        }
        if self.tailon != Some(true) {
            self.lensput(self.pval);
            self.phas = false;
            return;
        }
        let z = kst::zz(self.curbp as i64 - self.tailpred as i64);
        kst::put(&mut self.lens, z);
        if self.aln.is_some() {
            kst::put(&mut self.lalt[0], z);
            kst::put(&mut self.lalt[1], z);
        }
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
        let mut aoff = [0usize; 2];
        let mut alen = [0usize; 2];
        for v in 0..2 {
            aoff[v] = self.rawalt[v].len();
            alen[v] = self.lalt[v].len();
            let take = std::mem::take(&mut self.lalt[v]);
            self.rawalt[v].extend_from_slice(&take);
        }
        let off = self.rawlen;
        self.raw.seek(SeekFrom::Start(off))?;
        self.raw.write_all(&self.lens)?;
        self.raw.write_all(&self.ids)?;
        let len = self.lens.len() + self.ids.len();
        self.bm.push(BlkMeta { off, len, lens: self.lens.len(), rows: self.pend, bp: self.blkbp0, aoff, alen });
        self.rawlen += len as u64;
        self.nblk += 1;
        self.lens.clear();
        self.ids.clear();
        self.pend = 0;
        Ok(())
    }

    /// --sort: hold the record; the order is decided at close.
    fn sortpush(&mut self, name: &[u8], run_len: u64, label: &[u8], start: u64)
        -> Result<(), String> {
        let newrec = start == 0
            || self.srec.last().is_none_or(|r| r.0 != name);
        if newrec && start != 0 { return Err("a sequence does not start at 0".into()); }
        if !newrec && start != self.sortbp {
            return Err("gap or overlap: start is not the previous end".into());
        }
        self.sortbp = if newrec { run_len } else { self.sortbp + run_len };
        if newrec {
            if self.perm_on {
                let want = self.srec.len().to_string();
                if self.natord && want.as_bytes() != name { self.natord = false; }
                self.natbuf.extend_from_slice(name);
                self.natbuf.push(b'\n');
            }
            self.srec.push((name.to_vec(), Vec::new()));
        }
        let id = self.intern(label);
        self.srec.last_mut().unwrap().1.push((run_len, id));
        Ok(())
    }

    /// Order by label tuple, then run count, then lengths, then input rank.
    fn sortflush(&mut self) -> Result<(), String> {
        let lab = std::mem::take(&mut self.labtab);
        let mut idx: Vec<usize> = (0..self.srec.len()).collect();
        idx.sort_by(|&a, &b| {
            let (p, q) = (&self.srec[a].1, &self.srec[b].1);
            let n = p.len().min(q.len());
            for i in 0..n {
                if p[i].1 != q[i].1 { return p[i].1.cmp(&q[i].1); }
            }
            if p.len() != q.len() { return p.len().cmp(&q.len()); }
            for i in 0..n {
                if p[i].0 != q[i].0 { return p[i].0.cmp(&q[i].0); }
            }
            a.cmp(&b)
        });
        let recs = std::mem::take(&mut self.srec);
        self.labtab = lab;
        let n = recs.len();
        for &i in &idx {
            let (nm, runs) = &recs[i];
            let rank = i.to_string();
            let name: &[u8] = if self.perm_on { rank.as_bytes() } else { nm };
            let mut pos = 0u64;
            for &(len, id) in runs {
                let l = self.labtab[id as usize].clone();
                self.row_inner(name, len, &l, pos)?;
                pos += len;
            }
        }
        eprintln!("sorted {n} records by label tuple then lengths");
        Ok(())
    }

    pub fn row(&mut self, name: &[u8], run_len: u64, label: &[u8], start: u64)
        -> Result<(), String> {
        if self.opt.sort { return self.sortpush(name, run_len, label, start); }
        self.row_inner(name, run_len, label, start)
    }

    fn row_inner(&mut self, name: &[u8], run_len: u64, label: &[u8], start: u64)
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
            self.alnq = self.aln.as_ref().and_then(|a| a.find(name));
            self.currun = 0;
            self.curbp = 0;
            if self.pend > 0 && self.pend >= self.opt.block { self.flush().map_err(io)?; }
        }
        if start != self.curbp {
            return Err("gap or overlap: start is not the previous end".into());
        }
        // The C interns first and then decides, so the label joins the table
        // even when the reference answers for it.
        let mut id = self.intern(label);
        if let Some(a) = self.aln.take() {
            if label == b"=" { return Err("a label \"=\" is reserved under --aln".into()); }
            if let Some(q) = self.alnq {
                let rl = a.row(q, start, run_len);
                if rl != u32::MAX && rl == id {
                    id = self.intern(b"=");
                    self.alnhit += 1;
                }
            }
            self.aln = Some(a);
        }
        let id = id as u64;
        if self.phas { self.lensput(self.pval); self.phas = false; }
        if self.pend >= self.opt.block { self.flush().map_err(io)?; }
        if self.pend == 0 { self.blkbp0 = self.curbp; }
        self.plen = run_len;
        self.phas = true;
        self.pval = run_len;
        if let (Some(a), Some(q)) = (&self.aln, self.alnq) {
            let (e, rev) = a.end(q, start);
            self.prev_rev = rev;
            if e == start + run_len { self.pval = 0; }
        }
        kst::put(&mut self.ids, id);
        self.pend += 1;
        self.currun += 1;
        self.rows += 1;
        self.curbp += run_len;
        Ok(())
    }

    pub fn close(mut self, out: &mut impl Write) -> Result<Stats, String> {
        if self.opt.sort && !self.srec.is_empty() { self.sortflush()?; }
        self.tail();
        self.pushread();
        self.flush().map_err(io)?;
        let level = self.opt.level;

        // The codec search samples the first, middle and last block, which is
        // why it cannot run until every row has arrived. --slow reads them all.
        let rawtot: usize = self.bm.iter().map(|b| b.len).sum();
        let get = |i: usize| self.read_block(i).expect("staged block");
        let _s = crate::prof::span(8);
        let maxstream = self.bm.iter().map(|m| m.lens.max(m.len - m.lens)).max().unwrap_or(0);
        let p = pick(&get, self.bm.len(), rawtot, maxstream, level, &self.opt);
        drop(_s);
        let level = p.level;
        let (lvl_id, lvl_len) = (if p.lvl_id < 0 { level } else { p.lvl_id },
                                 if p.lvl_len < 0 { level } else { p.lvl_len });
        let aln_on = self.aln.is_some();
        let fast_aln = aln_on && !self.opt.slow;

        // Compress apart, write in order, a batch at a time. Holding every
        // compressed block would cost the whole file in memory; the C stages
        // them in a temp file for the same reason.
        let _s = crate::prof::span(1);
        let mut sizes = Vec::new();
        let mut body = tmpfile(&self.tmpdir).map_err(io)?;
        let (mut byt_lens, mut byt_labs) = (0u64, 0u64);
        let batch = rayon::current_num_threads() * 4;
        for chunk in (0..self.bm.len()).collect::<Vec<_>>().chunks(batch) {
            let done: Vec<Done> = chunk
                .par_iter()
                .map(|&i| {
                    let (lens, ids) = self.read_block(i).expect("staged block");
                    // Under --aln a block keeps the smallest of three length
                    // spellings: every sentinel, plus-strand sentinels only,
                    // or plain lengths. Zeros do not always compress better.
                    let m = &self.bm[i];
                    let mut vs: [Option<&[u8]>; 3] = [Some(&lens), None, None];
                    if aln_on {
                        for v in 0..2 {
                            let alt = &self.rawalt[v][m.aoff[v]..m.aoff[v] + m.alen[v]];
                            if alt.is_empty() || (alt.len() == lens.len() && alt == &lens[..]) { continue; }
                            if v == 1 {
                                if let Some(x) = vs[1] { if x == alt { continue; } }
                            }
                            vs[v + 1] = Some(alt);
                        }
                    }
                    const ORDER: [usize; 3] = [1, 2, 0];
                    let (mut best, mut bv, mut fa) = (0usize, usize::MAX, Vec::new());
                    for &v in &ORDER {
                        let Some(src) = vs[v] else { continue };
                        let f2 = comp_lens(src, p.len, lvl_len).expect("lens frame");
                        let a2 = f2.len();
                        if bv == usize::MAX || a2 + if fast_aln { a2 / 100 } else { 0 } < best {
                            best = a2; bv = v; fa = f2;
                        }
                    }
                    let nl = vs[bv].unwrap().len();
                    let b = comp_ids(&ids, p.id, lvl_id).expect("ids frame");
                    (fa, b, nl, ids.len(), m.rows, m.bp, bv as u64)
                })
                .collect();
            for (a, b, nl, ni, rows, bp, nosent) in &done {
                body.write_all(a).map_err(io)?;
                body.write_all(b).map_err(io)?;
                byt_lens += a.len() as u64;
                byt_labs += b.len() as u64;
                kst::put(&mut sizes, a.len() as u64);
                kst::put(&mut sizes, b.len() as u64);
                kst::put(&mut sizes, *ni as u64);
                kst::put(&mut sizes, if aln_on { (*nl as u64) << 2 | nosent } else { *nl as u64 });
                kst::put(&mut sizes, *rows);
                kst::put(&mut sizes, *bp);
                kst::put(&mut sizes, kst::fnv(kst::fnv(FNV0, a), b));
            }
        }

        drop(_s);
        let _s = crate::prof::span(2);
        let mut labs = Vec::new();
        for l in &self.labtab { labs.extend_from_slice(l); labs.push(b'\n'); }

        // The C's decision, step for step: which of the names are stored,
        // as text, as ordinals, in a sidecar, or as ranks plus a blob.
        let mut perm_on = self.perm_on;
        let mut ordinal = self.ordok && self.nidx > 0;
        if self.opt.keep_names { ordinal = false; }
        if !self.opt.drop_names && self.nidx > 1 {
            check_unique(if perm_on { &self.natbuf } else { &self.names })?;
        }
        let mut dropnames = self.opt.drop_names;
        if dropnames && perm_on { dropnames = false; ordinal = false; }
        else if dropnames { ordinal = true; }
        if !ordinal && self.ordok && self.names.is_empty() && self.nidx > 0 {
            for i in 0..self.nidx {
                self.names.extend_from_slice(i.to_string().as_bytes());
                self.names.push(b'\n');
            }
        }
        if perm_on && !self.natord && ordinal && !self.opt.drop_names {
            self.names = std::mem::take(&mut self.natbuf);
            ordinal = false;
            perm_on = false;
        }
        if ordinal { self.names.clear(); }

        let (mut nz, mut tokused, mut intused) = if ordinal {
            (Vec::new(), false, false)
        } else {
            let b = best_names(&self.names, self.nidx, level);
            (b.bytes, b.kind == 2, b.kind == 1)
        };
        // --sort keeps identity by storing a rank per record. Storing the
        // sorted names instead is sometimes smaller; the C compares and so
        // does this.
        let mut permblob: Vec<u8> = Vec::new();
        let mut permtok = false;
        if perm_on && self.opt.drop_names { self.natbuf.clear(); }
        else if perm_on && !self.natord && !self.natbuf.is_empty() {
            let nt: Vec<&[u8]> = self.natbuf.split(|&c| c == b'\n')
                .filter(|l| !l.is_empty()).collect();
            if nt.len() as u64 != self.nidx {
                return Err("input-order names do not match the records".into());
            }
            let mut srt = Vec::new();
            for line in self.names.split(|&c| c == b'\n') {
                if line.is_empty() { continue; }
                let v: u64 = std::str::from_utf8(line).ok()
                    .and_then(|t| t.parse().ok())
                    .ok_or("a stored rank is not a record number")?;
                if v >= self.nidx { return Err("a stored rank is not a record number".into()); }
                srt.extend_from_slice(nt[v as usize]);
                srt.push(b'\n');
            }
            let b = best_names(&srt, self.nidx, level);
            let p1 = ffi::zstd(&self.natbuf, level);
            let p2 = tok_best(&self.natbuf);
            let (pbest, ptok) = match p2 {
                Some(t) if t.len() < p1.len() => (t, true),
                _ => (p1, false),
            };
            if !b.bytes.is_empty() && b.bytes.len() <= nz.len() + pbest.len() {
                nz = b.bytes;
                tokused = b.kind == 2;
                intused = b.kind == 1;
                perm_on = false;
            } else {
                permblob = pbest;
                permtok = ptok;
            }
        }
        let npl = if perm_on { permblob.len() } else { 0 };

        drop(_s);
        let _s = crate::prof::span(3);
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
        h[1] = if self.opt.names_out.is_some() { 0 } else { nz.len() as u64 }
             | if tokused { TOKFLAG } else { 0 }
             | if dropnames && !(self.ordok && self.nidx > 0) { DROPPED } else { 0 }
             | if self.opt.names_out.is_some() { SIDEFLAG } else { 0 }
             | if intused { INTFLAG } else { 0 }
             | if npl > 0 || (perm_on && self.opt.drop_names) { PERMFLAG } else { 0 }
             | if permtok { PERMTOK } else { 0 }
             | if self.feathash.is_some() { DBFLAG } else { 0 }
             | if aln_on { ALNFLAG | LENFLAG } else { 0 }
             | if self.tailon == Some(true) { TAILFLAG } else { 0 };
        h[2] = ichunks.len() as u64;
        h[3] = sz.len() as u64;
        h[4] = self.rows;
        h[5] = p.len as u32 as u64 | self.feathash.map_or(0, |f| (f & 0xffff_ffff) << 32);
        h[6] = lz.len() as u64;
        h[7] = self.labtab.len() as u64;
        h[8] = self.opt.block;
        h[9] = p.id as u32 as u64 | self.feathash.map_or(0, |f| (f >> 32) << 32);
        h[10] = nichunk;
        h[11] = itz.len() as u64;
        h[12] = self.tailpred;
        h[13] = TAILSKIP as u64 | if aln_on { self.aln_hash() << 8 } else { 0 };
        h[14] = npl as u64;
        let mut hb = Vec::with_capacity(128);
        for v in h { hb.extend_from_slice(&v.to_le_bytes()); }
        let mut mh = kst::fnv(FNV0, &hb);
        mh = kst::fnv(mh, &lz);
        mh = kst::fnv(mh, &nz);
        mh = kst::fnv(mh, &ichunks);
        mh = kst::fnv(mh, &sz);
        mh = kst::fnv(mh, &itz);
        mh = kst::fnv(mh, &permblob[..npl]);
        h[15] = mh;
        hb.truncate(120);
        hb.extend_from_slice(&h[15].to_le_bytes());

        drop(_s);
        let _s = crate::prof::span(7);
        out.write_all(MAGIC).map_err(io)?;
        out.write_all(&hb).map_err(io)?;
        out.write_all(&lz).map_err(io)?;
        match &self.opt.names_out {
            Some(f) => {
                let mut nf = std::fs::File::create(f).map_err(io)?;
                nf.write_all(&nz).map_err(io)?;
            }
            None => out.write_all(&nz).map_err(io)?,
        }
        out.write_all(&ichunks).map_err(io)?;
        out.write_all(&sz).map_err(io)?;
        out.write_all(&itz).map_err(io)?;
        out.write_all(&permblob[..npl]).map_err(io)?;
        body.seek(SeekFrom::Start(0)).map_err(io)?;
        std::io::copy(&mut body, out).map_err(io)?;
        out.flush().map_err(io)?;
        drop(_s);

        Ok(Stats {
            rows: self.rows, reads: self.nidx, blocks: self.nblk,
            names: nz.len(), index: ichunks.len(), offs: sz.len(),
            lens: byt_lens, labs: byt_labs, raw: rawtot,
            tokused, intused, ordinal,
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

/// --slow: rANS order 0/1 x {PACK,RLE,X32,STRIPE}, arith order 0/1 x
/// {PACK,RLE}, bzip3 at six block sizes. A strict superset of the default.
pub const CAND_SLOW: [i32; 46] = [
    0x100, 0x80, 0x40, 0xc0, 0x8, 0x88, 0x48, 0xc8,
    0x4, 0x84, 0x44, 0xc4, 0xc, 0x8c, 0x4c, 0xcc,
    0x1, 0x81, 0x41, 0xc1, 0x9, 0x89, 0x49, 0xc9,
    0x5, 0x85, 0x45, 0xc5, 0xd, 0x8d, 0x4d, 0xcd,
    ffi::CODEC_ARITH, ffi::CODEC_ARITH | 0x80, ffi::CODEC_ARITH | 0x40, ffi::CODEC_ARITH | 0xc0,
    ffi::CODEC_ARITH | 1, ffi::CODEC_ARITH | 0x81, ffi::CODEC_ARITH | 0x41, ffi::CODEC_ARITH | 0xc1,
    ffi::CODEC_BZIP3 | 8, ffi::CODEC_BZIP3 | 16, ffi::CODEC_BZIP3 | 32,
    ffi::CODEC_BZIP3 | 64, ffi::CODEC_BZIP3 | 128, ffi::CODEC_BZIP3 | 256,
];

/// bzip3 decodes at about 8 MB/s against zstd's ~1000, so it is the wrong
/// pick for a file a query will seek into.
fn fast_enough(codec: i32) -> bool { codec & ffi::CODEC_BZIP3 == 0 }

/// One --sort record: the name, then its (length, label id) runs.
type SRec = (Vec<u8>, Vec<(u64, u32)>);

pub struct Pick { pub id: i32, pub len: i32, pub lvl_id: i32, pub lvl_len: i32, pub level: i32 }

/// `blocks` yields every block as (lens, ids); the default search reads three
/// of them, --slow reads them all, twice.
pub fn pick(all: &dyn Fn(usize) -> (Vec<u8>, Vec<u8>), nraw: usize, rawtot: usize,
            maxstream: usize, mut level: i32, o: &Options) -> Pick {
    let (block, slow, levelset) = (o.block, o.slow, o.levelset);
    let (forced_id, forced_len) = (o.id_codec, o.len_codec);
    let cand: &[i32] = if slow { &CAND_SLOW } else { &CAND_FAST };
    let nc = cand.len();
    if let (Some(i), Some(l)) = (forced_id, forced_len) {
        return Pick { id: i, len: l, lvl_id: -1, lvl_len: -1, level };
    }
    let seekable = nraw > 1 || rawtot as u64 > block;
    let mut picks: Vec<usize> = Vec::new();
    if nraw > 0 {
        picks.push(0);
        if nraw > 2 { picks.push(nraw / 2); }
        if nraw > 1 { picks.push(nraw - 1); }
    }
    // 4 MB per block when seekable; --slow uses whole blocks, since a 4 MB
    // sample cannot show what a BWT codec does to a 100 MB stream.
    let mut scap = if seekable { 4usize << 20 } else { 96 << 20 };
    if slow { scap = maxstream.clamp(4 << 20, 2048 << 20); }
    let big = maxstream >= 0x7fff_ffff;

    // A candidate the selection cannot take is not worth compressing. On a
    // seekable track that is both bzip3 entries, whose bz3_new alone asks
    // for its 16 MB or 64 MB block.
    let usable = |c: i32| (slow || !seekable || fast_enough(c))
                       && (!big || c & ffi::CODEC_BZIP3 != 0);
    let (mut zi, mut zl, mut z22i, mut z22l) = (0usize, 0usize, 0usize, 0usize);
    let mut ci = vec![0usize; nc];
    let mut cl = vec![0usize; nc];
    const BAD: usize = usize::MAX / 2;
    for &j in &picks {
        let (lens, ids) = all(j);
        let nl = lens.len().min(scap);
        let ni = ids.len().min(scap);
        let z22 = slow && !levelset;
        let step: Vec<(usize, usize)> = (0..nc + 2).into_par_iter()
            .map(|k| match k {
                _ if k < nc => if usable(cand[k]) {
                        (accumulate(0, &ids[..ni], cand[k], BAD),
                         accumulate(0, &lens[..nl], cand[k], BAD))
                    } else { (BAD, BAD) },
                _ if k == nc => (if ni > 0 { ffi::zstd(&ids[..ni], level).len() } else { 0 },
                                 if nl > 0 { ffi::zstd(&lens[..nl], level).len() } else { 0 }),
                _ if z22 => (if ni > 0 { ffi::zstd(&ids[..ni], 22).len() } else { 0 },
                             if nl > 0 { ffi::zstd(&lens[..nl], 22).len() } else { 0 }),
                _ => (0, 0),
            }).collect();
        zi += step[nc].0; zl += step[nc].1;
        z22i += step[nc + 1].0; z22l += step[nc + 1].1;
        for k in 0..nc {
            ci[k] = if ci[k] >= BAD || step[k].0 >= BAD { BAD } else { ci[k] + step[k].0 };
            cl[k] = if cl[k] >= BAD || step[k].1 >= BAD { BAD } else { cl[k] + step[k].1 };
        }
    }

    let (mut wi, mut wl) = (0i32, 0i32);
    let (mut bi, mut bl) = (zi, zl);
    if slow && !levelset && z22i < bi { bi = z22i; }
    if slow && !levelset && z22l < bl { bl = z22l; }
    let lvl_i = if slow && !levelset && z22i < zi { 22 } else { level };
    let lvl_l = if slow && !levelset && z22l < zl { 22 } else { level };
    for (k, &c) in cand.iter().enumerate() {
        let ok = (slow || !seekable || fast_enough(c)) && (!big || c & ffi::CODEC_BZIP3 != 0);
        let need_i = if slow { bi } else { (0.995_f64 * zi as f64) as usize };
        let need_l = if slow { bl } else { (0.995_f64 * zl as f64) as usize };
        let m = usize::from(slow);
        if zi > 0 && ci[k] < bi && ci[k] < need_i + m && ok { bi = ci[k]; wi = c; }
        if zl > 0 && cl[k] < bl && cl[k] < need_l + m && ok { bl = cl[k]; wl = c; }
    }
    if slow && !levelset && wi == 0 && wl == 0 {
        level = if lvl_i == 22 || lvl_l == 22 { 22 } else { level };
    } else if slow && !levelset && (wi == 0 || wl == 0) {
        level = if wi == 0 { lvl_i } else { lvl_l };
    }

    let (mut lvl_id, mut lvl_len) = (-1i32, -1i32);
    if slow {
        // Decide on the whole file, and drop any candidate that fails to
        // decode a block back.
        let mut ti = vec![0usize; nc + 2];
        let mut tl = vec![0usize; nc + 2];
        let mut live: Vec<bool> = (0..nc + 2).map(|k| if k < nc {
            (ci[k] < usize::MAX / 4 || cl[k] < usize::MAX / 4)
                && (!big || cand[k] & ffi::CODEC_BZIP3 != 0)
        } else { k == nc || !levelset }).collect();
        for b in 0..nraw {
            let (lens, ids) = all(b);
            let step: Vec<(usize, usize)> = (0..nc + 2).into_par_iter().map(|k| {
                if !live[k] { return (0, 0); }
                if k < nc {
                    (sz(ffi::paycomp(&ids, cand[k])), sz(ffi::paycomp(&lens, cand[k])))
                } else {
                    let lv = if k == nc { level } else { 22 };
                    (ffi::zstd(&ids, lv).len(), ffi::zstd(&lens, lv).len())
                }
            }).collect();
            for k in 0..nc + 2 { ti[k] += step[k].0; tl[k] += step[k].1; }
        }
        let (mut qi, mut ql);
        loop {
            qi = nc; ql = nc;
            for k in 0..nc + 2 {
                if !live[k] { continue; }
                if ti[k] < ti[qi] { qi = k; }
                if tl[k] < tl[ql] { ql = k; }
            }
            let (mut bad_i, mut bad_l) = (false, false);
            for b in 0..nraw {
                if bad_i && bad_l { break; }
                let (lens, ids) = all(b);
                if qi < nc && !ids.is_empty() && !bad_i && !verifies(&ids, cand[qi]) { bad_i = true; }
                if ql < nc && !lens.is_empty() && !bad_l && !verifies(&lens, cand[ql]) { bad_l = true; }
            }
            if !bad_i && !bad_l { break; }
            if bad_i { live[qi] = false; }
            if bad_l { live[ql] = false; }
        }
        wi = if qi < nc { cand[qi] } else { 0 };
        wl = if ql < nc { cand[ql] } else { 0 };
        if wi == 0 { lvl_id = if qi == nc + 1 { 22 } else { level }; }
        if wl == 0 { lvl_len = if ql == nc + 1 { 22 } else { level }; }
        eprintln!("codec: --slow ran {} live candidates over {nraw} blocks whole", nc + 2);
    }
    eprintln!("codec: len=0x{:x} id=0x{:x} from {} sample blocks, {}",
              forced_len.unwrap_or(wl), forced_id.unwrap_or(wi), picks.len(),
              if seekable { "seekable" } else { "one block" });
    Pick { id: forced_id.unwrap_or(wi), len: forced_len.unwrap_or(wl), lvl_id, lvl_len, level }
}

fn sz(v: Option<Vec<u8>>) -> usize { v.map_or(usize::MAX / 4, |x| x.len()) }

fn verifies(src: &[u8], codec: i32) -> bool {
    match ffi::paycomp(src, codec) {
        Some(t) => ffi::paydecomp(&t, codec, src.len()).as_deref() == Some(src),
        None => false,
    }
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