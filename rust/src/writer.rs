//! The KST writer. Rows arrive grouped by sequence, gapless, each sequence
//! starting at 0; blocks fill and are staged; the file is assembled at close.
//! Every layout decision here is c/karyozip.c's, because the two must produce
//! the same bytes.

use crate::ffi;
use crate::kst::{self, FNV0, TOKFLAG, INTFLAG, FLDFLAG, DROPPED, ALNFLAG, DBFLAG, SIDEFLAG, LENFLAG,
                 TAILFLAG, RUNFLAG, ILVFLAG, GRPFLAG, PERMFLAG, PERMTOK, CARTFLAG, ICHUNK};
use rayon::prelude::*;
use std::collections::HashMap;
use std::collections::HashSet;
use std::io::{Read, Seek, SeekFrom, Write};

pub const MAGIC: &[u8; 8] = b"KST00001";
const TAILSKIP: usize = 64;

pub struct Options {
    pub block: u64,
    pub level: i32,
    pub levelset: bool,
    pub slow: bool,
    pub whole: bool,       // -e high: decide on the whole file
    pub drop_names: bool,
    pub names_out: Option<String>,
    pub sort: bool,
    /// --sort=MEM: bytes held before a spill, 0 for none.
    pub sortmem: u64,
}

impl Default for Options {
    fn default() -> Self {
        Options { block: u64::MAX, level: 17, levelset: false, slow: true, whole: false, sortmem: 0,
                  drop_names: false, names_out: None,
                  sort: false }
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
    /// fnv of each label and of cur, for the row hash
    labh: Vec<u64>,
    curnh: u64,
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
    /// --sort writes ranks; the PAF still wants the real name.
    alnname: Option<Vec<u8>>,
    /// --sort=MEM: the bytes held, the runs already spilled, the next rank.
    held: u64,
    spill: Vec<std::fs::File>,
    natn: u64,
    natord: bool,
    feathash: Option<u64>,
    /// the input rows as text, stored so the track can serve as --ref
    rowhash: u64,
    cart: bool,
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
    pub fn set_feathash(&mut self, h: u64, cart: bool) { self.feathash = Some(h); self.cart = cart; }

    /// The reference's labels share the track's table: loadref interns them
    /// first, so the ids the reference answers with are the track's ids.
    pub fn set_aln(&mut self, a: crate::aln::Aln) {
        for l in a.labs.clone() { self.intern(&l); }
        self.aln = Some(a);
    }
    pub fn aln_hash(&self) -> u64 { self.aln.as_ref().map_or(0, |a| a.hash) }

    pub fn new(opt: Options, tmpdir: &std::path::Path) -> std::io::Result<Self> {
        let perm_on = opt.sort && opt.names_out.is_none();
        Ok(Writer {
            perm_on, alnname: None, held: 0, spill: Vec::new(), natn: 0,
            opt, lens: Vec::new(), ids: Vec::new(), names: Vec::new(), ibx: Vec::new(),
            nidx: 0, ordok: true, labtab: Vec::new(), labmap: HashMap::default(), lastid: None,
            cur: Vec::new(), currun: 0, curbp: 0, blkbp0: 0, pend: 0, nblk: 0, rows: 0,
            plen: 0, pval: 0, phas: false, tailsamp: Vec::new(), tailpred: 0, tailon: None,
            srec: Vec::new(), sortbp: 0, natbuf: Vec::new(), natord: true, feathash: None, cart: false, rowhash: kst::FNV0, labh: Vec::new(), curnh: 0,
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
        self.labh.push(kst::fnv(kst::FNV0, label));
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
                self.tailon = Some(dev * 12 <= self.tailpred * TAILSKIP as u64);
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
        if name.is_empty() { return Err("empty sequence name".into()); }
        if name.len() >= 512 { return Err("name too long".into()); }
        if run_len == 0 {
            return Err("end equals start: a run must cover at least one base".into());
        }
        if start > u64::MAX - run_len {
            return Err("a coordinate does not fit in 64 bits".into());
        }
        if run_len > 0xffff_ffff { return Err("run longer than 4 GB".into()); }
        let newrec = start == 0
            || self.srec.last().is_none_or(|r| r.nm != name);
        if newrec && start != 0 { return Err("a sequence does not start at 0".into()); }
        if !newrec && start != self.sortbp {
            return Err("gap or overlap: start is not the previous end".into());
        }
        self.sortbp = if newrec { run_len } else { self.sortbp + run_len };
        if newrec {
            if self.perm_on {
                let want = self.natn.to_string();
                if self.natord && want.as_bytes() != name { self.natord = false; }
                self.natbuf.extend_from_slice(name);
                self.natbuf.push(b'\n');
            }
            if self.natn > u32::MAX as u64 { return Err("--sort takes at most 2^32 records".into()); }
            self.srec.push(SRec { nm: name.to_vec(), runs: Vec::new(), idx: self.natn as u32 });
            self.natn += 1;
            self.held += name.len() as u64 + 1 + 24;   // the C's srec
            // Spill at a record boundary only, holding the one just opened.
            if self.opt.sortmem > 0 && self.held > self.opt.sortmem {
                let last = self.srec.pop().unwrap();
                self.spill_write()?;
                self.held = last.nm.len() as u64 + 1 + 24;
                self.srec.push(last);
            }
        }
        let id = self.intern(label);
        self.held += 8;
        self.srec.last_mut().unwrap().runs.push((run_len as u32, id));
        Ok(())
    }

    /// One sorted run to a temp file: nrun, name length, rank, name, runs.
    fn spill_write(&mut self) -> Result<(), String> {
        if self.srec.is_empty() { return Ok(()); }
        self.srec.sort_by(srec_cmp);
        let mut f = tmpfile(&self.tmpdir)
            .map_err(|_| "cannot create a sort temp file: check --tmp-dir".to_string())?;
        let mut buf = Vec::with_capacity(1 << 20);
        for r in &self.srec {
            buf.extend_from_slice(&(r.runs.len() as u32).to_le_bytes());
            buf.extend_from_slice(&(r.nm.len() as u32).to_le_bytes());
            buf.extend_from_slice(&r.idx.to_le_bytes());
            buf.extend_from_slice(&r.nm);
            for &(l, b) in &r.runs {
                buf.extend_from_slice(&l.to_le_bytes());
                buf.extend_from_slice(&b.to_le_bytes());
            }
            if buf.len() > (1 << 20) {
                f.write_all(&buf).map_err(|_| "sort spill write failed: is --tmp-dir full?".to_string())?;
                buf.clear();
            }
        }
        f.write_all(&buf).map_err(|_| "sort spill write failed: is --tmp-dir full?".to_string())?;
        f.flush().map_err(|_| "sort spill write failed: is --tmp-dir full?".to_string())?;
        eprintln!("spill {}: {} records", self.spill.len(), self.srec.len());
        self.spill.push(f);
        self.srec.clear();
        self.held = 0;
        Ok(())
    }

    /// The k-way merge over the spilled runs, in the same order.
    fn spillmerge(&mut self) -> Result<(), String> {
        let mut files = std::mem::take(&mut self.spill);
        let mut head: Vec<Option<SRec>> = Vec::with_capacity(files.len());
        for f in files.iter_mut() {
            f.seek(SeekFrom::Start(0)).map_err(io)?;
            head.push(srd_next(f)?);
        }
        let mut out = 0u64;
        loop {
            let mut best = usize::MAX;
            for i in 0..head.len() {
                if head[i].is_none() { continue; }
                if best == usize::MAX
                   || srec_cmp(head[i].as_ref().unwrap(), head[best].as_ref().unwrap())
                      == std::cmp::Ordering::Less { best = i; }
            }
            if best == usize::MAX { break; }
            let r = head[best].take().unwrap();
            self.sort_emit(&r)?;
            out += 1;
            head[best] = srd_next(&mut files[best])?;
        }
        eprintln!("external sort: {out} records merged from {} runs", files.len());
        Ok(())
    }

    /// One sorted record to the writer, under its rank when names are permuted.
    fn sort_emit(&mut self, r: &SRec) -> Result<(), String> {
        let rank = r.idx.to_string();
        let name: &[u8] = if self.perm_on { rank.as_bytes() } else { &r.nm };
        self.alnname = Some(r.nm.clone());
        let mut pos = 0u64;
        for &(len, id) in &r.runs {
            let l = self.labtab[id as usize].clone();
            self.row_inner(name, len as u64, &l, pos)?;
            pos += len as u64;
        }
        Ok(())
    }

    /// Order by label tuple, then run count, then lengths, then input rank.
    fn sortflush(&mut self) -> Result<(), String> {
        if self.opt.sortmem > 0 {   // the C always merges from disk
            self.spill_write()?;
            return if self.spill.is_empty() { Ok(()) } else { self.spillmerge() };
        }
        let mut recs = std::mem::take(&mut self.srec);
        recs.sort_by(srec_cmp);
        let n = recs.len();
        for r in &recs { self.sort_emit(r)?; }
        eprintln!("sorted {n} records by label tuple then lengths");
        Ok(())
    }

    /// The CLI gets these for free: its parser split on tabs and newlines,
    /// so a name or label can never hold one, and a zero-length run cannot
    /// reach the writer. A library caller can produce all of them, and a
    /// newline in a label makes a track that does not decode, because the
    /// label table is newline separated. Refusing here keeps encode
    /// lossless.
    pub fn row(&mut self, name: &[u8], run_len: u64, label: &[u8], start: u64)
        -> Result<(), String> {
        let bad = name.iter().find(|&&c| matches!(c, 0 | b'\n' | b'\t'))
            .or_else(|| label.iter().find(|&&c| matches!(c, 0 | b'\n')));
        match bad {
            Some(0) => return Err("a NUL byte inside a row".into()),
            Some(b'\n') => return Err("a newline inside a name or label".into()),
            Some(_) => return Err("a tab inside a sequence name".into()),
            None => {}
        }
        self.push(name, run_len, label, start)
    }

    /// The same row without those checks, for callers whose parser already
    /// split on tabs and newlines. This is the CLI's path, as it is the C's.
    pub fn push(&mut self, name: &[u8], run_len: u64, label: &[u8], start: u64)
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
            self.curnh = kst::fnv(kst::FNV0, name);
            // --sort writes ranks; the PAF still wants the real name.
            let ln: &[u8] = self.alnname.as_deref().unwrap_or(name);
            self.alnq = self.aln.as_ref().and_then(|a| a.find(ln));
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
        self.rowhash = kst::rowmix(self.rowhash, self.curnh, start, start + run_len, self.labh[id as usize]);
        if self.aln.is_some() {
            if label == b"=" { return Err("a label \"=\" is reserved under --aln".into()); }
            // borrowed, not taken: taking moved the whole Aln out and back per row
            let hit = self.alnq.is_some_and(|q| {
                let rl = self.aln.as_ref().unwrap().row(q, start, run_len);
                rl != u32::MAX && rl == id
            });
            if hit { id = self.intern(b"="); self.alnhit += 1; }
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
        // why it cannot run until every row has arrived. -e high reads them all.
        let rawtot: usize = self.bm.iter().map(|b| b.len).sum();
        let get = |i: usize| self.read_block(i).unwrap_or_else(|_| crate::bail("cannot read back a staged block"));
        let _s = crate::prof::span(8);
        let maxstream = self.bm.iter().map(|m| m.lens.max(m.len - m.lens)).max().unwrap_or(0);
        let onelay = self.bm.len() == 1 && !self.ibx.is_empty();   // the layout search ranks it all
        let p = pick(&get, self.bm.len(), rawtot, maxstream, level, &self.opt, onelay);
        drop(_s);
        let level = p.level;
        let (mut lvl_id, mut lvl_len) = (if p.lvl_id < 0 { level } else { p.lvl_id },
                                         if p.lvl_len < 0 { level } else { p.lvl_len });
        let mut idc = p.id;
        let aln_on = self.aln.is_some();

        // Compress apart, write in order, a batch at a time. Holding every
        // compressed block would cost the whole file in memory; the C stages
        // them in a temp file for the same reason.
        // The run count can head each read's lengths instead of living in the
        // index. One block only, since a split read has no count to head it,
        // and measured against the index it would replace.
        let mut altbuf: Vec<u8> = Vec::new();
        let (mut runhead, mut ilv, mut grp, mut alt_lens) = (false, false, false, 0usize);
        let mut ilvcodec = p.len;
        let mut altv: [Option<(Vec<u8>, usize)>; 3] = [None, None, None];
        let _l = crate::prof::span(9);
        if self.bm.len() == 1 && !self.ibx.is_empty() {
            let (lens0, ids0) = self.read_block(0).unwrap_or_else(|_| crate::bail("cannot read back a staged block"));
            // The layout is decided on the mode's sample, as the codec is: whole
            // reads up to 4 MB of lengths -e low, DEFSAMPLE MB the default, all
            // of them -e high. Each layout gets its own best codec there; only
            // the winner is built whole, compressed once and read back.
            let (slow, whole) = (self.opt.slow, self.opt.whole);
            let cap = if whole { usize::MAX } else { (if slow { DEFSAMPLE } else { 4 }) << 20 };
            let (mut pl, mut pi, mut pb, mut x, mut y, mut z) = (0usize, 0usize, 0usize, 0usize, 0usize, 0usize);
            while z < self.ibx.len() {
                let (k, n) = vget(&self.ibx, z); z = n;
                let mut j = 0;
                while j < k && x < lens0.len() && y < ids0.len() {
                    x = crate::reader::skip(&lens0, x); y = crate::reader::skip(&ids0, y); j += 1;
                }
                if x > cap && pl > 0 { break; }         // at least one read
                pl = x; pi = y; pb = z;
            }
            let (ls, is, bs) = (&lens0[..pl], &ids0[..pi], &self.ibx[..pb]);
            let idx = idxsize(bs, if level < 19 { level } else { 19 });
            let lay = [None, splice_heads(ls, bs), if slow { grp_pack(ls, is, bs) } else { None }, ilv_pack(ls, is, bs)];
            // the labels and the four length spellings, ranked in one pool; a
            // spelling's whole size is its sample's scaled by the lengths'
            let e: Vec<u8> = Vec::new();
            let src: Vec<&[u8]> = vec![is, ls, lay[1].as_deref().unwrap_or(&e), lay[2].as_deref().unwrap_or(&e), lay[3].as_deref().unwrap_or(&e)];
            let scale = |m: usize| if pl > 0 { (m as u128 * lens0.len() as u128 / pl as u128) as usize } else { m };
            let n5 = [ids0.len(), lens0.len(), scale(src[2].len()), scale(src[3].len()), scale(src[4].len())];
            let r = rank_all(&src, &n5, &self.opt);
            let idok = pi == 0 || r[0].2 > 0;
            let mut sz = [0usize; 4];
            if r[1].2 > 0 && idok { sz[0] = r[1].2 + r[0].2 + idx; }   // plain pays its index
            if r[2].2 > 0 && idok { sz[1] = r[2].2 + r[0].2; }
            if r[3].2 > 0 && idok { sz[2] = r[3].2 + r[0].2; }
            if r[4].2 > 0 { sz[3] = r[4].2; }                          // interleaved: no ids frame
            let mut l = 0;
            for k in 1..4 { if sz[k] > 0 && (sz[l] == 0 || sz[k] < sz[l]) { l = k; } }
            (idc, lvl_id) = (r[0].0, if r[0].1 < 0 { level } else { r[0].1 });
            (ilvcodec, lvl_len) = (r[l + 1].0, if r[l + 1].1 < 0 { level } else { r[l + 1].1 });
            let hb = if l > 0 { splice_heads(&lens0, &self.ibx) } else { None };
            match l {
                1 => { altbuf = hb.unwrap_or_default(); alt_lens = altbuf.len(); }
                2 => { altbuf = grp_pack(&lens0, &ids0, &self.ibx).unwrap_or_default(); alt_lens = altbuf.len(); grp = true; }
                3 => { altbuf = ilv_pack(&lens0, &ids0, &self.ibx).unwrap_or_default(); alt_lens = hb.map_or(0, |h| h.len()); ilv = true; }
                _ => {}
            }
            runhead = l > 0;
            if l > 0 && lens_ok(&altbuf, ilvcodec, lvl_len).is_none() { ilvcodec = 0; }   // zstd comes back
            let labels = if l == 3 { String::new() } else { format!(", labels {}", cname(idc)) };
            eprintln!("layout: {}, lengths {}{labels}, decided on {pl} of {} B of lengths",
                      ["plain", "headed", "grouped", "interleaved"][l], cname(ilvcodec), lens0.len());
            // --aln keeps three length spellings: respell every one the same way
            if runhead || ilv {
                altv[0] = Some((std::mem::take(&mut altbuf), alt_lens));
                for v in 0..2 {
                    let m = &self.bm[0];
                    if m.alen[v] == 0 { continue; }
                    let alt = &self.rawalt[v][m.aoff[v]..m.aoff[v] + m.alen[v]];
                    let Some(ah) = (if grp { grp_pack(alt, &ids0, &self.ibx) }
                                    else { splice_heads(alt, &self.ibx) }) else { continue };
                    let ahn = ah.len();
                    if !ilv { altv[v + 1] = Some((ah, ahn)); continue; }
                    if let Some(ai) = ilv_pack(alt, &ids0, &self.ibx) { altv[v + 1] = Some((ai, ahn)); }
                }
            }
        }
        if !(self.bm.len() == 1 && !self.ibx.is_empty()) {   // blocks: no layout search, the codec search's pick
            eprintln!("layout: plain, lengths {}, labels {}, decided on the codec search's sample", cname(p.len), cname(p.id));
        }
        drop(_l);
        let lencodec = ilvcodec;
        if runhead { self.ibx.clear(); }

        let _s = crate::prof::span(1);
        let mut sizes = Vec::new();
        let mut body = tmpfile(&self.tmpdir).map_err(io)?;
        let (mut byt_lens, mut byt_labs) = (0u64, 0u64);
        let batch = rayon::current_num_threads() * 4;
        for chunk in (0..self.bm.len()).collect::<Vec<_>>().chunks(batch) {
            let done: Vec<Done> = chunk
                .par_iter()
                .map(|&i| {
                    let (lens, ids) = self.read_block(i).unwrap_or_else(|_| crate::bail("cannot read back a staged block"));
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
                    if runhead || ilv {
                        for v in 0..3 {
                            match &altv[v] {
                                Some((b, _)) => vs[v] = Some(b),
                                None => if v > 0 { vs[v] = None; },
                            }
                        }
                    }
                    const ORDER: [usize; 3] = [1, 2, 0];
                    let (mut best, mut bv, mut fa) = (0usize, usize::MAX, Vec::new());
                    for &v in &ORDER {
                        let Some(src) = vs[v] else { continue };
                        let f = comp_frame(src, lencodec, lvl_len).unwrap_or_else(|| crate::bail("block compression failed"));
                        if bv == usize::MAX || f.len() < best { best = f.len(); bv = v; fa = f; }
                    }
                    let nl = if ilv { altv[bv].as_ref().unwrap().1 } else { vs[bv].unwrap().len() };
                    let b = if ilv { Vec::new() } else { comp_frame(&ids, idc, lvl_id).unwrap_or_else(|| crate::bail("block compression failed")) };
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
        if !self.opt.drop_names && self.nidx > 1 {
            check_unique(if perm_on { &self.natbuf } else { &self.names })?;
        }
        let mut dropnames = self.opt.drop_names;
        // Dropping names that are not dense ordinals leaves --aln with no way
        // back to the PAF, so they are kept as integer deltas instead.
        if dropnames && aln_on && !ordinal {
            eprintln!("--names drop: the names are not dense ordinals and --aln \
                       needs them to reach the PAF, so they are stored as \
                       integer deltas");
            dropnames = false;
        }
        if dropnames && perm_on { dropnames = false; ordinal = false; }
        else if dropnames { ordinal = true; }
        if !ordinal && self.ordok && self.names.is_empty() && self.nidx > 0 {
            for i in 0..self.nidx {
                self.names.extend_from_slice(i.to_string().as_bytes());
                self.names.push(b'\n');
            }
        }
        if ordinal && self.opt.names_out.is_some() && !self.opt.drop_names {
            return Err("--names drop:FILE: the names are already the record ordinals, so a \
                        file beside the track would be redundant: the track rebuilds them".into());
        }
        if perm_on && !self.natord && ordinal && !self.opt.drop_names {
            self.names = std::mem::take(&mut self.natbuf);
            ordinal = false;
            perm_on = false;
        }
        if ordinal { self.names.clear(); }

        // Rank the names codecs on whole names up to 4 MB -e low, DEFSAMPLE
        // MB the default, all -e high, then run the winner over the whole
        // buffer: a level 17 zstd of 9.4M ONT read names is 130 s on its own.
        let (mut nz, mut tokused, mut intused, mut fldused) = if ordinal {
            (Vec::new(), false, false, false)
        } else {
            let (mut q, mut smpn) = (0usize, 0u64);
            let ncap = if self.opt.whole { usize::MAX } else { (if self.opt.slow { DEFSAMPLE } else { 4 }) << 20 };
            while q < self.names.len() && q < ncap {
                match memchr::memchr(b'\n', &self.names[q..]) {
                    Some(k) => { q += k + 1; smpn += 1; }
                    None => break,
                }
            }
            let b = best_names(&self.names[..q], smpn, level, self.opt.slow);
            let r = if q == self.names.len() {
                (b.bytes, b.kind == 2, b.kind == 1, b.kind == 3)   // the sample was all of it
            } else {
                let whole = match b.kind {
                    1 => int_names(&self.names, self.nidx).map(|v| ffi::zstd(&v, level)),
                    2 => tok_chunked(&self.names, b.chunk),
                    3 => fld_names(&self.names, self.nidx),
                    _ => None,
                };
                match whole {
                    Some(v) if !v.is_empty() => (v, b.kind == 2, b.kind == 1, b.kind == 3),
                    _ => (ffi::zstd(&self.names, level), false, false, false),   // refused
                }
            };
            let k = if r.2 { "integers" } else if r.1 { "tokeniser" } else if r.3 { "fields" } else { "zstd" };
            let ck = if r.1 { format!(" chunk {}", b.chunk) } else { String::new() };
            eprintln!("names: {k}{ck}, decided on {q} of {} B", self.names.len());
            r
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
            let b = best_names(&srt, self.nidx, level, self.opt.slow);
            let p1 = ffi::zstd(&self.natbuf, level);
            let p2 = tok_best(&self.natbuf, self.opt.slow).map(|t| t.0);
            let (pbest, ptok) = match p2 {
                Some(t) if t.len() < p1.len() => (t, true),
                _ => (p1, false),
            };
            if !b.bytes.is_empty() && b.bytes.len() <= nz.len() + pbest.len() {
                nz = b.bytes;
                tokused = b.kind == 2;
                intused = b.kind == 1;
                fldused = b.kind == 3;
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
        {   // chunk bounds first, then every chunk compressed on its own thread
            let (mut q, mut cpend, mut seen) = (0usize, 0u64, 0u64);
            let (mut cb, mut cpe) = (vec![0usize], vec![0u64]);
            while q < self.ibx.len() {
                let (nruns, next) = vget(&self.ibx, q);
                q = next;
                cpend += nruns;
                seen += 1;
                if seen % ICHUNK == 0 || q >= self.ibx.len() { cb.push(q); cpe.push(cpend); }
            }
            let czs: Vec<Vec<u8>> = (0..cb.len() - 1).into_par_iter()
                .map(|i| ffi::zstd(&self.ibx[cb[i]..cb[i + 1]], meta)).collect();
            for (i, cz) in czs.iter().enumerate() {
                ichunks.extend_from_slice(cz);
                kst::put(&mut itab, cz.len() as u64);
                kst::put(&mut itab, cpe[i]);
                nichunk += 1;
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
             | if fldused { FLDFLAG } else { 0 }
             | if npl > 0 || (perm_on && self.opt.drop_names) { PERMFLAG } else { 0 }
             | if permtok { PERMTOK } else { 0 }
             | if self.feathash.is_some() { DBFLAG } else { 0 }
             | if self.cart { CARTFLAG } else { 0 }
             | if aln_on { ALNFLAG | LENFLAG } else { 0 }
             | if self.tailon == Some(true) { TAILFLAG } else { 0 }
             | if runhead { RUNFLAG } else { 0 }
             | if ilv { ILVFLAG } else { 0 }
             | if grp { GRPFLAG } else { 0 };
        h[2] = ichunks.len() as u64;
        h[3] = sz.len() as u64;
        h[4] = self.rows;
        h[5] = lencodec as u32 as u64 | self.feathash.map_or(0, |f| (f & 0xffff_ffff) << 32);
        h[6] = lz.len() as u64;
        h[7] = self.labtab.len() as u64;
        h[8] = self.opt.block;
        h[9] = idc as u32 as u64 | self.feathash.map_or(0, |f| (f >> 32) << 32);
        h[10] = nichunk;
        h[11] = itz.len() as u64;
        h[12] = self.tailpred;
        h[13] = TAILSKIP as u64 | (if aln_on { self.aln_hash() } else { self.rowhash } << 8);
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
            tokused, intused, fldused, ordinal, dropped: dropnames, alnhit: self.alnhit,
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
    pub tokused: bool, pub intused: bool, pub fldused: bool, pub ordinal: bool, pub dropped: bool,
    pub alnhit: u64,
}

use crate::kst::vget;

fn io(e: std::io::Error) -> String { e.to_string() }


/* ------------------------------------------------- choosing a codec */
/// One list. The default tries all of it and spends whatever that costs;
/// -e low is the subset marked true, the codecs whose encode and decode stay
/// reasonable, so it can never reach a candidate the default did not.
/// Codec 0 means zstd, so bare rANS order 0 is spelled 0x100.
pub const KCAND: [(i32, bool); 51] = [
    (0x100, true), (0x80, false), (0x40, false), (0xc0, false),
    (0x8, false), (0x88, false), (0x48, false), (0xc8, false),
    (0x4, false), (0x84, false), (0x44, false), (0xc4, false),
    (0xc, false), (0x8c, false), (0x4c, false), (0xcc, false),
    (0x1, true), (0x81, false), (0x41, true), (0xc1, false),
    (0x9, false), (0x89, false), (0x49, false), (0xc9, false),
    (0x5, false), (0x85, false), (0x45, false), (0xc5, false),
    (0xd, false), (0x8d, false), (0x4d, false), (0xcd, false),
    // STRIPE carries its stripe count in bits 8 and up
    (0x208, false), (0x209, true), (0x408, false), (0x409, false),
    (ffi::CODEC_ARITH, true), (ffi::CODEC_ARITH | 0x80, false),
    (ffi::CODEC_ARITH | 0x40, false), (ffi::CODEC_ARITH | 0xc0, false),
    (ffi::CODEC_ARITH | 1, true), (ffi::CODEC_ARITH | 0x81, false),
    (ffi::CODEC_ARITH | 0x41, false), (ffi::CODEC_ARITH | 0xc1, false),
    (ffi::CODEC_BZIP3 | 8, false), (ffi::CODEC_BZIP3 | 16, true), (ffi::CODEC_BZIP3 | 32, false),
    (ffi::CODEC_BZIP3 | 64, true), (ffi::CODEC_BZIP3 | 128, false), (ffi::CODEC_BZIP3 | 256, false),
    (ffi::CODEC_LZMA, false),
];

/// A bzip3 block bigger than the data is the same file as the smallest that
/// covers it, so only that one is worth a try. The mode's own list decides.
fn bz3_dup(c: i32, m: usize, n: usize, cand: &[i32]) -> bool {
    if c & ffi::CODEC_BZIP3 == 0 { return false; }
    let mb = ((c & 0xffff) as usize) << 20;
    let blocks = cand.iter().filter(|&&o| o & ffi::CODEC_BZIP3 != 0).map(|&o| ((o & 0xffff) as usize) << 20);
    let keep = blocks.clone().filter(|&o| o >= n).min().or_else(|| blocks.max()).unwrap_or(0);
    mb >= m && mb != keep
}

/// a codec id as a word, for the diagnostics
fn cname(c: i32) -> String {
    if c == 0 { return "zstd".into(); }
    if c & ffi::CODEC_LZMA != 0 { return "lzma".into(); }
    if c & ffi::CODEC_BZIP3 != 0 { return format!("bzip3-{}MB", c & 0xffff); }
    format!("{}-o{}{}{}", if c & ffi::CODEC_ARITH != 0 { "arith" } else { "rans" }, c & 1,
            if c & 0x80 != 0 { "+pack" } else { "" }, if c & 0x40 != 0 { "+rle" } else { "" })
}

/// the mode's candidates
/// MB the default ranks every codec on
pub const DEFSAMPLE: usize = 16;

pub fn candlist(slow: bool) -> Vec<i32> {
    KCAND.iter().filter(|(_, f)| slow || *f).map(|(c, _)| *c).collect()
}

/// bzip3 decodes at about 8 MB/s against zstd's ~1000, so it is the wrong
/// pick for a file a query will seek into.
// a query decodes a whole block: bzip3 at 8 MB/s, arith at 2.5x rANS
fn fast_enough(codec: i32) -> bool {
    codec & (ffi::CODEC_BZIP3 | ffi::CODEC_ARITH | ffi::CODEC_LZMA) == 0
}

/// One --sort record: the name, its (length, label id) runs, and the input
/// rank, which travels with it across a spill.
struct SRec { nm: Vec<u8>, runs: Vec<(u32, u32)>, idx: u32 }

/// Labels lead, then the run count, then the lengths, then the input rank.
fn srec_cmp(a: &SRec, b: &SRec) -> std::cmp::Ordering {
    let n = a.runs.len().min(b.runs.len());
    for i in 0..n {
        if a.runs[i].1 != b.runs[i].1 { return a.runs[i].1.cmp(&b.runs[i].1); }
    }
    if a.runs.len() != b.runs.len() { return a.runs.len().cmp(&b.runs.len()); }
    for i in 0..n {
        if a.runs[i].0 != b.runs[i].0 { return a.runs[i].0.cmp(&b.runs[i].0); }
    }
    a.idx.cmp(&b.idx)
}

pub struct Pick { pub id: i32, pub len: i32, pub lvl_id: i32, pub lvl_len: i32, pub level: i32 }

/// `blocks` yields every block as (lens, ids); the default search reads three
/// of them, -e high reads them all, twice.
pub fn pick(all: &dyn Fn(usize) -> (Vec<u8>, Vec<u8>), nraw: usize, rawtot: usize,
            maxstream: usize, level: i32, o: &Options, onelay: bool) -> Pick {
    let (block, slow, whole, levelset) = (o.block, o.slow, o.whole, o.levelset);
    let cv = candlist(slow); let cand: &[i32] = &cv;
    let nc = cand.len();
    if onelay { return Pick { id: 0, len: 0, lvl_id: -1, lvl_len: -1, level }; }   // the layout search ranks it all
    let seekable = nraw > 1 || rawtot as u64 > block;
    let mut picks: Vec<usize> = Vec::new();
    if nraw > 0 {
        picks.push(0);
        if nraw > 2 { picks.push(nraw / 2); }
        if nraw > 1 { picks.push(nraw - 1); }
    }
    // one rule, three samples: -e low its codecs on 4 MB, the default all of
    // them on DEFSAMPLE MB, -e high all of them on the whole file
    let mut scap = if slow { DEFSAMPLE << 20 } else { 4usize << 20 };
    if whole { scap = maxstream.clamp(4 << 20, 2048 << 20); }
    let big = maxstream >= 0x7fff_ffff;

    // A candidate the selection cannot take is not worth compressing. On a
    // seekable track that is both bzip3 entries, whose bz3_new alone asks
    // for its 16 MB or 64 MB block.
    let usable = |c: i32| (!seekable || fast_enough(c))
                       && (!big || c & (ffi::CODEC_BZIP3 | ffi::CODEC_LZMA) != 0);
    // bzip3 splits the input into blocks of its MB argument, so every size at
    // or above the stream gives one block and the same bytes: the smallest
    // such size stands for all of them
    // hts_pack and its arith twin refuse above 16 distinct bytes and the
    // encoder then clears the PACK bit, so such a candidate is its plain sibling
    let bz3dup = |c: i32, n: usize| bz3_dup(c, n, maxstream, &cand);
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
                _ if k < nc => if usable(cand[k]) && !bz3dup(cand[k], ni.max(nl)) {
                        (if pack_dup(cand[k], &ids[..ni]) { BAD }
                         else { packsize(&ids[..ni], cand[k], BAD) },
                         if pack_dup(cand[k], &lens[..nl]) { BAD }
                         else { packsize(&lens[..nl], cand[k], BAD) })
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

    let (wi, wl);
    let (mut lvl_id, mut lvl_len) = (-1i32, -1i32);
    {   // every mode: the smallest that round trips on what it was ranked on
        let mut ti = vec![0usize; nc + 2];
        let mut tl = vec![0usize; nc + 2];
        // a candidate lives per stream: one that does not round-trip the
        // ids may still be the smallest for the lengths
        let ok = |k: usize| !big || cand[k] & (ffi::CODEC_BZIP3 | ffi::CODEC_LZMA) != 0;
        let mut live_i: Vec<bool> = (0..nc + 2).map(|k| if k < nc {
            ci[k] < usize::MAX / 4 && ok(k)
        } else { k == nc || (!levelset && slow) }).collect();
        let mut live_l: Vec<bool> = (0..nc + 2).map(|k| if k < nc {
            cl[k] < usize::MAX / 4 && ok(k)
        } else { k == nc || (!levelset && slow) }).collect();
        // one block whose streams fit the sample: the sample already
        // compressed exactly this, so its sizes are the file's
        let reuse = !whole || (nraw == 1 && scap >= maxstream);
        if reuse {
            for k in 0..nc { ti[k] = ci[k]; tl[k] = cl[k]; }
            ti[nc] = zi; tl[nc] = zl;
            ti[nc + 1] = z22i; tl[nc + 1] = z22l;
        }
        for b in (0..nraw).take_while(|_| !reuse) {
            let (lens, ids) = all(b);
            let step: Vec<(usize, usize)> = (0..nc + 2).into_par_iter().map(|k| {
                if !live_i[k] && !live_l[k] { return (0, 0); }
                // An empty stream costs nothing and is not compressed, as
                // the C's `if (ni)` and `if (nl)` guards have it.
                let each = |s: &[u8], f: &dyn Fn(&[u8]) -> usize| {
                    if s.is_empty() { 0 } else { f(s) }
                };
                if k < nc {
                    (if live_i[k] { each(&ids, &|s| sz(ffi::paycomp(s, cand[k]))) } else { 0 },
                     if live_l[k] { each(&lens, &|s| sz(ffi::paycomp(s, cand[k]))) } else { 0 })
                } else {
                    let lv = if k == nc { level } else { 22 };
                    (each(&ids, &|s| ffi::zstd(s, lv).len()),
                     each(&lens, &|s| ffi::zstd(s, lv).len()))
                }
            }).collect();
            for k in 0..nc + 2 { ti[k] += step[k].0; tl[k] += step[k].1; }
        }
        let (mut qi, mut ql);
        loop {
            qi = nc; ql = nc;
            for k in 0..nc + 2 {
                if live_i[k] && ti[k] < ti[qi] { qi = k; }
                if live_l[k] && tl[k] < tl[ql] { ql = k; }
            }
            let (mut bad_i, mut bad_l) = (false, false);
            let over: Vec<usize> = if whole { (0..nraw).collect() } else { picks.clone() };
            for &b in &over {
                if bad_i && bad_l { break; }
                let (lens, ids) = all(b);
                let cap = if whole { usize::MAX } else { scap };   // the sample it was ranked on
                let (lens, ids) = (&lens[..lens.len().min(cap)], &ids[..ids.len().min(cap)]);
                if qi < nc && !ids.is_empty() && !bad_i && !verifies(ids, cand[qi]) { bad_i = true; }
                if ql < nc && !lens.is_empty() && !bad_l && !verifies(lens, cand[ql]) { bad_l = true; }
            }
            if !bad_i && !bad_l { break; }
            if bad_i { live_i[qi] = false; }
            if bad_l { live_l[ql] = false; }
        }
        wi = if qi < nc { cand[qi] } else { 0 };
        wl = if ql < nc { cand[ql] } else { 0 };
        if wi == 0 { lvl_id = if qi == nc + 1 { 22 } else { level }; }
        if wl == 0 { lvl_len = if ql == nc + 1 { 22 } else { level }; }
        eprintln!("codec: {} ranked {} candidates on {}", if !slow { "effort low" } else if whole { "effort high" } else { "effort medium" },
                  nc + 2, if whole { "the whole file" } else { "the sample" });
    }
    eprintln!("codec: len=0x{:x} id=0x{:x} from {} sample blocks, {}",
              wl, wi, picks.len(),
              if seekable { "seekable" } else { "one block" });
    Pick { id: wi, len: wl, lvl_id, lvl_len, level }
}

fn sz(v: Option<Vec<u8>>) -> usize { v.map_or(usize::MAX / 4, |x| x.len()) }

fn verifies(src: &[u8], codec: i32) -> bool {
    match ffi::paycomp(src, codec) {
        Some(t) => ffi::paydecomp(&t, codec, src.len()).as_deref() == Some(src),
        None => false,
    }
}

/// Size alone does not qualify a codec, so the sample round trips it; under
/// slow the whole-file pass round trips the winner and drops it on failure,
/// so the sample only sizes.
/// Labels counted down from 127 so they do not share bytes with the small
/// lengths beside them. Its own inverse, and the varint keeps its width.
pub fn ilv_flip(v: u64) -> u64 { if v < 128 { 127 - v } else { v } }

/// One stream a read: its run count, then a label beside its own length.
fn ilv_pack(lens: &[u8], ids: &[u8], ibx: &[u8]) -> Option<Vec<u8>> {
    let mut out = Vec::with_capacity(lens.len() + ids.len() + ibx.len() + 16);
    let (mut p, mut q, mut r) = (0usize, 0usize, 0usize);
    while r < ibx.len() {
        let (k, n) = vget(ibx, r); r = n;
        kst::put(&mut out, k);
        for _ in 0..k {
            if q >= ids.len() || p >= lens.len() { return None; }
            let (lb, n2) = vget(ids, q); kst::put(&mut out, ilv_flip(lb)); q = n2;
            let n2 = crate::reader::skip(lens, p); out.extend_from_slice(&lens[p..n2]); p = n2;
        }
    }
    if p != lens.len() || q != ids.len() { None } else { Some(out) }
}

/// Each read's lengths in label order: all of label 0's, then label 1's. The
/// label stream stays as it was, so it says where each one goes back.
fn grp_pack(lens: &[u8], ids: &[u8], ibx: &[u8]) -> Option<Vec<u8>> {
    let mut out = Vec::with_capacity(lens.len() + ibx.len() + 16);
    let (mut p, mut q, mut r) = (0usize, 0usize, 0usize);
    let mut kv: Vec<u64> = Vec::new();
    let mut lv: Vec<(usize, usize)> = Vec::new();
    while r < ibx.len() {
        let (k, n) = vget(ibx, r); r = n;
        kst::put(&mut out, k);
        if k >> 32 != 0 { return None; }
        kv.clear(); lv.clear();
        for j in 0..k {
            if p >= lens.len() || q >= ids.len() { return None; }
            let n2 = crate::reader::skip(lens, p); lv.push((p, n2)); p = n2;
            let (id, n3) = vget(ids, q); q = n3;
            if id >> 32 != 0 { return None; }
            kv.push((id << 32) | j);
        }
        kv.sort_unstable();
        for key in &kv { let (a, b) = lv[(key & 0xffff_ffff) as usize]; out.extend_from_slice(&lens[a..b]); }
    }
    if p != lens.len() || q != ids.len() { None } else { Some(out) }
}

/// the same lengths with each read's run count in front of its own
fn splice_heads(lens: &[u8], ibx: &[u8]) -> Option<Vec<u8>> {
    let mut out = Vec::with_capacity(lens.len() + ibx.len() + 16);
    let (mut p, mut q) = (0usize, 0usize);
    while q < ibx.len() {
        let (k, n) = vget(ibx, q); q = n;
        kst::put(&mut out, k);
        for _ in 0..k {
            if p >= lens.len() { return None; }
            let n2 = crate::reader::skip(lens, p);
            out.extend_from_slice(&lens[p..n2]); p = n2;
        }
    }
    if p != lens.len() { return None; } else { Some(out) }
}

/// hts_pack and its arith twin refuse above 16 distinct bytes and the
/// encoder then clears the PACK bit, so such a candidate is its plain sibling
fn pack_dup(c: i32, p: &[u8]) -> bool {
    if c & 0x80 == 0 || c & ffi::CODEC_BZIP3 != 0 { return false; }
    let (mut seen, mut d) = ([false; 256], 0);
    p.iter().any(|&b| { if !seen[b as usize] { seen[b as usize] = true; d += 1; } d > 16 })
}

/// Every stream by every codec of the mode in one pool, LZMA and bzip3 first
/// since they take longest. Per stream the smallest that round trips on its
/// sample wins, ties to the earlier of zstd, the table, zstd-22. `n` is each
/// stream's whole size; out: codec, zstd level or -1, size, 0 when none came back.
fn rank_all(src: &[&[u8]], n: &[usize], o: &Options) -> Vec<(i32, i32, usize)> {
    let cand = candlist(o.slow);
    let mut jobs: Vec<(u8, usize, usize, i32, i32)> = Vec::new();   // class, stream, rank, codec, level
    for (t, s) in src.iter().enumerate() {
        if s.is_empty() { continue; }
        jobs.push((if o.level >= 17 { 0 } else { 2 }, t, 0, 0, o.level));
        for (k, &c) in cand.iter().enumerate() {
            if pack_dup(c, s) || bz3_dup(c, s.len(), n[t], &cand) { continue; }
            if n[t] >= 0x7fff_ffff && c & (ffi::CODEC_BZIP3 | ffi::CODEC_LZMA) == 0 { continue; }   // rANS, arith: below INT_MAX
            let cls = if c & ffi::CODEC_LZMA != 0 { 0 } else if c & ffi::CODEC_BZIP3 != 0 { 1 } else { 2 };
            jobs.push((cls, t, k + 1, c, -1));
        }
        if o.slow && !o.levelset { jobs.push((0, t, cand.len() + 1, 0, 22)); }
    }
    // the slowest first, as the C: LZMA and zstd 17 and up, then bzip3, then
    // the rest, the larger stream first; taken in that order by every thread
    jobs.sort_unstable_by_key(|j| (j.0, std::cmp::Reverse(src[j.1].len()), j.1, j.2));
    use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};
    let (next, zs): (AtomicUsize, Vec<AtomicUsize>) = (AtomicUsize::new(0), jobs.iter().map(|_| AtomicUsize::new(0)).collect());
    // each stream's best so far (size, rank, job, bytes), kept for the read-back
    let best: Vec<std::sync::Mutex<Option<(usize, usize, usize, Vec<u8>)>>> = src.iter().map(|_| Default::default()).collect();
    let prof = crate::prof::ON.load(Relaxed);
    let jt: Vec<AtomicUsize> = jobs.iter().map(|_| AtomicUsize::new(0)).collect();   // ns, under --profile
    let tp = std::time::Instant::now();
    rayon::broadcast(|_| loop {
        let j = next.fetch_add(1, Relaxed);
        let Some(&(_, t, r, c, lv)) = jobs.get(j) else { break };
        let t0 = std::time::Instant::now();
        if c == 0 { zs[j].store(ffi::zstd(src[t], lv).len(), Relaxed); }
        else if let Some(v) = ffi::paycomp(src[t], c) {
            zs[j].store(v.len(), Relaxed);
            let mut b = best[t].lock().unwrap_or_else(|e| e.into_inner());
            if b.as_ref().is_none_or(|x| (v.len(), r) < (x.0, x.1)) { *b = Some((v.len(), r, j, v)); }
        }
        if prof { jt[j].store(t0.elapsed().as_nanos() as usize, Relaxed); }
    });
    if prof {   // the slowest jobs set the pool's critical path
        eprintln!("  rank pool: {} jobs, {:.3} s", jobs.len(), tp.elapsed().as_secs_f64());
        let mut o: Vec<usize> = (0..jobs.len()).collect();
        o.sort_by_key(|&j| std::cmp::Reverse(jt[j].load(Relaxed)));
        for &j in o.iter().take(5).filter(|&&j| jt[j].load(Relaxed) > 0) {
            eprintln!("  rank job: stream {} {} on {} B, {:.3} s", jobs[j].1, cname(jobs[j].3), src[jobs[j].1].len(),
                      jt[j].load(Relaxed) as f64 * 1e-9);
        }
    }
    let tp = std::time::Instant::now();
    let z: Vec<usize> = zs.into_iter().map(AtomicUsize::into_inner).collect();
    let mut kept: Vec<Option<(usize, Vec<u8>)>> = best.into_iter()
        .map(|m| m.into_inner().unwrap_or_else(|e| e.into_inner()).map(|x| (x.2, x.3))).collect();
    let picked: Vec<(i32, i32, usize)> = kept.par_iter_mut().enumerate().map(|(t, k)| {
        let mut live: Vec<usize> = (0..jobs.len()).filter(|&j| jobs[j].1 == t && z[j] > 0).collect();
        live.sort_by_key(|&j| (z[j], jobs[j].2));
        for j in live {
            let (_, _, _, c, lv) = jobs[j];
            // zstd is trusted, the rest read back; the pool's bytes when it kept them
            let ok = c == 0 || match k.take() {
                Some((kj, v)) if kj == j => ffi::paydecomp(&v, c, src[t].len()).as_deref() == Some(src[t]),
                _ => verifies(src[t], c),
            };
            if ok { return (c, lv, z[j]); }
        }
        (0, -1, 0)
    }).collect();
    if prof { eprintln!("  rank read-back: {:.3} s", tp.elapsed().as_secs_f64()); }
    picked
}

/// A respelt block is a different stream from the one the codec was chosen
/// on, so it counts only at the size it comes back at.
fn lens_ok(src: &[u8], codec: i32, level: i32) -> Option<usize> {
    let c = comp_frame(src, codec, level)?;
    let v = if codec != 0 { ffi::paydecomp(&c, codec, src.len())? } else { ffi::zstd_decompress(&c, src.len())? };
    if v == src { Some(c.len()) } else { None }
}

/// what the index blob costs, so the two spellings can be compared
fn idxsize(ibx: &[u8], meta: i32) -> usize {
    let (mut q, mut seen) = (0usize, 0u64);
    let mut cb = vec![0usize];
    while q < ibx.len() {
        let (_, next) = vget(ibx, q); q = next; seen += 1;
        if seen % ICHUNK == 0 || q >= ibx.len() { cb.push(q); }
    }
    (0..cb.len() - 1).into_par_iter().map(|i| ffi::zstd(&ibx[cb[i]..cb[i + 1]], meta).len()).sum()
}

fn packsize(src: &[u8], codec: i32, bad: usize) -> usize {
    if src.is_empty() { return bad; }
    ffi::paycomp(src, codec).map_or(bad, |t| t.len())
}

/// A lengths or labels frame: the chosen codec, else zstd.
pub fn comp_frame(src: &[u8], codec: i32, level: i32) -> Option<Vec<u8>> {
    if codec != 0 { ffi::paycomp(src, codec) } else { Some(ffi::zstd(src, level)) }
}

/* ------------------------------------------------- the name stream */


/// 0 zstd, 1 integer deltas, 2 the CRAM tokeniser.
pub struct BestNames {
    pub bytes: Vec<u8>,
    pub kind: u8,
    pub chunk: usize,
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
        if p >= names.len() || x > 1i64 << 50 { return None; }   // what the reader refuses
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
pub const TCAND: [usize; 4] = [250_000, 1_000_000, 3_000_000, 9_000_000];
/// The smallest tokeniser blob and its chunk.
pub fn tok_best(names: &[u8], slow: bool) -> Option<(Vec<u8>, usize)> {
    let mut best: Option<(Vec<u8>, usize)> = None;
    for (i, &chunk) in TCAND.iter().enumerate() {
        if !slow && i > 0 { continue; }   // -e low: the smallest chunk only
        // smallest wins; a tie keeps the smaller chunk
        if let Some(b) = tok_chunked(names, chunk) {
            if best.as_ref().is_none_or(|x| b.len() < x.0.len()) { best = Some((b, chunk)); }
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
    // One chunk per thread: 9.4M ONT read names is 38 chunks, and doing them
    // in turn was the whole of the Rust's names phase. The C threads it too.
    let blobs: Vec<Vec<u8>> = (0..nchunk).into_par_iter()
        .map(|i| ffi::tok_encode(&names[offs[i]..offs[i + 1]]))
        .collect::<Option<Vec<_>>>()?;
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
/// Read names split on ':' into fields that are either the same in every name
/// or a plain integer. The integers go down as zigzag deltas, a read's fields
/// together, and bzip3 takes the lot.

const FLDMAX: usize = 16;
/// ':' splits an Illumina name, '-' a UUID: try both, keep the smaller
pub fn fld_names(names: &[u8], n: u64) -> Option<Vec<u8>> {
    let a = fld_names1(names, n, b':');
    let b = fld_names1(names, n, b'-');
    match (a, b) {
        (Some(x), Some(y)) => Some(if y.len() < x.len() { y } else { x }),
        (x, y) => x.or(y),
    }
}

fn fld_names1(names: &[u8], n: u64, sep: u8) -> Option<Vec<u8>> {
    if n == 0 { return None; }
    let first: Vec<&[u8]> = names.split(|&c| c == b'\n').next()?.split(|&c| c == sep).collect();
    let nf = first.len();
    if nf == 0 || nf > FLDMAX { return None; }
    let (mut same, mut isnum) = (vec![true; nf], vec![true; nf]);
    let mut ishex: Vec<bool> = first.iter().map(|f| !f.is_empty() && f.len() % 2 == 0 && f.len() <= 64).collect();
    let mut seen = 0u64;
    for line in names.split(|&c| c == b'\n') {          // one pass to decide
        if seen == n { if line.is_empty() { break; } return None; }
        let f: Vec<&[u8]> = line.split(|&c| c == sep).collect();
        if f.len() != nf { return None; }
        for k in 0..nf {
            if same[k] && f[k] != first[k] { same[k] = false; }
            if isnum[k] {
                let g = f[k];
                if g.is_empty() || g.len() > 18 || (g.len() > 1 && g[0] == b'0')
                    || !g.iter().all(|c| c.is_ascii_digit()) { isnum[k] = false; }
            }
            if ishex[k] {
                let g = f[k];
                if g.len() != first[k].len()
                    || !g.iter().all(|c| c.is_ascii_digit() || (*c >= b'a' && *c <= b'f')) { ishex[k] = false; }
            }
            if !same[k] && !isnum[k] && !ishex[k] { return None; }
        }
        seen += 1;
    }
    if seen != n { return None; }
    let mut v = Vec::new();
    kst::put(&mut v, nf as u64);
    kst::put(&mut v, sep as u64);
    for k in 0..nf {
        let ty: u64 = if same[k] { 0 } else if isnum[k] { 1 } else { 2 };
        kst::put(&mut v, ty);
        if ty == 0 { kst::put(&mut v, first[k].len() as u64); v.extend_from_slice(first[k]); }
        else if ty == 2 { kst::put(&mut v, first[k].len() as u64); }
    }
    let mut prev = vec![0i64; nf];
    let mut left = n;
    for line in names.split(|&c| c == b'\n') {          // then one to write
        if left == 0 { break; }
        let mut k = 0usize;
        for f in line.split(|&c| c == sep) {
            if !same[k] && isnum[k] {
                let mut x: i64 = 0;
                for &c in f { x = x * 10 + (c - b'0') as i64; }
                let d = x - prev[k]; prev[k] = x;
                kst::put(&mut v, ((d << 1) ^ (d >> 63)) as u64);
            } else if !same[k] {
                let hx = |c: u8| if c <= b'9' { c - b'0' } else { c - b'a' + 10 };
                for j in (0..f.len() - 1).step_by(2) { v.push(hx(f[j]) << 4 | hx(f[j + 1])); }
            }
            k += 1;
        }
        left -= 1;
    }
    let z = ffi::paycomp(&v, ffi::CODEC_BZIP3 | 64)?;
    let mut out = Vec::with_capacity(z.len() + 10);
    kst::put(&mut out, v.len() as u64);
    out.extend_from_slice(&z);
    Some(out)
}

pub fn best_names(names: &[u8], n: u64, level: i32, slow: bool) -> BestNames {
    let t0 = std::time::Instant::now();
    let (z, t) = rayon::join(
        || ffi::zstd(names, level),
        || tok_best(names, slow),
    );
    let i = int_names(names, n).map(|v| ffi::zstd(&v, level));
    let mut out = BestNames { bytes: z, kind: 0, chunk: TCAND[0] };
    if let Some(v) = i {
        if !v.is_empty() && v.len() < out.bytes.len() { out = BestNames { bytes: v, kind: 1, chunk: TCAND[0] }; }
    }
    if let Some((v, c)) = t {
        if !v.is_empty() && v.len() < out.bytes.len() { out = BestNames { bytes: v, kind: 2, chunk: c }; }
    }
    if let Some(v) = fld_names(names, n) {
        if !v.is_empty() && v.len() < out.bytes.len() { out = BestNames { bytes: v, kind: 3, chunk: TCAND[0] }; }
    }
    if crate::prof::ON.load(std::sync::atomic::Ordering::Relaxed) {
        eprintln!("  names {} B: candidates {:.3} s", names.len(), t0.elapsed().as_secs_f64());
    }
    out
}

/// A repeated sequence name would make the index ambiguous, so it is refused.
pub fn check_unique(names: &[u8]) -> Result<(), String> {
    // As the C: hash every name, sort, and compare text only where two
    // hashes collide. A set of 10,000,000 names was a sixth of an encode.
    let lines = names.split(|&c| c == b'\n').filter(|l| !l.is_empty());
    let mut hs: Vec<u64> = lines.clone().map(|l| kst::fnv(FNV0, l)).collect();
    hs.sort_unstable();
    let rep: HashSet<u64> = hs.windows(2).filter(|w| w[0] == w[1]).map(|w| w[0]).collect();
    if rep.is_empty() { return Ok(()); }
    let mut exact: HashSet<&[u8]> = HashSet::new();
    for line in lines {
        if rep.contains(&kst::fnv(FNV0, line)) && !exact.insert(line) {
            let mut e = std::io::stderr().lock();
            let _ = e.write_all(b"karyozip: sequence ").and_then(|_| e.write_all(line)).and_then(|_| e.write_all(b" appears twice\n"));
            return Err("a sequence name repeats".into());
        }
    }
    Ok(())
}
/// One record back out of a spill file.
fn srd_next(f: &mut std::fs::File) -> Result<Option<SRec>, String> {
    let mut h = [0u8; 12];
    match f.read_exact(&mut h) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(_) => return Err("sort temp is corrupt".into()),
    }
    let nrun = u32::from_le_bytes(h[0..4].try_into().unwrap()) as usize;
    let nl = u32::from_le_bytes(h[4..8].try_into().unwrap()) as usize;
    let idx = u32::from_le_bytes(h[8..12].try_into().unwrap());
    if nl >= 512 { return Err("sort temp is corrupt".into()); }
    let mut nm = vec![0u8; nl];
    let mut rb = vec![0u8; nrun * 8];
    f.read_exact(&mut nm).map_err(|_| "sort temp is corrupt".to_string())?;
    f.read_exact(&mut rb).map_err(|_| "sort temp is corrupt".to_string())?;
    let runs = rb.chunks_exact(8)
        .map(|c| (u32::from_le_bytes(c[0..4].try_into().unwrap()),
                  u32::from_le_bytes(c[4..8].try_into().unwrap())))
        .collect();
    Ok(Some(SRec { nm, runs, idx }))
}
