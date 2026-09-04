//! The KST writer. Rows arrive grouped by sequence, gapless, each sequence
//! starting at 0; blocks fill and are staged; the file is assembled at close.
//! Every layout decision here is c/karyozip.c's, because the two must produce
//! the same bytes.

use crate::codec;
use crate::ffi;
use crate::names;
use crate::varint::{self, FNV0};
use rayon::prelude::*;
use std::collections::HashMap;
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
            varint::put(&mut self.lens, self.plen);
            self.phas = false;
            return;
        }
        if self.tailon != Some(true) {
            varint::put(&mut self.lens, self.plen);
            self.phas = false;
            return;
        }
        varint::put(&mut self.lens, varint::zz(self.curbp as i64 - self.tailpred as i64));
        self.phas = false;
    }

    fn pushread(&mut self) {
        if self.cur.is_empty() { return; }
        varint::put(&mut self.ibx, self.currun);
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
            varint::put(&mut self.lens, self.plen);
            self.phas = false;
        }
        if self.pend >= self.opt.block { self.flush().map_err(io)?; }
        if self.pend == 0 { self.blkbp0 = self.curbp; }
        self.plen = run_len;
        self.phas = true;
        varint::put(&mut self.ids, id);
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
        let p = codec::pick(&sref, seekable, level, self.opt.id_codec, self.opt.len_codec);

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
                    let a = codec::comp_lens(&lens, p.len, level).expect("lens frame");
                    let b = codec::comp_ids(&ids, p.id, level).expect("ids frame");
                    (a, b, lens.len(), ids.len(), self.bm[i].rows, self.bm[i].bp)
                })
                .collect();
            for (a, b, nl, ni, rows, bp) in &done {
                body.write_all(a).map_err(io)?;
                body.write_all(b).map_err(io)?;
                byt_lens += a.len() as u64;
                byt_labs += b.len() as u64;
                varint::put(&mut sizes, a.len() as u64);
                varint::put(&mut sizes, b.len() as u64);
                varint::put(&mut sizes, *ni as u64);
                varint::put(&mut sizes, *nl as u64);
                varint::put(&mut sizes, *rows);
                varint::put(&mut sizes, *bp);
                varint::put(&mut sizes, varint::fnv(varint::fnv(FNV0, a), b));
            }
        }

        let mut labs = Vec::new();
        for l in &self.labtab { labs.extend_from_slice(l); labs.push(b'\n'); }

        let mut ordinal = self.ordok && self.nidx > 0;
        if self.opt.keep_names { ordinal = false; }
        if self.nidx > 1 { names::check_unique(&self.names)?; }
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
            let b = names::best(&self.names, self.nidx, level);
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
                    varint::put(&mut itab, cz.len() as u64);
                    varint::put(&mut itab, spend);
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
        let mut mh = varint::fnv(FNV0, &hb);
        mh = varint::fnv(mh, &lz);
        mh = varint::fnv(mh, &nz);
        mh = varint::fnv(mh, &ichunks);
        mh = varint::fnv(mh, &sz);
        mh = varint::fnv(mh, &itz);
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
