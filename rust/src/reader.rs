//! Reading a KST track. Plain tracks only for now: no --aln, no --db, no
//! sidecar, no permutation. Every refusal message is the C's, because the
//! two are meant to be interchangeable.

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
    let mut v = Vec::with_capacity(n as usize * 8);
    for i in 0..n { v.extend_from_slice(i.to_string().as_bytes()); v.push(b'\n'); }
    v
}

/// --names FILE: a FASTA, a FASTQ, or one name per line. Rank N is record N.
fn read_name_file(path: &str, nreads: u64) -> Vec<u8> {
    let data = rd(std::fs::read(path), "cannot open the file given to --names");
    if data.len() > 1 && data[0] == 0x1f && data[1] == 0x8b {
        crate::bail("the file given to --names is gzipped; decompress it first");
    }
    let lines: Vec<&[u8]> = data.split(|&c| c == b'\n').collect();
    // 1 fastq, 0 fasta, 2 a plain list.
    let fq = if lines.len() >= 3 && lines[0].first() == Some(&b'@') && lines[2].first() == Some(&b'+') { 1 }
             else if lines.len() >= 2 && lines[0].first() == Some(&b'>') && lines[1].first() != Some(&b'>') { 0 }
             else { 2 };
    let mut out = Vec::new();
    let mut got = 0u64;
    for (i, ln) in lines.iter().enumerate() {
        let mut q: &[u8] = ln;
        if fq == 1 { if i % 4 != 0 { continue; } q = &q[1..]; }
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

#[derive(Default)]
pub struct ReadOpts {
    pub raw: bool,
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
    pub names: Vec<Vec<u8>>,
    ibx: Vec<u8>,
    hasidx: bool,
    mhash: u64,
}

fn rd<T>(r: Result<T, std::io::Error>, what: &str) -> T {
    r.unwrap_or_else(|_| crate::bail(what))
}

impl Reader {

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
        // An index trailer sits past the data: take it off, so every size
        // check below still bounds the data alone.
        if fsz >= 32 + 8 + 128 {
            let mut t = [0u8; 32];
            rd(f.seek(SeekFrom::Start(fsz - 32)), "seek");
            rd(f.read_exact(&mut t), "seek");
            let w: Vec<u64> = (0..4)
                .map(|i| u64::from_le_bytes(t[i * 8..i * 8 + 8].try_into().unwrap())).collect();
            if w[3] == IDXMAGIC && w[0] == h[0] && w[0] <= (fsz - 32) / 20
                && fsz - 32 - w[0] * 20 >= 8 + 128 { fsz -= 32 + w[0] * 20; hasidx = true; }
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
            if h[1] & ALNFLAG != 0 {
                crate::bail("this track was written with --aln: decode it with --aln PAF --ref R, or --raw to see what is stored");
            }
            if h[1] & DBFLAG != 0 {
                crate::bail("this track holds featureIDs: decode it with --db features.tsv --out DIR, or --raw to see what is stored");
            }
            if h[1] & LENFLAG != 0 {
                crate::bail("this track leaves out every length the reference implies: decode it with --aln PAF --ref R");
            }
        }
        if h[1] & PERMFLAG != 0 { crate::bail("this build does not read sorted tracks"); }
        rd(f.seek(SeekFrom::Start(8 + 128)), "seek");

        let lzl = h[6] as usize;
        let mut lzraw = vec![0u8; lzl];
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
        let mut nz = vec![0u8; nzl];
        if nzl > 0 { rd(f.read_exact(&mut nz), "short names"); }
        let mut iz = vec![0u8; h[2] as usize];
        rd(f.read_exact(&mut iz), "short index");
        let mut sz = vec![0u8; h[3] as usize];
        rd(f.read_exact(&mut sz), "short sizes");
        let mut itz = vec![0u8; h[11] as usize];
        rd(f.read_exact(&mut itz), "short index table");

        // A sidecar holds the names --names drop:FILE wrote.
        let mut nz = nz;
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
        let sbuf = zdec(&sz, "block sizes");
        let mut b = Vec::new();
        let (mut q, mut acc) = (0usize, 0u64);
        while q < sbuf.len() {
            let (a, n1) = vget(&sbuf, q);
            let (b2, n2) = vget(&sbuf, n1);
            let (c2, n3) = vget(&sbuf, n2);
            let (d2, n4) = vget(&sbuf, n3);
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
        let (mut q, mut acc2) = (0usize, 0u64);
        for _ in 0..h[10] {
            let (cl, n1) = vget(&tbuf, q);
            let (_sp, n2) = vget(&tbuf, n1);
            q = n2;
            acc2 += cl;
            coff.push(acc2);
        }
        if acc2 != h[2] { crate::bail("index chunk table does not match the index"); }
        if h[10] != h[0].div_ceil(ICHUNK) {
            crate::bail("corrupt file: index chunks do not match the sequence count");
        }
        let mut ibx = Vec::new();
        for c in 0..h[10] as usize {
            ibx.extend_from_slice(&zdec(&iz[coff[c] as usize..coff[c + 1] as usize], "index chunk"));
        }

        // The names, by whichever codec stored them.
        let namebuf = if let (Some(f), false) = (&o.names_file, was_side) {
            if !dropped && nzl > 0 {
                crate::bail("this track stores its own names: --names FILE would replace them, so it is refused. Unset it to read the names in the file");
            }
            read_name_file(f, h[0])
        } else if dropped && nzl == 0 {
            if !o.raw && !o.ordnames {
                crate::bail("the names of this track were dropped at encode time: give them back with --names FILE, or pass --names ordinal to accept the ordinals that are stored instead");
            }
            ordinals(h[0])
        } else if nzl == 0 {
            ordinals(h[0])
        } else if names_int {
            int_names_d(&nz, h[0])
        } else if names_tok {
            tok_d(&nz)
        } else {
            zdec(&nz, "names")
        };
        let mut names: Vec<Vec<u8>> = namebuf.split(|&c| c == b'\n').map(|s| s.to_vec()).collect();
        if names.last().map(|x| x.is_empty()) == Some(true) { names.pop(); }
        if names.len() as u64 != h[0] {
            crate::bail("names stream does not match the track's sequence count");
        }

        let data0 = 8 + 128 + lzl as u64 + nml + h[2] + h[3] + h[11];
        Reader {
            f, nreads: h[0], rows: h[4], lab, block: blk, b, data0,
            idcodec: h[9] as u32 as i32, lencodec: h[5] as u32 as i32,
            taildec, tailpred: h[12], tailskip: h[13] & 0xff, names, ibx,
            hasidx, mhash: h[15],
        }
    }

    /// Every sequence, in order. One forward pass over the blocks.
    pub fn decode(&mut self, out: &mut impl std::io::Write) {
        let mut buf: Vec<u8> = Vec::with_capacity(1 << 22);
        let mut c = Cursor { blk: usize::MAX, data: Vec::new(), lp: 0, ip: 0, left: 0 };
        let (mut ipos, mut emitted) = (0usize, 0u64);
        let tskip = self.tailskip;
        for i in 0..self.nreads {
            let (nruns, next) = vget(&self.ibx, ipos);
            ipos = next;
            let name = std::mem::take(&mut self.names[i as usize]);
            let mut st = 0u64;
            for k in 0..nruns {
                let (len, f) = self.step(&mut c, st, k + 1 == nruns && i >= tskip);
                let en = st + len;
                buf.extend_from_slice(&name);
                buf.push(b'\t');
                buf.extend_from_slice(st.to_string().as_bytes());
                buf.push(b'\t');
                buf.extend_from_slice(en.to_string().as_bytes());
                buf.push(b'\t');
                self.label(f, &mut buf);
                buf.push(b'\n');
                if buf.len() > (1 << 22) - 4096 {
                    rd(out.write_all(&buf), "write failed");
                    buf.clear();
                }
                st = en;
                emitted += 1;
            }
            self.names[i as usize] = name;
        }
        rd(out.write_all(&buf), "write failed");
        if emitted != self.rows {
            eprintln!("karyozip: emitted {emitted} rows, the header says {}", self.rows);
            crate::bail("file is truncated");
        }
    }

    /// Decodes one block: its two frames back to back, the split, its rows.
    fn load(&mut self, b: usize) -> (Vec<u8>, usize, u64) {
        if b >= self.b.len() { crate::bail("corrupt index: a block number is past the end"); }
        let m = &self.b[b];
        let mut cb = vec![0u8; m.size as usize];
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
    let mut o = Vec::with_capacity(cnt as usize * 8);
    let (mut q, mut cur) = (0usize, 0u64);
    for _ in 0..cnt {
        let (u, n) = vget(&v, q);
        q = n;
        cur = cur.wrapping_add(unzz(u));
        if cur > 1 << 50 { crate::bail("corrupt names: an integer name past 2^50"); }
        o.extend_from_slice(cur.to_string().as_bytes());
        o.push(b'\n');
    }
    if q != v.len() { crate::bail("corrupt names: the integer stream has more than the count"); }
    o
}

fn tok_d(b: &[u8]) -> Vec<u8> {
    if b.len() < 8 { crate::bail("truncated tokenised names"); }
    let nchunk = u32::from_le_bytes(b[4..8].try_into().unwrap()) as usize;
    let hdr = 8 + 8 * nchunk;
    if nchunk > (b.len() - 8) / 8 || hdr > b.len() {
        crate::bail("corrupt tokenised names header");
    }
    let mut coff = vec![0usize];
    for i in 0..nchunk {
        let cl = u32::from_le_bytes(b[8 + 4 * (nchunk + i)..12 + 4 * (nchunk + i)].try_into().unwrap());
        coff.push(coff[i] + cl as usize);
    }
    if coff[nchunk] > b.len() - hdr {
        crate::bail("corrupt tokenised names: chunks overrun the stream");
    }
    let mut o = Vec::new();
    for i in 0..nchunk {
        let s = &b[hdr + coff[i]..hdr + coff[i + 1]];
        let mut d = ffi::tok_decode(s)
            .unwrap_or_else(|| crate::bail("name tokeniser decode failed"));
        // The tokeniser separates with NUL; the format stores newlines.
        for c in d.iter_mut() { if *c == 0 { *c = b'\n'; } }
        o.append(&mut d);
    }
    o
}

pub fn vget(b: &[u8], mut i: usize) -> (u64, usize) {
    let (mut v, mut s) = (0u64, 0u32);
    loop {
        if i >= b.len() || s > 63 { crate::bail("corrupt varint"); }
        let c = b[i];
        i += 1;
        v |= ((c & 0x7f) as u64) << s;
        s += 7;
        if c < 0x80 { return (v, i); }
    }
}

fn unzz(v: u64) -> u64 { ((v >> 1) as i64 ^ -((v & 1) as i64)) as u64 }

/* ------------------------------------------------ one sequence at a time */

/// A position in the row stream: the block, its bytes, a cursor per frame,
/// and how many of its rows are left.
struct Cursor { blk: usize, data: Vec<u8>, lp: usize, ip: usize, left: u64 }

impl Reader {
    /// (block, offset in it, run count) for sequence `ri`.
    fn locate(&self, ri: u64) -> (usize, u64, u64) {
        let (mut q, mut cpend) = (0usize, 0u64);
        let mut nruns = 0u64;
        for _ in 0..=ri {
            let (n, next) = vget(&self.ibx, q);
            q = next;
            if cpend / self.block > u32::MAX as u64 {
                crate::bail("corrupt index: a block number past 32 bits");
            }
            nruns = n;
            cpend += n;
        }
        let start = cpend - nruns;
        (( start / self.block) as usize, start % self.block, nruns)
    }

    fn at(&mut self, blk: usize, off: u64) -> Cursor {
        let (data, split, rows) = self.load(blk);
        let mut c = Cursor { blk, data, lp: 0, ip: split, left: rows };
        for _ in 0..off { c.lp = skip(&c.data, c.lp); c.ip = skip(&c.data, c.ip); }
        c.left -= off;
        c
    }

    /// Next row: (length, label id). Steps to the next block when this one ends.
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
        } else {
            if v == 0 { crate::bail("corrupt track: a zero length outside a sentinel block"); }
            v
        };
        let (f, i2) = vget(&c.data, c.ip);
        c.ip = i2;
        c.left -= 1;
        (len, f)
    }

    fn label(&self, f: u64, out: &mut Vec<u8>) {
        match self.lab.get(f as usize) {
            Some(l) => out.extend_from_slice(l),
            None if self.lab.is_empty() => out.extend_from_slice(f.to_string().as_bytes()),
            None => crate::bail("corrupt block: a label id past the table"),
        }
    }

    fn find(&self, seq: &[u8]) -> i64 {
        self.names.iter().position(|n| n == seq).map_or(-1, |i| i as i64)
    }

    /// region F SEQ | SEQ POS | SEQ START END [HALO].
    pub fn region(&mut self, seq: &[u8], lo: u64, hi: u64, halo: u64,
                  out: &mut impl std::io::Write) {
        let ri = self.find(seq);
        if ri < 0 { crate::bail("sequence not found"); }
        let ri = ri as u64;
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
                buf.extend_from_slice(st.to_string().as_bytes());
                buf.push(b'\t');
                buf.extend_from_slice(en.to_string().as_bytes());
                buf.push(b'\t');
                self.label(f, &mut buf);
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
            .map(|i| (kst::fnv(FNV0, &self.names[i]), i as u32)).collect();
        pair.sort_unstable();
        let dup = pair.windows(2).filter(|w| w[0].0 == w[1].0).count() as u64;
        let mut t = Vec::with_capacity(n * 20 + 32);
        for (k, _) in &pair { t.extend_from_slice(&k.to_le_bytes()); }
        for (_, o) in &pair {
            t.extend_from_slice(&(h64b(&self.names[*o as usize])
                                  ^ (*o as u64).wrapping_mul(H64B_K)).to_le_bytes());
        }
        for (_, o) in &pair { t.extend_from_slice(&o.to_le_bytes()); }
        for v in [self.nreads, dup, self.mhash, IDXMAGIC] {
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
