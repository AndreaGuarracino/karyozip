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

/// One block: where it is, both frame sizes, both decoded sizes, its rows,
/// the bases before it, and its own hash.
struct Blk { off: u64, size: u64, l: u64, raw: u64, lraw: u64, rows: u64, hash: u64 }

pub struct Reader {
    f: std::fs::File,
    pub nreads: u64,
    pub rows: u64,
    lab: Vec<Vec<u8>>,
    b: Vec<Blk>,
    data0: u64,
    idcodec: i32,
    lencodec: i32,
    taildec: bool,
    tailpred: u64,
    tailskip: u64,
    pub names: Vec<Vec<u8>>,
    ibx: Vec<u8>,
}

fn rd<T>(r: Result<T, std::io::Error>, what: &str) -> T {
    r.unwrap_or_else(|_| crate::bail(what))
}

impl Reader {
    pub fn open(path: &str) -> Reader {
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
        // An index trailer sits past the data: take it off, so every size
        // check below still bounds the data alone.
        if fsz >= 32 + 8 + 128 {
            let mut t = [0u8; 32];
            rd(f.seek(SeekFrom::Start(fsz - 32)), "seek");
            rd(f.read_exact(&mut t), "seek");
            let w: Vec<u64> = (0..4)
                .map(|i| u64::from_le_bytes(t[i * 8..i * 8 + 8].try_into().unwrap())).collect();
            if w[3] == IDXMAGIC && w[0] == h[0] && w[0] <= (fsz - 32) / 20
                && fsz - 32 - w[0] * 20 >= 8 + 128 { fsz -= 32 + w[0] * 20; }
        }
        let nml = h[1] & !FLAGS;
        if nml > fsz || h[2] > fsz || h[3] > fsz || h[6] > fsz || h[11] > fsz
            || h[14] > fsz || h[11] + h[14] > fsz
            || 8 + 128 + nml + h[2] + h[3] + h[6] + h[11] + h[14] > fsz {
            crate::bail("corrupt header: the stream sizes do not fit the file");
        }
        if h[1] & (ALNFLAG | DBFLAG | SIDEFLAG | LENFLAG | PERMFLAG) != 0 {
            crate::bail("this build reads plain tracks only");
        }
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

        // The metadata hash covers the header with h[15] zeroed, then every
        // stream. A foreign or damaged file is caught here, not at a block.
        {
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
            if mh != h[15] { crate::bail("corrupt file: the metadata does not match its hash"); }
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
        let namebuf = if nzl == 0 {
            let mut v = Vec::with_capacity(h[0] as usize * 8);
            for i in 0..h[0] { v.extend_from_slice(i.to_string().as_bytes()); v.push(b'\n'); }
            v
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
            f, nreads: h[0], rows: h[4], lab, b, data0,
            idcodec: h[9] as u32 as i32, lencodec: h[5] as u32 as i32,
            taildec, tailpred: h[12], tailskip: h[13] & 0xff, names, ibx,
        }
    }

    /// The runs of every sequence, in order. One forward pass: a cursor per
    /// frame, advanced row by row, so a block is walked once however many
    /// sequences it holds. Restarting at the block head per sequence is what
    /// made the C's smooth 9x slower before it was fixed.
    pub fn decode(&mut self, out: &mut impl std::io::Write) {
        let mut emitted = 0u64;
        let mut buf: Vec<u8> = Vec::with_capacity(1 << 22);
        let mut blk = usize::MAX;
        let (mut data, mut lp, mut ip, mut left) = (Vec::new(), 0usize, 0usize, 0u64);
        let mut ipos = 0usize;
        for i in 0..self.nreads {
            let (nruns, next) = vget(&self.ibx, ipos);
            ipos = next;
            let name = std::mem::take(&mut self.names[i as usize]);
            let mut st = 0u64;
            for k in 0..nruns {
                if left == 0 {
                    blk = if blk == usize::MAX { 0 } else { blk + 1 };
                    let (d, split, rows) = self.load(blk);
                    data = d; lp = 0; ip = split; left = rows;
                }
                let last = k + 1 == nruns && i >= self.tailskip;
                let (v, l2) = vget(&data, lp);
                lp = l2;
                let len = if last && self.taildec {
                    let t = self.tailpred.wrapping_add(unzz(v));
                    if t <= st { crate::bail("corrupt track: sequence total is not past its last run"); }
                    t - st
                } else {
                    if v == 0 { crate::bail("corrupt track: a zero length outside a sentinel block"); }
                    v
                };
                let (f, i2) = vget(&data, ip);
                ip = i2;
                let en = st + len;
                buf.extend_from_slice(&name);
                buf.push(b'\t');
                buf.extend_from_slice(st.to_string().as_bytes());
                buf.push(b'\t');
                buf.extend_from_slice(en.to_string().as_bytes());
                buf.push(b'\t');
                match self.lab.get(f as usize) {
                    Some(l) => buf.extend_from_slice(l),
                    None if self.lab.is_empty() => buf.extend_from_slice(f.to_string().as_bytes()),
                    None => crate::bail("corrupt block: a label id past the table"),
                }
                buf.push(b'\n');
                if buf.len() > (1 << 22) - 4096 {
                    rd(out.write_all(&buf), "write failed");
                    buf.clear();
                }
                st = en;
                emitted += 1;
                left -= 1;
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
