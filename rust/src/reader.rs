//! Reading a KST track. Every refusal message is the C's, because the two
//! are meant to be interchangeable.

use crate::ffi;
use crate::kst::{self, FNV0, TOKFLAG, INTFLAG, FLDFLAG, DROPPED, ALNFLAG, DBFLAG, SIDEFLAG, LENFLAG,
                 TAILFLAG, RUNFLAG, ILVFLAG, GRPFLAG, PERMFLAG, PERMTOK, CARTFLAG, ICHUNK, FLAGS};
use std::io::{Read, Seek, SeekFrom};

const IDXMAGIC: u64 = 0x3130_5844_4954_534b;   // "KSTIDX01"

/// A decoded block: its number, bytes, where the labels start, rows.
type Cached = (usize, std::sync::Arc<Vec<u8>>, usize, u64);
pub type Shared = (Cached, std::sync::Arc<Vec<u32>>);

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
#[derive(Clone, Copy)]
struct Blk { off: u64, size: u64, l: u64, raw: u64, lraw: u64, rows: u64, bp: u64, hash: u64 }

/// seq, start, end, label, as decode hands each row out.
/// name, start, end, label, and the label's table index when the label is that entry's own text
pub type RowFn<'a> = dyn FnMut(&[u8], u64, u64, &[u8], Option<usize>) + 'a;

#[derive(Default)]
pub struct ReadOpts {
    pub raw: bool,
    /// the index command: names only, so no companion is asked for
    pub index: bool,
    pub db: bool,
    pub aln: Option<std::sync::Arc<crate::aln::Aln>>,
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
    names: std::cell::OnceCell<std::sync::Arc<(Vec<u8>, Vec<u64>)>>,
    nz: Vec<u8>,
    pz: Vec<u8>,
    names_file: Option<String>,
    /// side, perm_ext, dropped, int, tok, perm, perm_tok
    nf: u8,
    /// the appended index: offset of its key table, and its trailer's binding word
    idx: Option<(u64, u64, u64, u64, std::fs::File)>,   // index base, binding, bucket bits, payload, TRACK.kzi
    /// name -> rank, built on the first lookup, which is also where the C
    /// notices that two sequences share a name.
    /// as the C's buildhash: a power-of-two table of ordinals over the
    /// decoded names, open addressed, no key copied
    nameix: std::cell::OnceCell<(Vec<u32>, u64)>,
    /// the sequence index: chunk frames, their offsets, the run count each
    /// starts at, and the one chunk decoded now, as the C's loadichunk
    iz: Vec<u8>, coff: Vec<u64>, cpend: Vec<u64>, ibx: Vec<u8>, cur_chunk: u64,
    runhead: bool, ilv: bool, grp: bool, rtab: std::sync::Arc<Vec<u32>>,   // run-count prefix sums, lifted out of the lengths; 32 bits, as the C's offsets
    mhash: u64,
    cur_ri: u64,
    cur_q: usize,
    cur_pend: u64,
    cached: Option<Cached>,
    /// the rows as text, what binds this track when it serves as --ref
    pub rowhash: u64,
    /// not --aln, --db, --sort, --names drop, names in the file: fit for --ref
    pub plainref: bool,
    /// where the last seq_rows stopped: (block, rows into it, cursor), so a
    /// walk in order continues instead of re-skipping the block's rows
    sc: Option<(usize, u64, Cursor)>,
    feathash: u64,
    cart: bool,
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
    pub fn cart(&self) -> bool { self.cart }

    pub fn open_with(path: &str, o: &ReadOpts) -> Reader {
        let mut f = rd(std::fs::File::open(path), "open");
        let mut m = [0u8; 8];
        rd(f.read_exact(&mut m), "not a KST file: shorter than the magic");
        if &m != crate::writer::MAGIC {
            let sh: String = m.iter()
                .map(|&c| if (32..127).contains(&c) { c as char } else { '.' }).collect();
            eprintln!("karyozip: file starts with \"{sh}\", not \"KST00001\"");
            crate::bail("not a KST file");
        }
        let mut hb = [0u8; 128];
        rd(f.read_exact(&mut hb), "short header");
        let mut h = [0u64; 16];
        for (i, v) in h.iter_mut().enumerate() {
            *v = u64::from_le_bytes(hb[i * 8..i * 8 + 8].try_into().unwrap());
        }
        let h1_raw = h[1];
        let fsz = rd(f.seek(SeekFrom::End(0)), "seek");
        // the name index, TRACK.kzi beside the track, as the C reads it
        let mut idx = None;
        if let Ok(mut xf) = std::fs::File::open(format!("{path}.kzi")) {
            let xsz = xf.seek(SeekFrom::End(0)).unwrap_or(0);
            let mut t = [0u8; 40];
            if xsz >= 40 && xf.seek(SeekFrom::Start(xsz - 40)).is_ok() && xf.read_exact(&mut t).is_ok() {
                let w: Vec<u64> = (0..5).map(|i| u64::from_le_bytes(t[i * 8..i * 8 + 8].try_into().unwrap())).collect();
                if w[4] == IDXMAGIC && w[0] == h[0] && w[1] < 40 && w[3] <= xsz - 40 {
                    idx = Some((xsz - 40 - w[3], w[2], w[1], w[3], xf));
                }
            }
        }
        let nml = h[1] & !FLAGS;
        if nml > fsz || h[2] > fsz || h[3] > fsz || h[6] > fsz || h[11] > fsz
            || h[14] > fsz || h[11] + h[14] > fsz
            || 8 + 128 + nml + h[2] + h[3] + h[6] + h[11] + h[14] > fsz {
            crate::bail("corrupt header: the stream sizes do not fit the file");
        }
        rd(f.seek(SeekFrom::Start(8 + 128)), "seek");

        let lzl = h[6] as usize;
        let mut lzraw = vec![0u8; crate::ckalloc(lzl)];
        if lzl > 0 { rd(f.read_exact(&mut lzraw), "short labels"); }
        let mut lab: Vec<Vec<u8>> = Vec::new();
        if lzl > 0 {
            let lbuf = zdec(&lzraw, "label table");
            if h[7] > lbuf.len() as u64 { crate::bail("label count exceeds the label table"); }
            if h[7] == 0 { crate::bail("corrupt file: a label table with no labels"); }
            let mut q = 0;
            for _ in 0..h[7] {
                let nl = crate::ffi::find(b'\n', &lbuf[q..]).unwrap_or_else(|| crate::bail("bad label table"));
                lab.push(lbuf[q..q + nl].to_vec()); q += nl + 1;
            }
            if q != lbuf.len() { crate::bail("corrupt file: the label table has more lines than its count"); }
        } else if h[7] != 0 {
            crate::bail("corrupt file: a label count with no label table");
        }

        let names_tok = h[1] & TOKFLAG != 0;
        let names_int = h[1] & INTFLAG != 0;
        let names_fld = h[1] & FLDFLAG != 0;
        if names_tok as u8 + names_int as u8 + names_fld as u8 > 1 {
            crate::bail("corrupt file: names flagged more than one way");
        }
        let was_side = h[1] & SIDEFLAG != 0;
        let dropped = h[1] & DROPPED != 0;
        if !o.raw && !o.index {
            if h[1] & ALNFLAG != 0 && o.aln.is_none() {
                crate::bail("this track was written with --aln: read it with --aln PAF --ref R, or --raw to see what is stored");
            }
            if h[1] & ALNFLAG == 0 && o.aln.is_some() {
                crate::bail("this track was not written with --aln");
            }
            if h[1] & DBFLAG != 0 && !o.db {
                crate::bail("this track holds featureIDs: read it with --db features.tsv --out-dir DIR, or --raw to see what is stored");
            }
            if h[1] & DBFLAG == 0 && o.db {
                crate::bail("this track does not hold featureIDs");
            }
            if h[1] & LENFLAG != 0 && o.aln.is_none() {
                crate::bail("this track leaves out every length the reference implies: read it with --aln PAF --ref R");
            }
            if was_side && o.names_file.is_none()
               && std::fs::metadata(format!("{path}.names")).is_err() {
                crate::bail("the names of this track are in a sidecar: give it with --names FILE, or --raw to see the ordinals");
            }
        }
        if let Some(a) = &o.aln {
            if h[1] & ALNFLAG != 0 && (h[13] >> 8) != (a.hash & ((1u64 << 56) - 1)) {
                crate::bail("this PAF, reference or name map is not the one the track was encoded with");
            }
        }
        let was_perm = h[1] & PERMFLAG != 0;
        let perm_tok = h[1] & PERMTOK != 0;
        let taildec = h[1] & TAILFLAG != 0;
        let runhead = h[1] & RUNFLAG != 0;
        let ilv = h[1] & ILVFLAG != 0;
        let grp = h[1] & GRPFLAG != 0;
        if grp && (ilv || !runhead) { crate::bail("corrupt file: grouped lengths without their run counts"); }
        if runhead {   // the C allocates its u32 run table here
            crate::ckalloc(h[0].max(1).checked_mul(4).unwrap_or_else(|| crate::bail("allocation size overflows")) as usize);
        }
        if h[1] & CARTFLAG != 0 && h[1] & DBFLAG == 0 {
            crate::bail("corrupt file: vocabulary featureIDs without --db");
        }
        if taildec && h[13] & 0xff != 64 {
            crate::bail("corrupt file: tail skip is not the format constant");
        }
        let nzl = nml as usize;
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
            let (bp, n6) = vget(&sbuf, n5);
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
            b.push(Blk { off: acc, size: a + b2, l: a, raw: c2, lraw: d2, rows: e2, bp, hash: g2 });
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
            if cl > h[2] - acc2 { crate::bail("corrupt file: an index chunk runs past the index"); }   // lengths that wrap would pass the total
            acc2 += cl;
            coff.push(acc2);
            cpend.push(sp);
        }
        cpend.push(h[4]);
        if acc2 != h[2] { crate::bail("index chunk table does not match the index"); }
        if !runhead && h[10] != h[0].div_ceil(ICHUNK) {
            crate::bail("corrupt file: index chunks do not match the sequence count");
        }
        // as the C: a start out of order, or a first start past 0, shifts every row
        if cpend.windows(2).any(|w| w[0] > w[1]) || (h[10] > 0 && cpend[0] != 0) {
            crate::bail("corrupt file: index chunk starts are not in order");
        }

        // A permutation track with no stored blob keeps its names outside,
        // so --names FILE is the source there, not a replacement.
        let perm_ext = was_perm && pz.is_empty();

        // The names are decoded on first use: a query through the index
        // never needs them, and a short-read track has 20 million.
        let mut nf = 0u8;
        for (bit, on) in [(1, was_side), (2, perm_ext), (4, dropped), (8, names_int), (16, names_tok), (32, was_perm), (64, perm_tok), (128, names_fld)] {
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
            taildec, tailpred: h[12], tailskip: h[13] & 0xff, iz, coff, cpend, ibx: Vec::new(), cur_chunk: u64::MAX,
            runhead, ilv, grp, rtab: Default::default(),
            names: std::cell::OnceCell::new(), nz, pz, names_file: o.names_file.clone(),
            nf, idx,
            mhash: h[15], cur_ri: 0, cur_q: 0, cur_pend: 0, cached: None, sc: None,
            rowhash: if h1_raw & ALNFLAG != 0 { 0 } else { h[13] >> 8 },
            plainref: h1_raw & (ALNFLAG | DBFLAG | SIDEFLAG | DROPPED | PERMFLAG) == 0,
            feathash: (h[5] >> 32) | ((h[9] >> 32) << 32),
            cart: h[1] & CARTFLAG != 0,
            aln: o.aln.clone(), lenfile: h[1] & LENFLAG != 0 && !o.raw, raw: o.raw,
            alnq: None, blk_nos,
        }
    }

    /// Every row, in order, to a callback. decode and --db share it.
    pub fn for_each_row(&mut self, f: &mut RowFn) {
        let mut c = Cursor { blk: usize::MAX, data: std::sync::Arc::new(Vec::new()),
                             lp: 0, ip: 0, left: 0 };
        let mut emitted = 0u64;
        let tskip = self.tailskip;
        let (mut name, mut lab) = (Vec::new(), Vec::new());
        for i in 0..self.nreads {
            let nruns = self.locate(i).2;
            name.clear(); name.extend_from_slice(self.name_at(i));
            self.alnq = self.aln.as_ref().and_then(|a| a.find(&name));
            let mut st = 0u64;
            for k in 0..nruns {
                let (len, id) = self.step(&mut c, st, k + 1 == nruns && i >= tskip);
                lab.clear();
                self.label(id, &mut lab, st, len);
                let own = self.lab.get(id as usize).is_some_and(|l| self.aln.is_none() || l != b"=");
                f(&name, st, st + len, &lab, own.then_some(id as usize));
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
        let mut ocap = 1usize << 22;
        let mut buf: Vec<u8> = Vec::with_capacity(ocap);
        let mut c = Cursor { blk: usize::MAX, data: std::sync::Arc::new(Vec::new()), lp: 0, ip: 0, left: 0 };
        let (mut emitted, mut restored) = (0u64, 0u64);
        let tskip = self.tailskip;
        let mut name = Vec::new();
        for i in 0..self.nreads {
            let nruns = self.locate(i).2;
            name.clear(); name.extend_from_slice(self.name_at(i));
            self.alnq = self.aln.as_ref().and_then(|a| a.find(&name));
            let mut st = 0u64;
            for k in 0..nruns {
                let (len, f) = self.step(&mut c, st, k + 1 == nruns && i >= tskip);
                let en = st + len;
                if self.aln.is_some() && self.lab.get(f as usize).is_some_and(|l| l == b"=") { restored += 1; }
                let lab = self.label_at2(f, st, len);
                // the C's flush point, before a row that would not fit: a
                // damaged file then stops at the same byte in both
                let need = name.len() + lab.map_or(20, |l| l.len()) + 48;
                if buf.len() + need > ocap {
                    if !buf.is_empty() { rd(out.write_all(&buf), "write failed"); buf.clear(); }
                    ocap = ocap.max(need);
                }
                match lab {
                    Some(l) => kst::put_row(&mut buf, &name, st, en, l),
                    None => {
                        buf.extend_from_slice(&name);
                        buf.push(b'\t'); kst::putu64(&mut buf, st);
                        buf.push(b'\t'); kst::putu64(&mut buf, en);
                        buf.push(b'\t'); kst::putu64(&mut buf, f); buf.push(b'\n');
                    }
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
    /// A one-block track's decoded block, for another reader of the same file.
    pub fn share_block(&mut self) -> Option<Shared> {
        if self.b.len() != 1 { return None; }
        self.load(0);
        Some((self.cached.clone().unwrap(), self.rtab.clone()))
    }
    pub fn adopt_block(&mut self, s: Shared) { self.cached = Some(s.0); self.rtab = s.1; }

    fn load(&mut self, b: usize) -> (std::sync::Arc<Vec<u8>>, usize, u64) {
        if let Some((cb, d, sp, r)) = &self.cached {
            if *cb == b { return (d.clone(), *sp, *r); }
        }
        let (d, sp, r) = self.load_raw(b);
        let d = std::sync::Arc::new(d);
        if let Some((_, old, _, _)) = self.cached.replace((b, d.clone(), sp, r)) {
            if let Ok(v) = std::sync::Arc::try_unwrap(old) { ffi::recycle(v); }   // the cursor let it go first
        }
        (d, sp, r)
    }

    /// One stream a read, its run count then a label beside its own length.
    /// Split it back into the two regions the walks expect; the counts stay
    /// in front of the lengths for rmheads to lift out.
    fn ilv_split(&self, v: Vec<u8>, nl: usize, ni: usize) -> Vec<u8> {
        let (mut lens, mut ids) = (Vec::with_capacity(nl + ni), Vec::with_capacity(ni));   // the ids join lens: no regrowth
        let (mut i, mut seen) = (0usize, 0u64);
        while i < v.len() {
            let h0 = i;
            let (k, n) = vget(&v, i); i = n;
            if k == 0 || k > self.rows - seen { crate::bail("corrupt block: a run count is out of range"); }
            lens.extend_from_slice(&v[h0..i]);
            for _ in 0..k {
                if i >= v.len() { crate::bail("corrupt block: the interleaved frame is short"); }
                let (lb, n2) = vget(&v, i);
                kst::put(&mut ids, crate::writer::ilv_flip(lb)); i = n2;
                if i >= v.len() { crate::bail("corrupt block: the interleaved frame is short"); }
                let n2 = skip(&v, i); lens.extend_from_slice(&v[i..n2]); i = n2;
            }
            seen += k;
        }
        if seen != self.rows || lens.len() != nl || ids.len() != ni {
            crate::bail("corrupt block: the interleaved frame does not split into its two regions");
        }
        lens.extend_from_slice(&ids);
        lens
    }

    /// The run count heads every read on disk. Lift them into rtab and close
    /// the gap, so every walk below sees the layout it always saw.
    /// In place, as the C's memmove: a second copy of the lengths was the
    /// largest buffer of an assembly decode.
    fn rmheads(&mut self, lens: &mut [u8]) -> usize {
        let (mut w, mut i, mut seen) = (0usize, 0usize, 0u64);
        // 4 B a read, as the C's counts: 8 B was 295 MB more on 73.7M reads
        if self.rows > u32::MAX as u64 { crate::bail("corrupt index: a sequence offset past 32 bits"); }
        let mut tab: Vec<u32> = Vec::with_capacity(self.nreads as usize + 1);
        tab.push(0);
        while i < lens.len() {
            if tab.len() as u64 > self.nreads { crate::bail("corrupt block: more reads than the header claims"); }
            let (k, n) = vget(&lens, i); i = n;
            if k == 0 || k > self.rows - seen { crate::bail("corrupt block: a run count is out of range"); }
            seen += k; tab.push(seen as u32);
            let s0 = i;   // the read's lengths sit together: one move a read
            for _ in 0..k {
                if i >= lens.len() { crate::bail("corrupt block: a read runs past its lengths"); }
                i = skip(&lens, i);
                if i > lens.len() { crate::bail("corrupt block: a length runs past the lengths"); }
            }
            lens.copy_within(s0..i, w); w += i - s0;
        }
        if tab.len() as u64 != self.nreads + 1 || seen != self.rows {
            crate::bail("corrupt block: the run counts do not add up to the row count");
        }
        self.rtab = std::sync::Arc::new(tab);
        w
    }

    /// Undo grp_pack: the labels, still in their own order, say which group
    /// each length came from. Runs after rmheads, so rtab holds the counts.
    fn grp_unpack(&self, lens: &[u8], ids: &[u8]) -> Vec<u8> {
        let mut out = Vec::with_capacity(lens.len() + ids.len());   // the ids join it next: no regrowth
        let (mut p, mut q) = (0usize, 0usize);
        // 8 B a run in each, as the C's: positions as pairs were 16 B in two
        let (mut kv, mut gv, mut lv, mut cnt): (Vec<u64>, Vec<u64>, Vec<u64>, Vec<u64>) = Default::default();
        for w in self.rtab.windows(2) {
            let k = (w[1] - w[0]) as u64;
            let mut top = 0u64;
            kv.clear(); gv.clear();
            for j in 0..k {
                if p >= lens.len() || q >= ids.len() { crate::bail("corrupt block: a grouped read runs past its frames"); }
                let (v, n2) = vget(lens, p); gv.push(v); p = n2;
                let (id, n3) = vget(ids, q); q = n3;
                if p > lens.len() || q > ids.len() { crate::bail("corrupt block: a grouped value runs past its frame"); }
                if id >> 32 != 0 { crate::bail("corrupt block: a grouped label is out of range"); }
                kv.push((id << 32) | j); top = top.max(id);
            }
            lv.clear(); lv.resize(k as usize, 0);
            // label order, stable in j, as the C's: a few runs sort, many count
            if k <= 32 || top > 4 * k + 1024 {
                kv.sort_unstable();
                for (j, key) in kv.iter().enumerate() { lv[(key & 0xffff_ffff) as usize] = gv[j]; }
            } else {
                cnt.clear(); cnt.resize(top as usize + 2, 0);
                for key in &kv { cnt[(key >> 32) as usize + 1] += 1; }
                for c in 1..=top as usize { cnt[c] += cnt[c - 1]; }
                for (j, key) in kv.iter().enumerate() {
                    let g = &mut cnt[(key >> 32) as usize];
                    lv[j] = gv[*g as usize]; *g += 1;   // j's length is the next of its label
                }
            }
            for &v in lv.iter() { kst::put(&mut out, v); }
        }
        if p != lens.len() || q != ids.len() { crate::bail("corrupt block: the grouped frames do not line up"); }
        if out.len() != lens.len() { crate::bail("corrupt block: the grouped lengths change size"); }   // as the C: re-spelled, not copied
        out
    }

    fn load_raw(&mut self, b: usize) -> (Vec<u8>, usize, u64) {
        if b >= self.b.len() { crate::bail("corrupt index: a block number is past the end"); }
        let m = self.b[b].clone();
        let mut cb = ffi::raw_buf(crate::ckalloc(m.size as usize));
        cb.resize(m.size as usize, 0);
        rd(self.f.seek(SeekFrom::Start(self.data0 + m.off)), "seek");
        rd(self.f.read_exact(&mut cb), "file is truncated: a block is missing or short");
        if kst::fnv(FNV0, &cb) != m.hash {
            crate::bail("corrupt block: its bytes do not match the stored hash");
        }
        let (la, lb) = (m.l as usize, (m.size - m.l) as usize);
        let (lraw, raw) = (m.lraw as usize, m.raw as usize);
        // both frames into one buffer, the labels behind the lengths, as the C's blkbuf
        let mut data = ffi::raw_buf(crate::ckalloc(lraw + raw));
        frame(&cb[..la], self.lencodec, &mut data, if self.ilv { lraw + raw } else { lraw });
        if !self.ilv { frame(&cb[la..la + lb], self.idcodec, &mut data, raw); }
        // as the C frees it: not needed past its frames; kept for the next
        // block only when there is one, or a one-block decode holds it to the end
        if self.b.len() > 1 { ffi::recycle(cb); } else { drop(cb); }
        let mut data = if self.ilv { self.ilv_split(data, lraw, raw) } else { data };
        let mut split = lraw;
        if self.runhead {   // close the gap the run counts leave, as the C's memmove
            split = self.rmheads(&mut data[..lraw]);
            data.copy_within(lraw.., split);
            data.truncate(split + raw); data.shrink_to_fit();   // the heads' touched tail goes back
        }
        if self.grp {
            let out = self.grp_unpack(&data[..split], &data[split..]);
            data[..split].copy_from_slice(&out);
        }
        (data, split, m.rows)
    }
}

fn frame(src: &[u8], codec: i32, out: &mut Vec<u8>, raw: usize) {
    let ok = if codec != 0 { ffi::paydecomp_into(src, codec, out, raw) } else { ffi::zstd_decompress_into(src, out, raw) };
    if !ok { crate::bail("frame decode failed"); }
}

fn zdec(src: &[u8], what: &str) -> Vec<u8> {
    let n = ffi::zstd_frame_size(src)
        .unwrap_or_else(|| crate::bail(&format!("{what} has no usable frame size")));
    ffi::zstd_decompress(src, n)
        .unwrap_or_else(|| crate::bail(&format!("{what} decompress failed")))
}

/// The field-delta names, as the C's fldnamesd rebuilds them.
fn fld_names_d(b: &[u8], cnt: u64) -> Vec<u8> {
    let (raw, hp) = vget(b, 0);
    if raw == 0 || raw > 1 << 40 { crate::bail("corrupt names: a bad field blob size"); }
    let v = match crate::ffi::paydecomp(&b[hp..], crate::ffi::CODEC_BZIP3 | 64, raw as usize) {
        Some(v) => v, None => crate::bail("field names decompress failed"),
    };
    let (nf, mut q) = vget(&v, 0);
    if nf == 0 || nf > 16 { crate::bail("corrupt names: a bad field count"); }
    let (sep, n1) = vget(&v, q); q = n1;
    let sep = sep as u8;
    let (mut cf, mut var, mut hw): (Vec<&[u8]>, Vec<u64>, Vec<usize>) = (Vec::new(), Vec::new(), Vec::new());
    for _ in 0..nf {
        let (t, n2) = vget(&v, q); q = n2;
        if t > 2 { crate::bail("corrupt names: a bad field type"); }
        var.push(t);
        if t == 1 { cf.push(&[]); hw.push(0); }
        else if t == 2 {
            let (l, n3) = vget(&v, q); q = n3;
            if l == 0 || l % 2 == 1 || l > 64 { crate::bail("corrupt names: a bad hex field width"); }
            cf.push(&[]); hw.push(l as usize);
        } else {
            let (l, n3) = vget(&v, q); q = n3;
            if q + l as usize > v.len() { crate::bail("corrupt names: a field runs past the blob"); }
            cf.push(&v[q..q + l as usize]); hw.push(0); q += l as usize;
        }
    }
    let nvar = var.iter().filter(|&&x| x == 1).count();
    let mut cur = vec![0i64; nvar];
    let mut o = Vec::with_capacity(crate::ckalloc(cnt as usize * 32));
    for _ in 0..cnt {
        let mut k = 0usize;
        for i in 0..nf as usize {
            if i > 0 { o.push(sep); }
            if var[i] == 1 {
                let (u, n2) = vget(&v, q); q = n2;
                cur[k] = cur[k].wrapping_add(unzz(u) as i64);
                if cur[k] < 0 { crate::bail("corrupt names: a negative field"); }
                kst::putu64(&mut o, cur[k] as u64);
                k += 1;
            } else if var[i] == 2 {
                const HEX: &[u8; 16] = b"0123456789abcdef";
                for _ in 0..hw[i] / 2 {
                    if q >= v.len() { crate::bail("corrupt names: a hex field runs past the blob"); }
                    o.push(HEX[(v[q] >> 4) as usize]); o.push(HEX[(v[q] & 15) as usize]); q += 1;
                }
            } else { o.extend_from_slice(cf[i]); }
        }
        o.push(b'\n');
    }
    if q != v.len() { crate::bail("corrupt names: the field stream has more than the count"); }
    o
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
    // raw and coded length per chunk, as the C's tokd
    let (mut woff, mut coff) = (vec![0usize], vec![0usize]);
    for i in 0..nchunk {
        let at = |o: usize| u32::from_le_bytes(b[o..o + 4].try_into().unwrap()) as usize;
        woff.push(woff[i] + at(8 + 4 * i));
        coff.push(coff[i] + at(8 + 4 * (nchunk + i)));
    }
    if coff[nchunk] > b.len() - hdr {
        crate::bail("corrupt tokenised names: chunks overrun the stream");
    }
    // One chunk per thread into one preallocated buffer, so the peak is the
    // names once, not twice. The tokeniser separates with NUL; the format
    // stores newlines.
    use rayon::prelude::*;
    let mut out = vec![0u8; crate::ckalloc(woff[nchunk])];
    let mut slots: Vec<&mut [u8]> = Vec::with_capacity(nchunk);
    let mut rest = out.as_mut_slice();
    for i in 0..nchunk { let (a, r) = rest.split_at_mut(woff[i + 1] - woff[i]); slots.push(a); rest = r; }
    // At most TOKDEC chunks at once, as the C: each costs 127 MB in the decoder.
    let pool = rayon::ThreadPoolBuilder::new().num_threads(rayon::current_num_threads().min(TOKDEC)).build()
        .unwrap_or_else(|_| crate::bail("cannot start the name decode threads"));
    let ok = pool.install(|| slots.into_par_iter().enumerate().all(|(i, slot)| {
        match ffi::tok_decode(&b[hdr + coff[i]..hdr + coff[i + 1]]) {
            Some(d) if d.len() == slot.len() => {
                for (o, &c) in slot.iter_mut().zip(&d) { *o = if c == 0 { b'\n' } else { c }; }
                true
            }
            _ => false,
        }
    }));
    if !ok { crate::bail("name tokeniser decode failed"); }
    out
}

pub use crate::kst::{unzz, vget};

/// Tokenised name chunks decoded at once; 2026-10-02, 3000000 names, 8 threads:
/// 8 at once 1.54 s 1136 MB, 4 1.60 s 570 MB, 2 2.28 s 352 MB.
const TOKDEC: usize = 4;

/* ------------------------------------------------ one sequence at a time */

/// A position in the row stream: the block, its bytes, a cursor per frame,
/// and how many of its rows are left.
struct Cursor { blk: usize, data: std::sync::Arc<Vec<u8>>, lp: usize, ip: usize, left: u64 }

impl Reader {
    /// (block, offset in it, run count) for sequence `ri`. A forward cursor,
    /// so walking every sequence is linear and not quadratic.
    /// One index chunk decoded and checked, the cursor at its first sequence.
    fn load_ichunk(&mut self, c: u64) {
        if c >= self.coff.len() as u64 - 1 { crate::bail("index chunk out of range"); }
        let d = zdec(&self.iz[self.coff[c as usize] as usize..self.coff[c as usize + 1] as usize], "index chunk");
        let want = if c + 1 < self.cpend.len() as u64 - 1 { ICHUNK } else { self.nreads - c * ICHUNK };
        let (mut i, mut sum, mut k) = (0usize, 0u64, 0u64);
        while k < want && i < d.len() {
            let (v, n) = vget(&d, i);
            if v > u32::MAX as u64 { crate::bail("corrupt index: a run count past 32 bits"); }
            sum += v; i = n; k += 1;
        }
        if k != want || i != d.len() || self.cpend[c as usize] + sum != self.cpend[c as usize + 1] {
            crate::bail("corrupt index: an index chunk does not match its table entry");
        }
        self.ibx = d; self.cur_chunk = c;
        self.cur_ri = c * ICHUNK; self.cur_q = 0; self.cur_pend = self.cpend[c as usize];
    }

    fn locate(&mut self, ri: u64) -> (usize, u64, u64) {
        if self.runhead {
            if ri >= self.nreads { crate::bail("corrupt index: a sequence past the end"); }
            if self.rtab.is_empty() { self.load(0); }
            let start = self.rtab[ri as usize] as u64;
            let n = self.rtab[ri as usize + 1] as u64 - start;
            if start + n > self.rows { crate::bail("corrupt index: a sequence points outside the block table"); }
            return (0, start, n);
        }
        let want = ri / ICHUNK;
        if want != self.cur_chunk || ri < self.cur_ri { self.load_ichunk(want); }
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
    // inline, its rare branches cold: a call a row cost the Rust decode
    // more than the C's, whose loop has them in line
    #[inline]
    fn step(&mut self, c: &mut Cursor, st: u64, last: bool) -> (u64, u64) {
        if c.left == 0 { self.next_block(c); }
        let (v, l2) = vget(&c.data, c.lp);
        c.lp = l2;
        let len = if last && self.taildec {
            let t = self.tailpred.wrapping_add(unzz(v));
            if t <= st { crate::bail("corrupt track: sequence total is not past its last run"); }
            t - st
        } else if v != 0 {
            v
        } else {
            self.zero_len(c, st)
        };
        let (f, i2) = vget(&c.data, c.ip);
        c.ip = i2;
        c.left -= 1;
        (len, f)
    }

    #[cold]
    fn next_block(&mut self, c: &mut Cursor) {
        let b = if c.blk == usize::MAX { 0 } else { c.blk + 1 };
        c.data = std::sync::Arc::new(Vec::new());   // let go of the old block, so load can reuse it
        let (data, split, rows) = self.load(b);
        *c = Cursor { blk: b, data, lp: 0, ip: split, left: rows };
    }

    /// A zero length: --raw keeps it, a sentinel block asks the reference.
    #[cold]
    fn zero_len(&self, c: &Cursor, st: u64) -> u64 {
        if self.raw { return 0; }
        let lendec = self.lenfile && self.blk_nos.get(c.blk).copied().unwrap_or(0) != 1;
        if !lendec { crate::bail("corrupt track: a zero length outside a sentinel block"); }
        let q = self.alnq.unwrap_or_else(|| crate::bail("a length sentinel has no reference answer: wrong PAF or --ref"));
        let (e, _) = self.aln.as_ref().unwrap().end(q, st);
        if e <= st { crate::bail("a length sentinel has no reference answer: wrong PAF or --ref"); }
        e - st
    }

    /// The label bytes, or None when the table is empty and the id is the label.
    #[inline]
    fn label_at2(&self, f: u64, st: u64, len: u64) -> Option<&[u8]> {
        let l: &[u8] = match self.lab.get(f as usize) {
            Some(l) => l,
            None if self.lab.is_empty() => return None,
            None => crate::bail("corrupt block: a label id past the table"),
        };
        if self.aln.is_some() && l == b"=" { return Some(self.ref_label(st, len)); }
        Some(l)
    }

    #[cold]
    fn ref_label(&self, st: u64, len: u64) -> &[u8] {
        let a = self.aln.as_ref().unwrap();
        let q = self.alnq.unwrap_or_else(|| crate::bail("a delta row has no reference answer: wrong PAF or reference track"));
        let rl = a.row(q, st, len);
        if rl == u32::MAX { crate::bail("a delta row has no reference answer: wrong PAF or reference track"); }
        &self.lab[rl as usize]
    }

    fn label(&self, f: u64, out: &mut Vec<u8>, st: u64, len: u64) {
        match self.label_at2(f, st, len) {
            Some(l) => out.extend_from_slice(l),
            None => kst::putu64(out, f),
        }
    }

    /// What the C's ensurenames does: the names, by whichever codec stored
    /// them, or from --names FILE, or the ordinals.
    fn decode_names(&self) -> (Vec<u8>, Vec<u64>) {
        let (side, perm_ext, dropped) = (self.nf & 1 != 0, self.nf & 2 != 0, self.nf & 4 != 0);
        let (names_int, names_tok, was_perm, perm_tok) = (self.nf & 8 != 0, self.nf & 16 != 0, self.nf & 32 != 0, self.nf & 64 != 0);
        let names_fld = self.nf & 128 != 0;
        let (nz, pz, n) = (&self.nz, &self.pz, self.nreads);
        let mut names = if let (Some(f), false) = (&self.names_file, side || perm_ext) {
            if !dropped && !nz.is_empty() {
                crate::bail("this track stores its own names: --names FILE would replace them, so it is refused. Unset it to read the names in the file");
            }
            read_name_file(f, n)
        } else if dropped && nz.is_empty() {
            if !self.raw {
                crate::bail("the names of this track were dropped at compress time: give them back with --names FILE");
            }
            ordinals(n)
        } else if nz.is_empty() {
            ordinals(n)
        } else if names_fld {
            fld_names_d(nz, n)
        } else if names_int {
            int_names_d(nz, n)
        } else if names_tok {
            tok_d(nz)
        } else {
            zdec(nz, "names")
        };
        if names.last() == Some(&b'\n') { names.pop(); }
        let mut noff: Vec<u64> = vec![0];
        if !names.is_empty() {
            let mut at = 0usize;
            while let Some(k) = crate::ffi::find(b'\n', &names[at..]) { at += k + 1; noff.push(at as u64); }
            noff.push(names.len() as u64 + 1);
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
                    None if self.raw => None,
                    None => crate::bail("--sort --names drop keeps the rank of every record and no name text: give the original FASTA/FASTQ in --names FILE"),
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
                let (mut nn, mut no) = (Vec::new(), vec![0u64]);
                for i in 0..noff.len() - 1 {
                    let nm = &names[noff[i] as usize..noff[i + 1] as usize - 1];
                    let v: u64 = std::str::from_utf8(nm).ok().and_then(|t| t.parse().ok())
                        .unwrap_or_else(|| crate::bail("a stored rank is not a record number"));
                    if v >= n { crate::bail("a stored rank is not a record number"); }
                    nn.extend_from_slice(tab[v as usize]); nn.push(b'\n'); no.push(nn.len() as u64);
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

    /// The appended index, as the C reads it: the name's bucket, its
    /// fingerprints, and a name decoded only when one fingerprint is not enough.
    fn idxlookup(&self, nm: &[u8]) -> i64 {
        use std::os::unix::fs::FileExt;
        let Some((base, bind, k, pay, xf)) = &self.idx else { return -1 };
        let (base, bind, k, pay) = (*base, *bind, *k, *pay);
        let n = self.nreads;
        let (nb, w) = (1u64 << k, ixw(n.max(1)));
        let ns = nb.div_ceil(1024);
        if bind != self.mhash ^ self.namesrc() || pay != nb + ns * 8 + (n * (32 + w)).div_ceil(8) + 8 { return -1; }
        // only the name's bucket is read, as the C's mmap touches only it:
        // reading the whole index was 481 MB a query on 64M short reads
        let h = h64b(nm);
        let b = if k > 0 { h >> (64 - k) } else { 0 };
        let get = |off: u64, len: u64| { let mut v = vec![0u8; len as usize]; rd(xf.read_exact_at(&mut v, base + off), "index read"); v };
        let b0 = b & !1023;
        let cnt = get(b0, b - b0 + 1);
        let mut at = u64::from_le_bytes(get(nb + (b >> 10) * 8, 8)[..8].try_into().unwrap());
        for &c in &cnt[..cnt.len() - 1] { at += c as u64; }
        let (e0, e1) = (at, (at + cnt[cnt.len() - 1] as u64).min(n));
        let byte0 = (e0 * (32 + w)) >> 3;
        let ent = if e1 > e0 { get(nb + ns * 8 + byte0, ((e1 - e0) * (32 + w) + 7) / 8 + 8) } else { Vec::new() };
        let hits: Vec<u64> = (e0..e1).filter_map(|e| {
            let bit = e * (32 + w) - byte0 * 8;
            let o = ixget(&ent, bit + 32, w);
            (o < n && ixget(&ent, bit, 32) == ixfp(h, o)).then_some(o)
        }).collect();
        if hits.len() == 1 && self.names.get().is_none() { return hits[0] as i64; }   // one fingerprint: no name to decode
        hits.into_iter().find(|&o| self.name_at(o) == nm).map_or(-1, |o| o as i64)
    }

    /// As find, but a one-block track decodes its block before the names, as
    /// the C: the codec's working memory is gone when the names arrive.
    pub fn find_seq(&mut self, seq: &[u8]) -> i64 {
        let hit = self.idxlookup(seq);
        if hit >= 0 { return hit; }
        if self.b.len() == 1 { self.load(0); }
        self.find(seq)
    }

    pub fn find(&self, seq: &[u8]) -> i64 {
        let hit = self.idxlookup(seq);
        if hit >= 0 { return hit; }
        let (tab, mask) = self.nameix.get_or_init(|| {
            if self.nreads > u32::MAX as u64 { crate::bail("more than 2^32 sequences: the name table stores the ordinal in 32 bits"); }
            let mut cap = 1u64; while cap < self.nreads * 2 { cap <<= 1; }
            let (mut tab, mask) = (vec![u32::MAX; cap as usize], cap - 1);
            for i in 0..self.nreads {
                let nm = self.name_at(i);
                let mut s = (kst::fnv(FNV0, nm) & mask) as usize;
                while tab[s] != u32::MAX {
                    if self.name_at(tab[s] as u64) == nm { crate::bail("a sequence name repeats: the names are not usable"); }
                    s = (s + 1) & mask as usize;
                }
                tab[s] = i as u32;
            }
            (tab, mask)
        });
        let mut s = (kst::fnv(FNV0, seq) & mask) as usize;
        while tab[s] != u32::MAX {
            if self.name_at(tab[s] as u64) == seq { return tab[s] as i64; }
            s = (s + 1) & *mask as usize;
        }
        -1
    }

    /// region F SEQ | SEQ POS | SEQ START END [HALO].
    /// --db: the sequence's rows over the picked sets, touching equal tuples
    /// merged, so one set gives what its own track gives. Needs the whole sequence.
    #[allow(clippy::too_many_arguments)]
    fn region_set(&mut self, ri: u64, seq: &[u8], lo: u64, hi: u64,
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
            a = f; b = l + 1;
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

    #[allow(clippy::too_many_arguments)]
    pub fn region(&mut self, seq: &[u8], lo: u64, hi: u64,
                  feats: Option<&crate::db::Feats>, cols: &[usize], out: &mut impl std::io::Write) {
        let ri = self.find_seq(seq);
        if ri < 0 { crate::bail("sequence not found"); }
        let ri = ri as u64;
        self.alnq = self.aln.as_ref().and_then(|a| a.find(seq));
        if let Some(fe) = feats { return self.region_set(ri, seq, lo, hi, fe, cols, out); }
        let (blk, off, nruns) = self.locate(ri);
        let tskip = self.tailskip;
        let last_of = |k: u64| k + 1 == nruns && ri >= tskip;
        // As the C: a block that opens inside this sequence records where,
        // so with a START the walk begins at the last block opening at or
        // before it, not at the sequence's first run.
        let (start, bs) = (blk as u64 * self.block + off, self.block);
        let (mut jump, mut kskip, mut st0) = (blk, 0u64, 0u64);
        if lo > 0 {
            while jump + 1 < self.b.len() && (jump as u64 + 1) * bs < start + nruns
                  && self.b[jump + 1].bp != 0 && self.b[jump + 1].bp <= lo { jump += 1; }
        }
        if jump != blk { kskip = jump as u64 * bs - start; st0 = self.b[jump].bp; }

        let mut c = if jump != blk { self.at(jump, 0) } else { self.at(blk, off) };
        let mut buf = Vec::with_capacity(1 << 20);
        let mut st = st0;
        for k in kskip..nruns {
            let (len, f) = self.step(&mut c, st, last_of(k));
            let en = st + len;
            if en > lo && (hi == 0 || st < hi) {
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

pub(crate) fn skip(b: &[u8], mut i: usize) -> usize {
    while i < b.len() && b[i] >= 0x80 { i += 1; }
    i + 1
}

/* ------------------------------------------------------------ the index */

/// the name hash the index buckets and fingerprints by
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
        if self.nreads > u32::MAX as u64 {
            crate::bail("more than 2^32 sequences: the index stores the ordinal in 32 bits");
        }
        let n = self.nreads;
        let w = ixw(n.max(1));
        let mut k = 0u32;
        while k < 30 && (1u64 << (k + 2)) < n { k += 1; }   // about 4 reads a bucket
        let hs: Vec<u64> = (0..n).map(|i| h64b(self.name_at(i))).collect();
        let key = loop {                                   // a count must fit its byte
            let mut key: Vec<(u64, u64)> = hs.iter().enumerate()
                .map(|(i, &h)| ((if k > 0 { (h >> (64 - k)) << 32 } else { 0 }) | ixfp(h, i as u64), i as u64)).collect();
            key.sort_unstable();
            let top = key.chunk_by(|a, b| a.0 >> 32 == b.0 >> 32).map(|c| c.len()).max().unwrap_or(0);
            if top <= 255 { break key; }
            k += 1;
        };
        let (nb, dup) = (1u64 << k, key.windows(2).filter(|p| p[0].0 == p[1].0).count() as u64);
        let ns = nb.div_ceil(1024);
        let pay = nb + ns * 8 + (n * (32 + w)).div_ceil(8) + 8;
        let mut t = vec![0u8; crate::ckalloc(pay as usize)];
        let eo = (nb + ns * 8) as usize;
        for (i, &(kk, o)) in key.iter().enumerate() {
            t[(kk >> 32) as usize] += 1;
            let (mut bit, mut v) = (i as u64 * (32 + w), (kk & 0xffff_ffff) | (o << 32));
            let mut left = 32 + w;
            while left != 0 {
                let sh = bit & 7;
                let take = left.min(8 - sh);
                t[eo + (bit >> 3) as usize] |= ((v & ((1u64 << take) - 1)) << sh) as u8;
                v >>= take; bit += take; left -= take;
            }
        }
        let mut sum = 0u64;
        for g in 0..ns {
            t[(nb + g * 8) as usize..][..8].copy_from_slice(&sum.to_le_bytes());
            for b in g * 1024..nb.min((g + 1) * 1024) { sum += t[b as usize] as u64; }
        }
        for v in [n, k as u64, self.mhash ^ self.namesrc(), pay, IDXMAGIC] { t.extend_from_slice(&v.to_le_bytes()); }
        // beside the track in TRACK.kzi, written whole, then renamed into place
        use std::io::Write;
        let (ip, tp) = (format!("{path}.kzi"), format!("{path}.kzi.{}", std::process::id()));
        let ok = std::fs::File::create(&tp).and_then(|mut f| { f.write_all(&t)?; f.flush() }).is_ok()
            && std::fs::rename(&tp, &ip).is_ok();
        if !ok { let _ = std::fs::remove_file(&tp); crate::bail("index write failed"); }
        eprintln!("index: {n} sequences, {} B in {ip}, {dup} fingerprint collisions", pay + 40);
    }
}

/// the name's 32 bits, bound to its ordinal: a damaged entry fails to match
fn ixfp(h: u64, o: u64) -> u64 { (h ^ o.wrapping_mul(H64B_K)) & 0xffff_ffff }
/// bits for an ordinal below n, at least 1
fn ixw(n: u64) -> u64 { let mut w = 1; while w < 32 && (n - 1) >> w != 0 { w += 1; } w }
fn ixget(p: &[u8], bit: u64, nb: u64) -> u64 {
    let v = u64::from_le_bytes(p[(bit >> 3) as usize..][..8].try_into().unwrap());
    (v >> (bit & 7)) & ((1u64 << nb) - 1)
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
    pub fn names_shared(&self) -> std::sync::Arc<(Vec<u8>, Vec<u64>)> {
        self.names.get_or_init(|| std::sync::Arc::new(self.decode_names())).clone()
    }
    pub fn adopt_names(&mut self, n: std::sync::Arc<(Vec<u8>, Vec<u64>)>) { let _ = self.names.set(n); }
    pub fn labels(&self) -> &[Vec<u8>] { &self.lab }
    /// run count of sequence `ri`, from the index alone
    pub fn nruns_of(&mut self, ri: u64) -> u64 { self.locate(ri).2 }

    pub fn taildec(&self) -> bool { self.taildec }
    pub fn nblk(&self) -> usize { self.b.len() }
    pub fn blocksz(&self) -> u64 { self.block }
    /// where block `b` opens inside the sequence its first row belongs to
    pub fn blk_bp(&self, b: usize) -> u64 { self.b[b].bp }
    /// (first block, first row, run count) of sequence `ri`
    pub fn seq_span(&mut self, ri: u64) -> (usize, u64, u64) {
        let (blk, off, n) = self.locate(ri);
        (blk, blk as u64 * self.block + off, n)
    }
    /// The rows of sequence `ri` that lie in block `b`, as (start, label id),
    /// and the position after the last of them. The sequence's own
    /// coordinates: block `b` opens at blk_bp(b) unless it is the first.
    pub fn blk_runs(&mut self, ri: u64, b: usize) -> (Vec<(u64, u32)>, u64) {
        let (fb, row0, n) = self.seq_span(ri);
        let (mut row, mut pos) = if b == fb { (row0, 0) } else { (b as u64 * self.block, self.b[b].bp) };
        let mut c = self.at(b, row - b as u64 * self.block);
        let end = (b as u64 * self.block + self.b[b].rows).min(row0 + n);
        let mut out = Vec::with_capacity((end - row) as usize);
        while row < end {
            let (len, f) = self.step(&mut c, pos, false);
            if len == 0 { crate::bail("corrupt track: a zero length outside a sentinel block"); }
            if f as usize >= self.lab.len() { crate::bail("corrupt block: a label id past the table"); }
            out.push((pos, f as u32));
            pos += len;
            row += 1;
        }
        (out, pos)
    }

    /// One sequence of a plain track as run starts, label ids and its end:
    /// what a --ref contig is, decoded when first asked.
    pub fn seq_runs(&mut self, ri: u64) -> (crate::aln::Starts, Vec<u32>, u64) {
        use crate::aln::Starts;
        let (blk, off, nruns) = self.locate(ri);
        let mut c = self.at(blk, off);
        // starts into 4 bytes from the first row; a start past 32 bits, which
        // no chromosome reaches, moves what is there into 8 and goes on
        let (mut st, mut lab) = (Starts::W(Vec::with_capacity(nruns as usize)), Vec::with_capacity(nruns as usize));
        let mut s = 0u64;
        for k in 0..nruns {
            let (len, f) = self.step(&mut c, s, k + 1 == nruns && ri >= self.tailskip);
            if f as usize >= self.lab.len() { crate::bail("corrupt block: a label id past the table"); }
            if len == 0 { crate::bail("reference track: a run of no length"); }
            if let Starts::W(w) = &st { if s >> 32 != 0 { st = Starts::L(w.iter().map(|&x| x as u64).collect()); } }
            match &mut st { Starts::W(w) => w.push(s as u32), Starts::L(l) => l.push(s) }
            lab.push(f as u32); s += len;
        }
        self.sc = None;
        (st, lab, s)
    }

    /// One sequence as arrays: starts, ends, hierarchy node per run, and
    /// whether the label was "novel". What smoothing works on.
    /// The four arrays are the caller's and are cleared here, so a walk over
    /// many sequences allocates once, as the C's r->sst and friends do.
    pub fn seq_rows(&mut self, ri: u64, labmap: &[i32], isnovel: &[u8], root: i32,
                    out: &mut SeqRows) {
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
        let SeqRows { st, en, feat, nov } = out;
        for v in [&mut *st, &mut *en] { v.clear(); v.reserve(n); }
        feat.clear(); feat.reserve(n);
        nov.clear(); nov.reserve(n);
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
    }
}

/// One sequence as arrays, reused across sequences.
#[derive(Default)]
pub struct SeqRows { pub st: Vec<u64>, pub en: Vec<u64>, pub feat: Vec<i32>, pub nov: Vec<u8> }
