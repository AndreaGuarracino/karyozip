//! The KaryoScope Track (KST) format: writer and reader, the same bytes as the C in `../c`

pub mod aln;
pub mod db;
pub mod ffi;
pub mod kst;
pub mod reader;
pub mod smooth;
pub mod writer;

thread_local! {
    /// set on a decode-ahead thread: a refusal is kept for its block
    pub static DEFER: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}
/// A decode-ahead thread's refusal, carried to the walk
pub struct Deferred(pub String);

/// Every refusal is the C's own wording, and exits 1 as the C does
pub fn bail(msg: &str) -> ! {
    if DEFER.with(|d| d.get()) { std::panic::panic_any(Deferred(msg.to_string())); }
    eprintln!("karyozip: {msg}");
    db::drop_tmp();
    ffi::drop_out();
    prof::dump();
    std::process::exit(1);
}

/// The C's ALLOCMAX: a forged length is a refusal, not an allocator panic
pub fn ckalloc(n: usize) -> usize {
    if n > 1 << 38 { bail("corrupt file: an impossible allocation size"); }
    n
}

/// A BED coordinate: digits only, and it must fit in 64 bits
pub fn bedcoord_r(b: &[u8], what: &'static str) -> Result<u64, &'static str> {
    if b.is_empty() { return Err("empty start or end column"); }
    let mut v: u64 = 0;
    for &c in b {
        if !c.is_ascii_digit() { return Err(what); }
        if v > (u64::MAX - (c - b'0') as u64) / 10 {
            return Err("a coordinate does not fit in 64 bits");
        }
        v = v * 10 + (c - b'0') as u64;
    }
    Ok(v)
}

pub fn bedcoord(b: &[u8], what: &'static str) -> u64 {
    bedcoord_r(b, what).unwrap_or_else(|e| bail(e))
}

/// A whole text input, plain or gzipped, through the same reader as a stream
pub fn slurp(path: &str, what: &'static str) -> Vec<u8> {
    let f = std::fs::File::open(path).unwrap_or_else(|_| bail(what));
    let mut v = Vec::new();
    std::io::Read::read_to_end(&mut unzip(f), &mut v).unwrap_or_else(|_| bail("read error on a companion file"));
    v
}

/// Every text input, as the C's zfopen: plain, BGZF blocks inflated side by side and CRC-checked, or any gzip
pub fn unzip(r: impl std::io::Read + Send + 'static) -> Box<dyn std::io::Read + Send> {
    Box::new(Zr { r: Box::new(r), mode: 0, reof: false, raw: Vec::new(), rp: 0, out: Vec::new(), op: 0, gz: None })
}

const ZCHUNK: usize = 8 << 20;
type Gz = flate2::read::MultiGzDecoder<std::io::Chain<std::io::Cursor<Vec<u8>>, Box<dyn std::io::Read + Send>>>;
struct Zr { r: Box<dyn std::io::Read + Send>, mode: u8, reof: bool, raw: Vec<u8>, rp: usize, out: Vec<u8>, op: usize, gz: Option<Gz> }

/// One BGZF block's size, 0 if it is not one
fn bgzf_size(p: &[u8]) -> usize {
    if p.len() < 18 || p[0] != 0x1f || p[1] != 0x8b || p[2] != 8 || p[3] & 4 == 0 { return 0; }
    let xl = p[10] as usize | (p[11] as usize) << 8;
    let mut i = 12;
    while i + 4 <= 12 + xl && i + 4 <= p.len() {
        let sl = p[i + 2] as usize | (p[i + 3] as usize) << 8;
        if p[i] == b'B' && p[i + 1] == b'C' && sl == 2 { return (p[i + 4] as usize | (p[i + 5] as usize) << 8) + 1; }
        i += 4 + sl;
    }
    0
}
fn le32(p: &[u8]) -> usize { u32::from_le_bytes([p[0], p[1], p[2], p[3]]) as usize }
fn zbad() -> std::io::Error { std::io::Error::other("bad or cut gzip input") }

impl Zr {
    /// At least `want` unused bytes in raw, or the end
    fn fill(&mut self, want: usize) -> std::io::Result<()> {
        self.raw.drain(..self.rp); self.rp = 0;
        let mut b = vec![0u8; 1 << 20];
        while self.raw.len() < want && !self.reof {
            match self.r.read(&mut b) {
                Ok(0) => self.reof = true,
                Ok(k) => self.raw.extend_from_slice(&b[..k]),
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
                Err(e) => return Err(e),
            }
        }
        Ok(())
    }
    /// The next run of whole blocks into out: 1, 0 at the end, or 2 when what follows is not a whole BGZF block
    fn bgzf(&mut self) -> std::io::Result<u8> {
        use rayon::prelude::*;
        self.fill(ZCHUNK)?;
        let (mut bo, mut q) = (Vec::new(), self.rp);
        while q < self.raw.len() {
            let bs = bgzf_size(&self.raw[q..]);
            if bs < 26 || q + bs > self.raw.len() { break; }
            bo.push(q); q += bs;
        }
        bo.push(q);
        if bo.len() == 1 { return Ok(if self.rp < self.raw.len() { 2 } else { 0 }); }
        let mut oo = vec![0usize];
        for w in bo.windows(2) { oo.push(oo.last().unwrap() + le32(&self.raw[w[1] - 4..])); }
        self.out = vec![0u8; *oo.last().unwrap()]; self.op = 0;
        let mut slots: Vec<&mut [u8]> = Vec::with_capacity(bo.len() - 1);
        let mut rest = self.out.as_mut_slice();
        for w in oo.windows(2) { let (a, r) = rest.split_at_mut(w[1] - w[0]); slots.push(a); rest = r; }
        let raw = &self.raw;
        // libdeflate, one decompressor a thread, as the C's
        let ok = slots.into_par_iter().enumerate().with_min_len(16).map_init(ffi::Inflater::new, |d, (i, o)| {
            let b = &raw[bo[i]..bo[i + 1]];
            let xl = b[10] as usize | (b[11] as usize) << 8;
            if o.len() > 65536 || 12 + xl + 8 > b.len() { return false; }
            d.inflate(&b[12 + xl..b.len() - 8], o) && ffi::crc32(o) as usize == le32(&b[b.len() - 8..])
        }).all(|x| x);
        self.rp = q;
        if ok { Ok(1) } else { Err(zbad()) }
    }
    fn to_stream(&mut self) {
        let rest = std::io::Cursor::new(self.raw.split_off(self.rp));
        let r = std::mem::replace(&mut self.r, Box::new(std::io::empty()));
        self.gz = Some(flate2::read::MultiGzDecoder::new(std::io::Read::chain(rest, r)));
        self.mode = 3;
    }
}

impl std::io::Read for Zr {
    fn read(&mut self, b: &mut [u8]) -> std::io::Result<usize> {
        if self.mode == 0 {
            self.fill(18)?;
            let p = &self.raw[self.rp..];
            self.mode = if p.len() < 2 || p[0] != 0x1f || p[1] != 0x8b { 1 } else if bgzf_size(p) > 0 { 2 } else { 3 };
            if self.mode == 3 { self.to_stream(); }
        }
        if self.mode == 1 {                        // plain: what was peeked, then the reader
            if self.rp < self.raw.len() {
                let k = (self.raw.len() - self.rp).min(b.len());
                b[..k].copy_from_slice(&self.raw[self.rp..self.rp + k]); self.rp += k;
                return Ok(k);
            }
            return self.r.read(b);
        }
        while self.mode == 2 && self.op == self.out.len() {
            match self.bgzf()? {
                0 => return Ok(0),
                2 => { if self.reof && self.raw.len() - self.rp < 18 { return Err(zbad()); } self.to_stream(); }
                _ => {}
            }
        }
        if self.mode == 2 {
            let k = (self.out.len() - self.op).min(b.len());
            b[..k].copy_from_slice(&self.out[self.op..self.op + k]); self.op += k;
            return Ok(k);
        }
        self.gz.as_mut().unwrap().read(b)
    }
}

/// A text input one line at a time, plain or gzipped, without holding the file
pub struct Lines { r: std::io::BufReader<Box<dyn std::io::Read + Send>>, buf: Vec<u8>, pend: usize, pub nl: bool, pub err: &'static str }
impl Lines {
    pub fn open(path: &str, what: &'static str) -> Lines {
        let f = std::fs::File::open(path).unwrap_or_else(|_| bail(what));
        Lines::new(f)
    }
    pub fn new(r: impl std::io::Read + Send + 'static) -> Lines {
        Lines { r: std::io::BufReader::with_capacity(1 << 20, unzip(r)), buf: Vec::new(), pend: 0, nl: true, err: "read error" }
    }
    /// The next line without its newline, or None at the end; `nl` says whether the line had one
    pub fn next_line(&mut self) -> Option<&[u8]> { self.next_line_nl().map(|(l, _)| l) }
    /// The line and whether a newline ended it, in place when it lies inside the buffer
    pub fn next_line_nl(&mut self) -> Option<(&[u8], bool)> {
        use std::io::BufRead;
        self.r.consume(std::mem::take(&mut self.pend));
        let hit = crate::ffi::find(b'\n', self.r.fill_buf().unwrap_or_else(|_| bail(self.err)));
        if let Some(i) = hit {
            self.pend = i + 1; self.nl = true;
            return Some((&self.r.buffer()[..i], true));
        }
        self.buf.clear();
        match self.r.read_until(b'\n', &mut self.buf) {
            Ok(0) => None,
            Ok(_) => { self.nl = self.buf.last() == Some(&b'\n'); if self.nl { self.buf.pop(); } Some((&self.buf, self.nl)) }
            Err(_) => bail(self.err),
        }
    }
}

pub mod prof {
    // --profile: the C's phases and closing line; times are the machine's, never compared
    use std::sync::atomic::{AtomicBool, AtomicU64, Ordering::Relaxed};
    use std::sync::OnceLock;
    use std::time::Instant;

    pub static ON: AtomicBool = AtomicBool::new(false);
    static NS: [AtomicU64; 10] = [const { AtomicU64::new(0) }; 10];
    static HW: [AtomicU64; 10] = [const { AtomicU64::new(0) }; 10];
    static THR: AtomicU64 = AtomicU64::new(0);
    static T0: OnceLock<Instant> = OnceLock::new();
    static NAME: [&str; 10] = ["parse", "compress", "names", "index", "smooth",
                              "parse:lines", "parse:append", "write", "codec search", "layout search"];

    pub fn start(on: bool, threads: usize) {
        ON.store(on, Relaxed);
        THR.store(threads as u64, Relaxed);
        let _ = T0.set(Instant::now());
    }

    /// Charges its phase for as long as it lives
    pub struct Span(usize, Instant);
    pub fn span(i: usize) -> Option<Span> {
        (ON.load(Relaxed) && i < NAME.len()).then(|| Span(i, Instant::now()))   // 10 is "not this phase"
    }
    /// The high-water mark is read from /proc at most every 10 ms, as the C.
    static LAST: AtomicU64 = AtomicU64::new(0);
    impl Drop for Span {
        fn drop(&mut self) {
            NS[self.0].fetch_add(self.1.elapsed().as_nanos() as u64, Relaxed);
            let now = T0.get().map_or(0, |t| t.elapsed().as_nanos() as u64);
            if now - LAST.load(Relaxed) > 10_000_000 {
                LAST.store(now, Relaxed);
                HW[self.0].store(field("/proc/self/status", "VmHWM:"), Relaxed);
            }
        }
    }

    /// One kB number from a "Key: N kB" line, or 0 when /proc is not there
    fn field(path: &str, key: &str) -> u64 {
        std::fs::read_to_string(path).ok()
            .and_then(|t| t.lines().find(|l| l.starts_with(key))
                           .and_then(|l| l.split_whitespace().nth(1)?.parse().ok()))
            .unwrap_or(0)
    }

    /// utime plus stime in seconds, from /proc/self/stat in 100 Hz ticks, split after the comm field's ')'
    fn cpu() -> f64 {
        let Ok(t) = std::fs::read_to_string("/proc/self/stat") else { return 0.0 };
        let Some(i) = t.rfind(')') else { return 0.0 };
        let f: Vec<&str> = t[i + 1..].split_whitespace().collect();
        let get = |k: usize| f.get(k).and_then(|v| v.parse::<f64>().ok()).unwrap_or(0.0);
        (get(11) + get(12)) / 100.0
    }

    pub fn dump() {
        if !ON.load(Relaxed) { return; }
        let tot: f64 = NS.iter().map(|n| n.load(Relaxed) as f64).sum::<f64>() / 1e9;
        for i in 0..NAME.len() {
            let s = NS[i].load(Relaxed) as f64 / 1e9;
            if s > 0.0 {
                eprintln!("  phase {:<12} {:8.3} s  {:5.1}%   hwm {:5} MB",
                          NAME[i], s, 100.0 * s / if tot > 0.0 { tot } else { 1.0 },
                          HW[i].load(Relaxed) / 1024);
            }
        }
        let wall = T0.get().map_or(0.0, |t| t.elapsed().as_secs_f64());
        let (c, hwm) = (cpu(), field("/proc/self/status", "VmHWM:"));
        eprintln!("  cpu {c:.3} s over {wall:.3} s wall = {:.2} cores busy \
                   ({} allowed), peak rss {:.2} GB",
                  if wall > 0.0 { c / wall } else { 0.0 },
                  THR.load(Relaxed), hwm as f64 / 1048576.0);
    }
}
