//! The KaryoScope Track (KST) format: a lossless, seekable container for
//! k-mer annotation tracks. `writer::Writer` writes one, `reader::Reader`
//! reads one, and both produce the same bytes as the C in `../c`.

pub mod aln;
pub mod db;
pub mod ffi;
pub mod kst;
pub mod prof;
pub mod reader;
pub mod smooth;
pub mod writer;

/// Every refusal is the C's own wording, and exits 1 as the C does.
pub fn bail(msg: &str) -> ! {
    eprintln!("karyozip: {msg}");
    db::drop_tmp();
    prof::dump();
    std::process::exit(1);
}

/// The C's ALLOCMAX. A length read out of a file can say anything, and the
/// answer to a forged one is a refusal, not an allocator panic.
pub fn ckalloc(n: usize) -> usize {
    if n > 1 << 38 { bail("corrupt file: an impossible allocation size"); }
    n
}

/// A BED coordinate: digits only, and it must fit in 64 bits.
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

/// A whole text input, inflated when it is gzipped. zlib and this agree on
/// the bytes, so the C and this read the same file the same way.
pub fn slurp(path: &str, what: &'static str) -> Vec<u8> {
    let d = std::fs::read(path).unwrap_or_else(|_| bail(what));
    if d.len() > 1 && d[0] == 0x1f && d[1] == 0x8b { inflate(&d[..], what) } else { d }
}

/// The same, for a stream: stdin may be a pipe, so the first bytes are
/// peeked and put back rather than seeked over.
pub fn unzip(r: impl std::io::Read + 'static) -> Box<dyn std::io::Read> {
    let mut b = std::io::BufReader::with_capacity(1 << 16, r);
    let gz = matches!(std::io::BufRead::fill_buf(&mut b), Ok(p) if p.len() > 1 && p[0] == 0x1f && p[1] == 0x8b);
    if gz { Box::new(flate2::read::MultiGzDecoder::new(b)) } else { Box::new(b) }
}

/// A text input one line at a time, plain or gzipped, without holding the
/// file: a reference track is a gigabyte of text for 24 million runs.
pub struct Lines { r: std::io::BufReader<Box<dyn std::io::Read>>, buf: Vec<u8>, pub nl: bool, pub err: &'static str }
impl Lines {
    pub fn open(path: &str, what: &'static str) -> Lines {
        let f = std::fs::File::open(path).unwrap_or_else(|_| bail(what));
        Lines::new(f)
    }
    pub fn new(r: impl std::io::Read + 'static) -> Lines {
        Lines { r: std::io::BufReader::with_capacity(1 << 20, unzip(r)), buf: Vec::new(), nl: true, err: "read error" }
    }
    /// The next line without its newline, or None at the end; `nl` says
    /// whether the line had one.
    pub fn next_line(&mut self) -> Option<&[u8]> {
        use std::io::BufRead;
        self.buf.clear();
        match self.r.read_until(b'\n', &mut self.buf) {
            Ok(0) => None,
            Ok(_) => { self.nl = self.buf.last() == Some(&b'\n'); if self.nl { self.buf.pop(); } Some(&self.buf) }
            Err(_) => bail(self.err),
        }
    }
}

fn inflate(d: &[u8], what: &'static str) -> Vec<u8> {
    use std::io::Read;
    let mut v = Vec::new();
    flate2::read::MultiGzDecoder::new(d).read_to_end(&mut v).unwrap_or_else(|_| bail(what));
    v
}
