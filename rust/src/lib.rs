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

/// A BED coordinate: digits only, and it must fit in 64 bits.
pub fn bedcoord_r(b: &[u8], what: &'static str) -> Result<u64, &'static str> {
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
