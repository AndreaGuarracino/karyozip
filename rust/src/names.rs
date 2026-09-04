//! The name stream: which of three codecs stores it, and the refusal when a
//! sequence name repeats.

use crate::ffi;
use crate::varint::{self, FNV0};
use std::collections::HashSet;

/// 0 zstd, 1 integer deltas, 2 the CRAM tokeniser.
pub struct Best {
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
        varint::put(&mut v, varint::zz(d));
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
pub fn best(names: &[u8], n: u64, level: i32) -> Best {
    let (z, t) = rayon::join(
        || ffi::zstd(names, level),
        || tok_best(names),
    );
    let i = int_names(names, n).map(|v| ffi::zstd(&v, level));
    let mut out = Best { bytes: z, kind: 0 };
    if let Some(v) = i {
        if !v.is_empty() && v.len() < out.bytes.len() { out = Best { bytes: v, kind: 1 }; }
    }
    if let Some(v) = t {
        if !v.is_empty() && v.len() < out.bytes.len() { out = Best { bytes: v, kind: 2 }; }
    }
    out
}

/// A repeated sequence name would make the index ambiguous, so it is refused.
pub fn check_unique(names: &[u8]) -> Result<(), String> {
    let mut seen: HashSet<u64> = HashSet::new();
    let mut exact: HashSet<&[u8]> = HashSet::new();
    for line in names.split(|&c| c == b'\n') {
        if line.is_empty() { continue; }
        let h = varint::fnv(FNV0, line);
        if !seen.insert(h) && !exact.insert(line) {
            return Err(format!("sequence {} appears twice",
                               String::from_utf8_lossy(line)));
        }
        exact.insert(line);
    }
    Ok(())
}
