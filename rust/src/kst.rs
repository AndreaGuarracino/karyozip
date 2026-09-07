//! Format primitives shared by the writer and, later, the reader.
//! The format's number encoding: LEB128 unsigned, and zigzag for signed.

pub fn put(out: &mut Vec<u8>, mut v: u64) {
    while v >= 0x80 {
        out.push((v as u8 & 0x7f) | 0x80);
        v >>= 7;
    }
    out.push(v as u8);
}

/// LEB128 back out, with the bound the forged-file fuzzing needs.
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

pub fn unzz(v: u64) -> u64 { ((v >> 1) as i64 ^ -((v & 1) as i64)) as u64 }

pub fn zz(v: i64) -> u64 {
    ((v as u64) << 1) ^ ((v >> 63) as u64)
}

/// FNV-1a, the hash the format uses everywhere.
pub const FNV0: u64 = 0xcbf2_9ce4_8422_2325;
/// one row into the row hash: name hash, start, end, label hash, as words
pub fn rowmix(h: u64, nh: u64, st: u64, en: u64, lh: u64) -> u64 {
    let m = |h: u64, v: u64| (h ^ v).wrapping_mul(0x0100_0000_01b3);
    m(m(m(m(h, nh), st), en), lh)
}

pub fn fnv(mut h: u64, b: &[u8]) -> u64 {
    for &x in b {
        h ^= x as u64;
        h = h.wrapping_mul(0x100_0000_01b3);
    }
    h
}

/// Decimal digits straight into the buffer. `to_string()` allocates a String
/// per number, which is a quarter of decode's time on a 40M-row track.
pub fn putu64(out: &mut Vec<u8>, v: u64) {
    let mut b = [0u8; 20];
    let n = putu64_at(&mut b, v);
    out.extend_from_slice(&b[..n]);
}

/// The digits at the front of `b`, and how many: one row of a BED goes into
/// a stack buffer and reaches the output in one call, as the C writes it
/// through one pointer.
pub fn putu64_at(b: &mut [u8], v: u64) -> usize {
    // the digit count by compares, then one division a digit, in place:
    // a count by division doubled the divisions and a scratch array cost
    // a memcpy, and this is every coordinate of every row written
    const P: [u64; 19] = [10, 100, 1_000, 10_000, 100_000, 1_000_000, 10_000_000, 100_000_000,
        1_000_000_000, 10_000_000_000, 100_000_000_000, 1_000_000_000_000, 10_000_000_000_000,
        100_000_000_000_000, 1_000_000_000_000_000, 10_000_000_000_000_000,
        100_000_000_000_000_000, 1_000_000_000_000_000_000, 10_000_000_000_000_000_000];
    let n = 1 + P.iter().take_while(|&&p| v >= p).count();
    let (mut i, mut t) = (n, v);
    loop {
        i -= 1;
        b[i] = b'0' + (t % 10) as u8;
        t /= 10;
        if i == 0 { return n; }
    }
}

/// FNV-1a as a `HashMap` hasher. The default is SipHash, which on a
/// one-byte label costs more than the lookup it protects, and these keys
/// come from the file the caller already chose to read.
#[derive(Default, Clone, Copy)]
pub struct Fnv;
pub struct FnvHasher(u64);

impl std::hash::BuildHasher for Fnv {
    type Hasher = FnvHasher;
    fn build_hasher(&self) -> FnvHasher { FnvHasher(FNV0) }
}
impl std::hash::Hasher for FnvHasher {
    fn finish(&self) -> u64 { self.0 }
    fn write(&mut self, b: &[u8]) { self.0 = fnv(self.0, b); }
}

pub type FnvMap<K, V> = std::collections::HashMap<K, V, Fnv>;
