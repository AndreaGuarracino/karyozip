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
    let mut i = b.len();
    let mut n = v;
    loop {
        i -= 1;
        b[i] = b'0' + (n % 10) as u8;
        n /= 10;
        if n == 0 { break; }
    }
    out.extend_from_slice(&b[i..]);
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
