//! The format's number encoding: LEB128 unsigned, and zigzag for signed.

pub fn put(out: &mut Vec<u8>, mut v: u64) {
    while v >= 0x80 {
        out.push((v as u8 & 0x7f) | 0x80);
        v >>= 7;
    }
    out.push(v as u8);
}

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
