//! Format primitives shared by the writer and, later, the reader.
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
