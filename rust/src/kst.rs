//! Format primitives shared by the writer and the reader: LEB128 unsigned, zigzag for signed

/// The header flags in h[1], above the 49-bit names length
pub const TOKFLAG: u64 = 1 << 63;   // names through the tokeniser
pub const INTFLAG: u64 = 1 << 62;   // names as integer deltas
pub const RUNFLAG: u64 = 1 << 61;   // run counts head each read's lengths
pub const DROPPED: u64 = 1 << 60;   // names not stored
pub const ALNFLAG: u64 = 1 << 59;   // rows deltaed against a reference
pub const DBFLAG: u64 = 1 << 58;    // labels are featureIDs
pub const LENFLAG: u64 = 1 << 56;   // lengths the reference implies are 0
pub const TAILFLAG: u64 = 1 << 55;  // last run lengths predicted
pub const CARTFLAG: u64 = 1 << 52;  // vocabulary featureIDs
pub const ILVFLAG: u64 = 1 << 51;   // labels interleaved with lengths
pub const FLDFLAG: u64 = 1 << 50;   // names as field deltas
pub const GRPFLAG: u64 = 1 << 49;   // lengths grouped by label
pub const VARFLAG: u64 = 1 << 57;   // blocks of whole reads, a read count each
pub const FLAGS: u64 = TOKFLAG | INTFLAG | RUNFLAG | DROPPED | ALNFLAG | DBFLAG | LENFLAG
                      | TAILFLAG | CARTFLAG | ILVFLAG | FLDFLAG | GRPFLAG | VARFLAG;
pub const NMLBITS: u32 = 49;        // the names length below the flags
pub const ICHUNK: u64 = 65536;      // sequences per index chunk

#[inline]
pub fn put(out: &mut Vec<u8>, mut v: u64) {
    while v >= 0x80 {
        out.push((v as u8 & 0x7f) | 0x80);
        v >>= 7;
    }
    out.push(v as u8);
}

/// LEB128 back out, with the bound the forged-file fuzzing needs
#[inline]
pub fn vget(b: &[u8], i: usize) -> (u64, usize) {
    match b.get(i) {   // one byte is most values: no loop for it
        Some(&c) if c < 0x80 => (c as u64, i + 1),
        _ => vget_long(b, i),
    }
}
#[cold]
fn vget_long(b: &[u8], mut i: usize) -> (u64, usize) {
    let (mut v, mut s) = (0u64, 0u32);
    loop {
        if i >= b.len() { crate::bail("corrupt varint"); }
        let c = b[i];
        i += 1;
        if s == 63 && c > 1 { crate::bail("corrupt varint"); }   // the tenth byte holds bit 63 alone and ends it
        v |= ((c & 0x7f) as u64) << s;
        s += 7;
        if c < 0x80 { return (v, i); }
    }
}

/// `vget` with the refusal returned, for the loaders of the write path
pub fn vget_r(b: &[u8], mut i: usize) -> crate::KzResult<(u64, usize)> {
    let (mut v, mut s) = (0u64, 0u32);
    let bad = || crate::Error::input("corrupt varint");
    loop {
        let c = *b.get(i).ok_or_else(bad)?;
        i += 1;
        if s == 63 && c > 1 { return Err(bad()); }   // the tenth byte holds bit 63 alone and ends it
        v |= ((c & 0x7f) as u64) << s;
        s += 7;
        if c < 0x80 { return Ok((v, i)); }
    }
}

/// The writer's own LEB128, which it wrote whole: decoded with no refusal (a cut one reads to the end), so
/// nothing on the write path exits; the same values as `vget` on what `put` writes
#[inline]
pub fn vget_w(b: &[u8], i: usize) -> (u64, usize) {
    match b.get(i) {
        Some(&c) if c < 0x80 => (c as u64, i + 1),
        _ => vget_soft(b, i),
    }
}
#[cold]
fn vget_soft(b: &[u8], mut i: usize) -> (u64, usize) {
    let (mut v, mut s) = (0u64, 0u32);
    while i < b.len() {
        let c = b[i];
        i += 1;
        if s < 64 { v |= ((c & 0x7f) as u64) << s; }
        s += 7;
        if c < 0x80 { break; }
    }
    (v, i)
}
/// `reader::skip` for the writer's own varints: past one, never past the end, no refusal
#[inline]
pub fn skip_w(b: &[u8], i: usize) -> usize {
    let mut j = i;
    while j < b.len() && b[j] >= 0x80 { j += 1; }
    (j + 1).min(b.len().max(i))
}

pub fn unzz(v: u64) -> u64 { ((v >> 1) as i64 ^ -((v & 1) as i64)) as u64 }

pub fn zz(v: i64) -> u64 {
    ((v as u64) << 1) ^ ((v >> 63) as u64)
}

/// FNV-1a, the hash the format uses everywhere
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

/// Decimal digits straight into the buffer, as the C's putu64: to_string() allocates a String per number
#[inline]
pub fn put_row(buf: &mut Vec<u8>, name: &[u8], st: u64, en: u64, lab: &[u8]) {
    buf.reserve(name.len() + 44 + lab.len() + 1);
    // SAFETY: reserved above, the name plus 44 bytes for three tabs and two 20-digit numbers, the label and the newline
    unsafe {
        let n0 = buf.len();
        let p = buf.as_mut_ptr().add(n0);
        std::ptr::copy_nonoverlapping(name.as_ptr(), p, name.len());
        let row = std::slice::from_raw_parts_mut(p.add(name.len()), 44);
        row[0] = b'\t';
        let mut nl = 1 + putu64_at(&mut row[1..], st);
        row[nl] = b'\t'; nl += 1;
        nl += putu64_at(&mut row[nl..], en);
        row[nl] = b'\t'; nl += 1;
        let q = p.add(name.len() + nl);
        std::ptr::copy_nonoverlapping(lab.as_ptr(), q, lab.len());
        *q.add(lab.len()) = b'\n';
        buf.set_len(n0 + name.len() + nl + lab.len() + 1);
    }
}

pub fn putu64(out: &mut Vec<u8>, v: u64) {
    let mut b = [0u8; 20];
    let n = putu64_at(&mut b, v);
    out.extend_from_slice(&b[..n]);
}

/// The digits at the front of `b`, and how many: a row goes out in one call, as the C writes it
pub fn putu64_at(b: &mut [u8], v: u64) -> usize {
    // the digit count by compares, then one division a digit, in place: this writes every coordinate
    const P: [u64; 19] = [10, 100, 1_000, 10_000, 100_000, 1_000_000, 10_000_000, 100_000_000,
        1_000_000_000, 10_000_000_000, 100_000_000_000, 1_000_000_000_000, 10_000_000_000_000,
        100_000_000_000_000, 1_000_000_000_000_000, 10_000_000_000_000_000,
        100_000_000_000_000_000, 1_000_000_000_000_000_000, 10_000_000_000_000_000_000];
    const D2: &[u8; 200] = b"00010203040506070809101112131415161718192021222324\
2526272829303132333435363738394041424344454647484950515253545556575859\
60616263646566676869707172737475767778798081828384858687888990919293949596979899";
    let mut n = (((63 - (v | 1).leading_zeros()) * 1233) >> 12) as usize + 1;
    if n <= 19 && v >= P[n - 1] { n += 1; }   // the estimate is off by one
    // n <= 20 and the table index is under 200, so the digit writes need no check
    assert!(b.len() >= n && n <= 20);
    let (mut i, mut t) = (n, v);
    unsafe {
        while t >= 100 {
            let q = (t % 100) as usize; t /= 100; i -= 2;
            *b.get_unchecked_mut(i) = *D2.get_unchecked(2 * q);
            *b.get_unchecked_mut(i + 1) = *D2.get_unchecked(2 * q + 1);
        }
        if t < 10 { i -= 1; *b.get_unchecked_mut(i) = b'0' + t as u8; }
        else { i -= 2;
            *b.get_unchecked_mut(i) = *D2.get_unchecked(2 * t as usize);
            *b.get_unchecked_mut(i + 1) = *D2.get_unchecked(2 * t as usize + 1); }
    }
    n
}

/// FNV-1a as a `HashMap` hasher: SipHash costs more than the lookup on one-byte labels
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
    fn write_u32(&mut self, v: u32) { self.0 = (self.0 ^ v as u64).wrapping_mul(0x100_0000_01b3); }
}

/// A tuple of label ids hashed a word at a time, as the C's tuphash
#[derive(Clone, Default, PartialEq, Eq)]
pub struct Tup(pub Vec<u32>);
impl std::hash::Hash for Tup {
    fn hash<H: std::hash::Hasher>(&self, h: &mut H) { for &v in &self.0 { h.write_u32(v); } }
}

pub type FnvMap<K, V> = std::collections::HashMap<K, V, Fnv>;

/// A label table: ids in first-seen order and each label's FNV; the last id, a small cache, then the map, as the C's intern
pub struct Intern { pub tab: Vec<Vec<u8>>, pub h: Vec<u64>, map: FnvMap<Vec<u8>, u32>, last: Option<u32>, cache: [u32; 512] }

impl Default for Intern {
    fn default() -> Self { Intern { tab: Vec::new(), h: Vec::new(), map: FnvMap::default(), last: None, cache: [u32::MAX; 512] } }
}

impl Intern {
    #[inline(always)]
    pub fn get(&mut self, label: &[u8]) -> u32 {
        if let Some(id) = self.last {
            if self.tab[id as usize] == label { return id; }
        }
        let n = label.len();
        let mut ck = (n as u32).wrapping_mul(2654435761);
        if n > 0 { ck ^= (label[0] as u32) << 7 ^ (label[n - 1] as u32) << 15 ^ (label[n >> 1] as u32) << 23; }
        let ck = ((ck ^ (ck >> 13)) & 511) as usize;
        let cid = self.cache[ck];
        if cid != u32::MAX && self.tab[cid as usize] == label { self.last = Some(cid); return cid; }
        let id = match self.map.get(label) {
            Some(&id) => id,
            None => {
                let id = self.tab.len() as u32;
                self.tab.push(label.to_vec());
                self.h.push(fnv(FNV0, label));
                self.map.insert(label.to_vec(), id);
                id
            }
        };
        self.last = Some(id);
        self.cache[ck] = id;
        id
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// a refusal as an Err, not an exit: bail panics on a thread with DEFER set
    fn refused<T>(f: impl FnOnce() -> T + std::panic::UnwindSafe) -> Option<String> {
        crate::DEFER.with(|d| d.set(true));
        let r = std::panic::catch_unwind(f);
        crate::DEFER.with(|d| d.set(false));
        match r {
            Ok(_) => None,
            Err(e) => Some(e.downcast::<crate::Deferred>().expect("a refusal, not a bug").0),
        }
    }

    #[test]
    fn varint_round_trips() {
        let mut vals = vec![0u64, 1, 0x7f, 0x80, 0x3fff, 0x4000, u32::MAX as u64, 1 << 49, 1 << 56, (1 << 63) - 1, 1 << 63, u64::MAX];
        for k in 0..64 { vals.push(1 << k); vals.push((1u64 << k).wrapping_sub(1)); }
        let mut buf = Vec::new();
        for &v in &vals { put(&mut buf, v); }
        let mut i = 0;
        for &v in &vals { let (got, n) = vget(&buf, i); assert_eq!(got, v); i = n; }
        assert_eq!(i, buf.len());
        // the widths: 7 bits a byte, 10 bytes for the top bit
        for (v, w) in [(0u64, 1usize), (127, 1), (128, 2), (u64::MAX >> 1, 9), (1 << 63, 10), (u64::MAX, 10)] {
            let mut b = Vec::new(); put(&mut b, v); assert_eq!(b.len(), w, "width of {v:#x}");
        }
    }

    #[test]
    fn varint_refuses_corrupt_input() {
        let msg = Some("corrupt varint".to_string());
        assert_eq!(refused(|| vget(&[], 0)), msg);
        assert_eq!(refused(|| vget(&[0x80], 0)), msg);            // a continuation and then nothing
        assert_eq!(refused(|| vget(&[0x05], 1)), msg);            // past the end
        assert_eq!(refused(|| vget(&[0xff; 11], 0)), msg);        // longer than 64 bits
        assert_eq!(refused(|| vget(&[0xff, 0xff, 0x01], 0)), None);
    }

    /// ten bytes: nine continuations, then a last byte that may carry bit 63 and nothing above it
    fn ten(top: u8) -> Vec<u8> { let mut b = vec![0x80u8; 9]; b.push(top); b }

    #[test]
    fn varint_refuses_a_tenth_byte_past_bit_63() {
        let msg = Some("corrupt varint".to_string());
        for top in 2u8..=127 {   // 2^64 and up, which a u64 cannot hold
            assert_eq!(refused(|| vget(&ten(top), 0)), msg, "tenth byte {top:#x}");
            assert_eq!(refused(|| crate::reader::skip(&ten(top), 0)), msg, "skip, tenth byte {top:#x}");
        }
        for top in [0x80u8, 0x81, 0xff] {   // a tenth byte that does not end it: an eleventh, or the cut
            let mut b = ten(top); b.push(0x00);
            assert_eq!(refused(|| vget(&b, 0)), msg, "tenth byte {top:#x}, then an eleventh");
            assert_eq!(refused(|| crate::reader::skip(&b, 0)), msg, "skip, tenth byte {top:#x}, then an eleventh");
            assert_eq!(refused(|| vget(&ten(top), 0)), msg, "tenth byte {top:#x}, then the end");
        }
        let mut b = vec![0xffu8; 9]; b.push(0x02);   // all 64 bits and one over
        assert_eq!(refused(|| vget(&b, 0)), msg);
        for n in 1..10 { assert_eq!(refused(|| vget(&ten(0x01)[..n], 0)), msg, "cut after {n} bytes"); }
    }

    #[test]
    fn varint_ten_bytes_that_fit() {
        let mut max = vec![0xffu8; 9]; max.push(0x01);
        assert_eq!(vget(&max, 0), (u64::MAX, 10));
        assert_eq!(vget(&ten(0x01), 0), (1 << 63, 10));
        assert_eq!(vget(&ten(0x00), 0), (0, 10));   // not the shortest, but it fits: still read
        assert_eq!(vget(&[0x80, 0x00], 0), (0, 2));
        for b in [max, ten(0x01), ten(0x00)] { assert_eq!(crate::reader::skip(&b, 0), 10); }
        // each width's ends, written and read back, with something after them
        for w in 1..=10u32 {
            let lo = if w == 1 { 0 } else { 1u64 << (7 * (w - 1)) };
            let hi = if w == 10 { u64::MAX } else { (1u64 << (7 * w)) - 1 };
            for v in [lo, hi] {
                let mut b = Vec::new(); put(&mut b, v); assert_eq!(b.len(), w as usize, "width of {v:#x}");
                b.push(0x05);
                assert_eq!(vget(&b, 0), (v, w as usize));
                assert_eq!(crate::reader::skip(&b, 0), w as usize);
                assert_eq!(vget(&b, w as usize), (5, w as usize + 1));
            }
        }
    }

    /// vget_r and the writer's vget_w and skip_w read what put writes as vget does; vget_r returns what vget refuses
    #[test]
    fn varint_twins_agree() {
        let mut vals = vec![0u64, 1, 0x7f, 0x80, 0x3fff, 0x4000, u32::MAX as u64, 1 << 63, u64::MAX];
        for k in 0..64 { vals.push(1 << k); vals.push((1u64 << k).wrapping_sub(1)); }
        let mut buf = Vec::new();
        for &v in &vals { put(&mut buf, v); }
        let (mut i, mut j, mut k) = (0, 0, 0);
        for &v in &vals {
            let (a, n) = vget(&buf, i);
            assert_eq!(vget_r(&buf, j).ok(), Some((a, n)));
            assert_eq!(vget_w(&buf, k), (a, n));
            assert_eq!(skip_w(&buf, k), n);
            assert_eq!(a, v);
            (i, j, k) = (n, n, n);
        }
        for b in [&[][..], &[0x80], &[0xff; 11], &ten(0x02)] {
            assert_eq!(vget_r(b, 0).map_err(|e| e.message().to_string()), Err("corrupt varint".to_string()));
        }
    }

    #[test]
    fn zigzag_round_trips() {
        for v in [0i64, 1, -1, 2, -2, 63, -64, i32::MAX as i64, i32::MIN as i64, i64::MAX, i64::MIN] {
            assert_eq!(unzz(zz(v)) as i64, v);
        }
        assert_eq!((zz(0), zz(-1), zz(1), zz(-2)), (0, 1, 2, 3));
        assert_eq!(zz(i64::MIN), u64::MAX);
    }

    #[test]
    fn putu64_at_every_digit_count() {
        let mut vals = vec![0u64, u64::MAX];
        let mut p = 1u64;
        for _ in 0..20 {
            vals.extend([p, p.saturating_sub(1), p.saturating_add(1), p.saturating_mul(5)]);
            p = p.saturating_mul(10);
        }
        for v in vals {
            let mut b = [0u8; 20];
            let n = putu64_at(&mut b, v);
            assert_eq!(&b[..n], v.to_string().as_bytes(), "{v}");
            let mut o = Vec::new(); putu64(&mut o, v); assert_eq!(o, v.to_string().as_bytes());
        }
    }

    #[test]
    fn put_row_writes_one_bed_line() {
        let mut b = b"x".to_vec();
        put_row(&mut b, b"chr1", 0, u64::MAX, b"lab");
        assert_eq!(b, format!("xchr1\t0\t{}\tlab\n", u64::MAX).as_bytes());
    }
}
