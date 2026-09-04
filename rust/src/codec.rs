//! Choosing a codec for the two frames, exactly as c/karyozip.c does.
//! The choice decides the stored bytes, so every rule here is the C's rule.

use crate::ffi;

/// The eight candidates the default search tries, beside zstd. Codec 0 means
/// zstd, so bare rANS order 0 is spelled 0x100.
pub const CAND_FAST: [i32; 8] = [
    0x100, 1, 65, 521,
    ffi::CODEC_ARITH, ffi::CODEC_ARITH | 1,
    ffi::CODEC_BZIP3 | 16, ffi::CODEC_BZIP3 | 64,
];

/// bzip3 decodes at about 8 MB/s against zstd's ~1000, so it is the wrong
/// pick for a file a query will seek into.
fn fast_enough(codec: i32) -> bool {
    codec & ffi::CODEC_BZIP3 == 0
}

pub struct Pick {
    pub id: i32,
    pub len: i32,
}

/// `samples` are whole blocks: (lens frame, ids frame).
pub fn pick(samples: &[(&[u8], &[u8])], seekable: bool, level: i32,
            forced_id: Option<i32>, forced_len: Option<i32>) -> Pick {
    if let (Some(i), Some(l)) = (forced_id, forced_len) {
        return Pick { id: i, len: l };
    }
    // 4 MB per block when seekable: a smaller sample cannot show what a BWT
    // codec does to a 100 MB stream.
    let scap = if seekable { 4 << 20 } else { 96 << 20 };
    let mut maxstream = 0usize;
    for (l, i) in samples {
        maxstream = maxstream.max(l.len()).max(i.len());
    }
    let big = maxstream >= 0x7fff_ffff;

    let mut zi = 0usize;
    let mut zl = 0usize;
    let mut ci = [0usize; CAND_FAST.len()];
    let mut cl = [0usize; CAND_FAST.len()];
    const BAD: usize = usize::MAX / 2;

    for (lens, ids) in samples {
        let nl = lens.len().min(scap);
        let ni = ids.len().min(scap);
        if ni > 0 { zi += ffi::zstd(&ids[..ni], level).len(); }
        if nl > 0 { zl += ffi::zstd(&lens[..nl], level).len(); }
        for (k, &c) in CAND_FAST.iter().enumerate() {
            // Size alone does not qualify a codec: it must decode back.
            ci[k] = accumulate(ci[k], &ids[..ni], c, BAD);
            cl[k] = accumulate(cl[k], &lens[..nl], c, BAD);
        }
    }

    let (mut wi, mut wl) = (0i32, 0i32);
    let (mut bi, mut bl) = (zi, zl);
    // The 0.5% margin, in the C's exact arithmetic.
    let need_i = (0.995_f64 * zi as f64) as usize;
    let need_l = (0.995_f64 * zl as f64) as usize;
    for (k, &c) in CAND_FAST.iter().enumerate() {
        let ok = (!seekable || fast_enough(c)) && (!big || c & ffi::CODEC_BZIP3 != 0);
        if zi > 0 && ci[k] < bi && ci[k] < need_i && ok { bi = ci[k]; wi = c; }
        if zl > 0 && cl[k] < bl && cl[k] < need_l && ok { bl = cl[k]; wl = c; }
    }
    Pick { id: forced_id.unwrap_or(wi), len: forced_len.unwrap_or(wl) }
}

fn accumulate(acc: usize, src: &[u8], codec: i32, bad: usize) -> usize {
    if src.is_empty() { return acc; }
    if acc >= bad { return bad; }
    match ffi::paycomp(src, codec) {
        Some(t) => match ffi::paydecomp(&t, codec, src.len()) {
            Some(v) if v == src => acc + t.len(),
            _ => bad,
        },
        None => bad,
    }
}

/// The lengths frame: a forced codec, else zstd.
pub fn comp_lens(src: &[u8], codec: i32, level: i32) -> Option<Vec<u8>> {
    if codec != 0 { ffi::paycomp(src, codec) } else { Some(ffi::zstd(src, level)) }
}

/// The ids frame, same rule.
pub fn comp_ids(src: &[u8], codec: i32, level: i32) -> Option<Vec<u8>> {
    if codec != 0 { ffi::paycomp(src, codec) } else { Some(ffi::zstd(src, level)) }
}
