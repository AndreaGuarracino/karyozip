//! The only unsafe in the crate, besides the reader's mmap. Every raw
//! pointer and every length check lives here; callers see slices and
//! `Option<Vec<u8>>`.
//!
//! The codecs are the C ones from ../c. htscodecs has no Rust port, and
//! bzip3's crate lacks the two fixes our copy carries, one of which loses
//! half the data on an input that is an exact multiple of the block size.

use std::os::raw::{c_char, c_int, c_uchar, c_uint, c_void};

pub const CODEC_ARITH: i32 = 0x10000;
pub const CODEC_BZIP3: i32 = 0x20000;

unsafe extern "C" {
    fn rans_compress_bound_4x16(size: c_uint, order: c_int) -> c_uint;
    fn rans_compress_to_4x16(inp: *mut c_uchar, in_size: c_uint,
                             out: *mut c_uchar, out_size: *mut c_uint,
                             order: c_int) -> *mut c_uchar;
    fn rans_uncompress_to_4x16(inp: *mut c_uchar, in_size: c_uint,
                               out: *mut c_uchar, out_size: *mut c_uint) -> *mut c_uchar;
    fn arith_compress_to(inp: *mut c_uchar, in_size: c_uint,
                         out: *mut c_uchar, out_size: *mut c_uint,
                         order: c_int) -> *mut c_uchar;
    fn arith_uncompress_to(inp: *mut c_uchar, in_size: c_uint,
                           out: *mut c_uchar, out_size: *mut c_uint) -> *mut c_uchar;
    fn tok3_encode_names(blk: *mut c_char, len: c_int, level: c_int, use_arith: c_int,
                         out_len: *mut c_int, last_start_p: *mut c_int) -> *mut c_uchar;
    fn tok3_decode_names(inp: *mut c_uchar, sz: c_uint, out_len: *mut c_uint) -> *mut c_uchar;
    fn bz3_bound(size: usize) -> usize;
    fn bz3_compress(block_size: c_uint, inp: *const c_uchar, out: *mut c_uchar,
                    in_size: usize, out_size: *mut usize) -> c_int;
    fn bz3_decompress(inp: *const c_uchar, out: *mut c_uchar,
                      in_size: usize, out_size: *mut usize) -> c_int;
    fn bz3_version() -> *const c_char;
    fn free(p: *mut c_void);
}

/// libbz3's version string, as the C `--version` prints it.
pub fn bzip3_version() -> String {
    // SAFETY: bz3_version returns a static NUL-terminated string.
    unsafe { std::ffi::CStr::from_ptr(bz3_version()).to_string_lossy().into_owned() }
}

/// zstd, at `level`. The bound and the call are zstd's own.
pub fn zstd(src: &[u8], level: i32) -> Vec<u8> {
    // SAFETY: the bound is what zstd asks for, and both slices are live.
    unsafe {
        let cap = zstd_sys::ZSTD_compressBound(src.len());
        let mut out = vec![0u8; cap];
        let n = zstd_sys::ZSTD_compress(out.as_mut_ptr() as *mut c_void, cap,
                                        src.as_ptr() as *const c_void, src.len(), level);
        assert!(zstd_sys::ZSTD_isError(n) == 0, "zstd compress failed");
        out.truncate(n);
        out
    }
}

/// The decoded size a zstd frame declares, or None when it does not.
pub fn zstd_frame_size(src: &[u8]) -> Option<usize> {
    // SAFETY: reads the frame header of a live slice.
    unsafe {
        let n = zstd_sys::ZSTD_getFrameContentSize(src.as_ptr() as *const c_void, src.len());
        if n == u64::MAX || n == u64::MAX - 1 { None } else { Some(n as usize) }
    }
}

pub fn zstd_decompress(src: &[u8], raw: usize) -> Option<Vec<u8>> {
    // SAFETY: out is sized to the caller's expected length.
    unsafe {
        let mut out = vec![0u8; raw];
        let n = zstd_sys::ZSTD_decompress(out.as_mut_ptr() as *mut c_void, raw,
                                          src.as_ptr() as *const c_void, src.len());
        if zstd_sys::ZSTD_isError(n) != 0 || n != raw { return None; }
        Some(out)
    }
}

pub fn zstd_version() -> String {
    // SAFETY: a static string from the linked libzstd.
    unsafe { std::ffi::CStr::from_ptr(zstd_sys::ZSTD_versionString()).to_string_lossy().into_owned() }
}

/// One payload codec, chosen by `codec`, exactly as the C `paycomp` does.
/// Returns None when the codec refuses, which the caller treats as "did not
/// win" rather than as an error.
pub fn paycomp(src: &[u8], codec: i32) -> Option<Vec<u8>> {
    let n = src.len();
    if codec & CODEC_BZIP3 != 0 {
        // SAFETY: out is bz3_bound(n), which is what bz3_compress needs.
        unsafe {
            let bb = bz3_bound(n);
            let mut out = vec![0u8; bb];
            let mut got = bb;
            let bs = ((codec & 0xffff) as c_uint) << 20;
            if bz3_compress(bs, src.as_ptr(), out.as_mut_ptr(), n, &mut got) != 0 { return None; }
            out.truncate(got);
            return Some(out);
        }
    }
    if n >= 0x7fff_ffff { return None; }   // htscodecs is 32-bit throughout
    if codec & CODEC_ARITH != 0 {
        let ob = (n as u64) * 2 + 1024;
        if ob > u32::MAX as u64 { return None; }   // the C wrapped here once
        // SAFETY: ol is the buffer's true length and is updated in place.
        unsafe {
            let mut out = vec![0u8; ob as usize];
            let mut ol = ob as c_uint;
            let mut inp = src.to_vec();
            let r = arith_compress_to(inp.as_mut_ptr(), n as c_uint,
                                      out.as_mut_ptr(), &mut ol, codec & 0xffff);
            if r.is_null() { return None; }
            out.truncate(ol as usize);
            Some(out)
        }
    } else {
        // SAFETY: the bound comes from rans itself.
        unsafe {
            let mut ol = rans_compress_bound_4x16(n as c_uint, codec);
            let mut out = vec![0u8; ol as usize];
            let mut inp = src.to_vec();
            let r = rans_compress_to_4x16(inp.as_mut_ptr(), n as c_uint,
                                          out.as_mut_ptr(), &mut ol, codec);
            if r.is_null() { return None; }
            out.truncate(ol as usize);
            Some(out)
        }
    }
}

/// The matching decompressor. `raw` is the decoded length the file claims.
pub fn paydecomp(src: &[u8], codec: i32, raw: usize) -> Option<Vec<u8>> {
    if codec & CODEC_BZIP3 != 0 {
        // SAFETY: out is the claimed decoded length; got is checked against it.
        unsafe {
            let mut out = vec![0u8; raw];
            let mut got = raw;
            if bz3_decompress(src.as_ptr(), out.as_mut_ptr(), src.len(), &mut got) != 0
                || got != raw { return None; }
            return Some(out);
        }
    }
    if raw > u32::MAX as usize { return None; }
    // SAFETY: out is the claimed length and ul is checked against it.
    unsafe {
        let mut out = vec![0u8; raw];
        let mut ul = raw as c_uint;
        let mut inp = src.to_vec();
        let q = if codec & CODEC_ARITH != 0 {
            arith_uncompress_to(inp.as_mut_ptr(), src.len() as c_uint, out.as_mut_ptr(), &mut ul)
        } else {
            rans_uncompress_to_4x16(inp.as_mut_ptr(), src.len() as c_uint, out.as_mut_ptr(), &mut ul)
        };
        if q.is_null() || ul as usize != raw { return None; }
        Some(out)
    }
}

/// The CRAM name tokeniser over one chunk. It overwrites every separator
/// with NUL and does not put it back, so it gets a copy, never the caller's
/// buffer: that trap has bitten this project before.
pub fn tok_encode(names: &[u8]) -> Option<Vec<u8>> {
    if names.is_empty() || names.len() > c_int::MAX as usize { return None; }
    let mut copy = names.to_vec();
    // SAFETY: out_len is written by the callee; the returned buffer is malloc'd
    // by C and freed here after copying.
    unsafe {
        let mut ol: c_int = 0;
        let p = tok3_encode_names(copy.as_mut_ptr() as *mut c_char, copy.len() as c_int,
                                  9, 0, &mut ol, std::ptr::null_mut());
        if p.is_null() || ol <= 0 { return None; }
        let v = std::slice::from_raw_parts(p, ol as usize).to_vec();
        free(p as *mut c_void);
        Some(v)
    }
}

pub fn tok_decode(src: &[u8]) -> Option<Vec<u8>> {
    // SAFETY: as above; the callee reports the decoded length.
    unsafe {
        let mut ol: c_uint = 0;
        let mut inp = src.to_vec();
        let p = tok3_decode_names(inp.as_mut_ptr(), src.len() as c_uint, &mut ol);
        if p.is_null() { return None; }
        let v = std::slice::from_raw_parts(p, ol as usize).to_vec();
        free(p as *mut c_void);
        Some(v)
    }
}
