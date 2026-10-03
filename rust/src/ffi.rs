//! The only unsafe in the crate. Every raw
//! pointer and every length check lives here; callers see slices and
//! `Option<Vec<u8>>`.
//!
//! The codecs are the C ones from ../c. htscodecs has no Rust port, and
//! bzip3's crate lacks the two fixes our copy carries, one of which loses
//! half the data on an input that is an exact multiple of the block size.

use std::os::raw::{c_char, c_int, c_uchar, c_uint, c_void};

pub const CODEC_ARITH: i32 = 0x10000;
pub const CODEC_BZIP3: i32 = 0x20000;
pub const CODEC_LZMA: i32 = 0x40000;

/// Stores no size, so a size check cannot be made.
pub fn nosz(c: i32) -> bool {
    let u = c as u32;
    u != 0 && u & CODEC_BZIP3 as u32 == 0 && u & 0x10 != 0
}

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

/// Let a write past the file-size limit return EFBIG instead of killing us,
/// which is what the C does. SIGXFSZ is 25 on Linux, the only target this
/// tree builds for: the reader already reads /proc.
pub fn ignore_sigxfsz() {
    unsafe extern "C" { fn signal(sig: c_int, handler: usize) -> usize; }
    // SAFETY: SIG_IGN is 1, and signal takes no argument we own.
    unsafe { signal(25, 1); }
}

/// glibc raises its mmap threshold to the size of the first large mapped
/// chunk it frees, after which every big buffer comes from a per-thread
/// heap and stays resident once freed: 55 MB against the C's 37 on a
/// 300000-row encode. Pinned at the default, so big buffers are mapped and
/// unmapped, as the C's allocation pattern leaves them.
/// The OpenMP default team, which libbz3's block loop takes: -t, as the C.
pub fn omp_threads(n: usize) {
    unsafe extern "C" { fn omp_set_num_threads(n: c_int); }
    // SAFETY: takes one int and only sets libgomp's default team size.
    unsafe { omp_set_num_threads(n as c_int); }
}

pub fn pin_mmap_threshold() {
    unsafe extern "C" { fn mallopt(param: c_int, value: c_int) -> c_int; }
    // SAFETY: M_MMAP_THRESHOLD is -3; mallopt takes two ints.
    unsafe { mallopt(-3, 128 << 10); }
}

/// Whether fds 0, 1 and 2 were open when the process started. The runtime
/// reopens a closed one on /dev/null before main, so the check runs from
/// .init_array, ahead of it, as the C's fcntl check at the top of main.
pub static STDFDS_OPEN: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(true);
#[used]
#[unsafe(link_section = ".init_array")]
static CHECK_STDFDS: extern "C" fn() = check_stdfds;
extern "C" fn check_stdfds() {
    unsafe extern "C" { fn fcntl(fd: c_int, cmd: c_int, ...) -> c_int; }
    // SAFETY: F_GETFD is 1 and takes no argument.
    let open = (0..3).all(|fd| unsafe { fcntl(fd, 1) } >= 0);
    STDFDS_OPEN.store(open, std::sync::atomic::Ordering::Relaxed);
}

/// The vendored htscodecs, from its own version.h at build time.
pub const HTSCODECS_VERSION: &str = env!("HTSCODECS_VERSION");

/// libbz3's version string, as the C `--version` prints it.
pub fn bzip3_version() -> String {
    // SAFETY: bz3_version returns a static NUL-terminated string.
    unsafe { std::ffi::CStr::from_ptr(bz3_version()).to_string_lossy().into_owned() }
}

/// zstd, at `level`. The bound and the call are zstd's own.
/// An output buffer the callee fills: capacity only, no zero fill, so the
/// pages past what it writes are never touched. The caller sets the length
/// after the call. The codecs take `unsigned char *in` but never write it,
/// which the C relies on too, so the input is not copied either.
fn raw_buf(cap: usize) -> Vec<u8> { Vec::with_capacity(cap) }

pub fn zstd(src: &[u8], level: i32) -> Vec<u8> {
    // SAFETY: the bound is what zstd asks for, and both slices are live.
    unsafe {
        let cap = zstd_sys::ZSTD_compressBound(src.len());
        let mut out = raw_buf(cap);
        let n = zstd_sys::ZSTD_compress(out.as_mut_ptr() as *mut c_void, cap,
                                        src.as_ptr() as *const c_void, src.len(), level);
        assert!(zstd_sys::ZSTD_isError(n) == 0, "zstd compress failed");
        out.set_len(n);
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
        let mut out = raw_buf(crate::ckalloc(raw));
        let n = zstd_sys::ZSTD_decompress(out.as_mut_ptr() as *mut c_void, raw,
                                          src.as_ptr() as *const c_void, src.len());
        if zstd_sys::ZSTD_isError(n) != 0 || n != raw { return None; }
        out.set_len(raw);
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
/// every field pinned: a preset that moved between liblzma releases would
/// change the bytes, and the C and the Rust must write the same ones
fn lzma_filt(o: &mut lzma_sys::lzma_options_lzma) -> [lzma_sys::lzma_filter; 2] {
    o.dict_size = 64 << 20; o.lc = 3; o.lp = 0; o.pb = 0;
    o.mode = lzma_sys::LZMA_MODE_NORMAL; o.nice_len = 273;
    o.mf = lzma_sys::LZMA_MF_BT4; o.depth = 0;
    [lzma_sys::lzma_filter { id: lzma_sys::LZMA_FILTER_LZMA2,
                             options: o as *mut _ as *mut std::ffi::c_void },
     lzma_sys::lzma_filter { id: lzma_sys::LZMA_VLI_UNKNOWN, options: std::ptr::null_mut() }]
}

pub fn paycomp(src: &[u8], codec: i32) -> Option<Vec<u8>> {
    let n = src.len();
    if codec & CODEC_LZMA != 0 {
        let mut o: lzma_sys::lzma_options_lzma = unsafe { std::mem::zeroed() };
        let f = lzma_filt(&mut o);
        let ob = n + n / 16 + 4096;
        // SAFETY: out is ob bytes and op is where the encoder stopped.
        unsafe {
            let mut out = raw_buf(ob);
            let mut op = 0usize;
            if lzma_sys::lzma_raw_buffer_encode(f.as_ptr(), std::ptr::null(),
                    src.as_ptr(), n, out.as_mut_ptr(), &mut op, ob) != lzma_sys::LZMA_OK {
                return None;
            }
            out.set_len(op);
            return Some(out);
        }
    }
    if codec & CODEC_BZIP3 != 0 {
        // SAFETY: out is bz3_bound(n), which is what bz3_compress needs.
        unsafe {
            let bb = bz3_bound(n);
            let mut out = raw_buf(bb);
            let mut got = bb;
            let bs = ((codec & 0xffff) as c_uint) << 20;
            if bz3_compress(bs, src.as_ptr(), out.as_mut_ptr(), n, &mut got) != 0 { return None; }
            out.set_len(got);
            return Some(out);
        }
    }
    if n >= 0x7fff_ffff { return None; }   // htscodecs is 32-bit throughout
    if codec & CODEC_ARITH != 0 {
        let ob = (n as u64) * 2 + 1024;
        if ob > u32::MAX as u64 { return None; }   // the C wrapped here once
        // SAFETY: ol is the buffer's true length and is updated in place.
        unsafe {
            let mut out = raw_buf(ob as usize);
            let mut ol = ob as c_uint;
            let r = arith_compress_to(src.as_ptr() as *mut u8, n as c_uint,
                                      out.as_mut_ptr(), &mut ol, codec & 0xffff);
            if r.is_null() { return None; }
            out.set_len(ol as usize);
            Some(out)
        }
    } else {
        // SAFETY: the bound comes from rans itself.
        unsafe {
            let mut ol = rans_compress_bound_4x16(n as c_uint, codec);
            let mut out = raw_buf(ol as usize);
            let r = rans_compress_to_4x16(src.as_ptr() as *mut u8, n as c_uint,
                                          out.as_mut_ptr(), &mut ol, codec);
            if r.is_null() { return None; }
            out.set_len(ol as usize);
            Some(out)
        }
    }
}

/// The matching decompressor. `raw` is the decoded length the file claims.
pub fn paydecomp(src: &[u8], codec: i32, raw: usize) -> Option<Vec<u8>> {
    crate::ckalloc(raw);
    if codec & CODEC_LZMA != 0 {
        let mut o: lzma_sys::lzma_options_lzma = unsafe { std::mem::zeroed() };
        let f = lzma_filt(&mut o);
        // SAFETY: out is the claimed length and op is checked against it.
        unsafe {
            let mut out = raw_buf(raw);
            let (mut ip, mut op) = (0usize, 0usize);
            if lzma_sys::lzma_raw_buffer_decode(f.as_ptr(), std::ptr::null(),
                    src.as_ptr(), &mut ip, src.len(), out.as_mut_ptr(), &mut op, raw)
                    != lzma_sys::LZMA_OK || op != raw { return None; }
            out.set_len(raw);
            return Some(out);
        }
    }
    if codec & CODEC_BZIP3 != 0 {
        // SAFETY: out is the claimed decoded length; got is checked against it.
        unsafe {
            let mut out = raw_buf(raw);
            let mut got = raw;
            if bz3_decompress(src.as_ptr(), out.as_mut_ptr(), src.len(), &mut got) != 0
                || got != raw { return None; }
            out.set_len(raw);
            return Some(out);
        }
    }
    if raw > u32::MAX as usize { return None; }
    // SAFETY: out is the claimed length and ul is checked against it.
    unsafe {
        let mut out = raw_buf(raw);
        let mut ul = raw as c_uint;
        let q = if codec & CODEC_ARITH != 0 {
            arith_uncompress_to(src.as_ptr() as *mut u8, src.len() as c_uint, out.as_mut_ptr(), &mut ul)
        } else {
            rans_uncompress_to_4x16(src.as_ptr() as *mut u8, src.len() as c_uint, out.as_mut_ptr(), &mut ul)
        };
        if q.is_null() || ul as usize != raw { return None; }
        out.set_len(raw);
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
        let p = tok3_decode_names(src.as_ptr() as *mut u8, src.len() as c_uint, &mut ol);
        if p.is_null() { return None; }
        let v = std::slice::from_raw_parts(p, ol as usize).to_vec();
        free(p as *mut c_void);
        Some(v)
    }
}

unsafe extern "C" { fn dup2(a: i32, b: i32) -> i32; }

/// -i and -o: the file takes the place of stdin or stdout
pub fn redirect(path: &str, fd: i32) {
    use std::os::{fd::AsRawFd, unix::fs::OpenOptionsExt};
    let f = if fd == 0 { std::fs::File::open(path) }
            else { std::fs::OpenOptions::new().write(true).create(true).truncate(true).mode(0o644).open(path) };
    if !f.is_ok_and(|f| unsafe { dup2(f.as_raw_fd(), fd) } >= 0) { crate::bail(&format!("cannot open {path}")) }
}
