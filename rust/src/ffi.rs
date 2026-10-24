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
    // c/xz/kzlzma.c: LZMA2 with every setting pinned there, for both binaries
    fn kz_lzma_encode(src: *const c_uchar, n: usize, out: *mut c_uchar, ob: usize) -> usize;
    fn kz_lzma_decode(src: *const c_uchar, n: usize, dst: *mut c_uchar, raw: usize) -> c_int;
    fn kz_lzma_version() -> *const c_char;
    fn bz3_version() -> *const c_char;
    fn free(p: *mut c_void);
}

/// Let a write past the file-size limit return EFBIG instead of killing us,
/// which is what the C does. SIGXFSZ is 25 on Linux, the only target this
/// tree builds for: the reader already reads /proc.
pub fn ignore_sigxfsz() {
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
pub const LIBDEFLATE_VERSION: &str = env!("LIBDEFLATE_VERSION");

/// libbz3's version string, as the C `--version` prints it.
pub fn bzip3_version() -> String {
    // SAFETY: bz3_version returns a static NUL-terminated string.
    unsafe { std::ffi::CStr::from_ptr(bz3_version()).to_string_lossy().into_owned() }
}

thread_local! {
    /// Spent buffers this thread may reuse, as the C decodes every block
    /// into one blkbuf: a fresh vector a block was an mmap a block, and
    /// three times the C's page faults on an -e low decode.
    static SPARE: std::cell::RefCell<Vec<Vec<u8>>> = const { std::cell::RefCell::new(Vec::new()) };
}
/// An output buffer the callee fills: capacity only, no zero fill, so the
/// pages past what it writes are never touched. The caller sets the length
/// after the call. The codecs take `unsigned char *in` but never write it,
/// which the C relies on too, so the input is not copied either.
/// A buffer with room for cap bytes, length 0: a spare when one is big enough.
pub fn raw_buf(cap: usize) -> Vec<u8> {
    SPARE.with(|s| {
        let mut s = s.borrow_mut();
        match s.iter().position(|v| v.capacity() >= cap) {
            Some(i) => { let mut v = s.swap_remove(i); v.clear(); v }
            None => Vec::with_capacity(cap),
        }
    })
}
/// Hand back a buffer no one reads any more; four are kept, the largest.
pub fn recycle(v: Vec<u8>) {
    if v.capacity() < 1 << 17 { return; }   // below the mmap threshold: malloc reuses it already
    SPARE.with(|s| {
        let mut s = s.borrow_mut();
        s.push(v);
        if s.len() > 4 { let i = (0..s.len()).min_by_key(|&i| s[i].capacity()).unwrap(); s.swap_remove(i); }
    })
}

/// zstd, at `level`. The bound and the call are zstd's own.
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
    let mut out = raw_buf(crate::ckalloc(raw));
    zstd_decompress_into(src, &mut out, raw).then_some(out)
}

/// zstd onto the end of `out`, as paydecomp_into.
pub fn zstd_decompress_into(src: &[u8], out: &mut Vec<u8>, raw: usize) -> bool {
    crate::ckalloc(raw);
    out.reserve(raw);
    let n0 = out.len();
    // SAFETY: out has room for raw bytes past n0, and zstd is told so.
    unsafe {
        let n = zstd_sys::ZSTD_decompress(out.as_mut_ptr().add(n0) as *mut c_void, raw,
                                          src.as_ptr() as *const c_void, src.len());
        if zstd_sys::ZSTD_isError(n) != 0 || n != raw { return false; }
        out.set_len(n0 + raw);
        true
    }
}

pub fn lzma_version() -> String {
    // SAFETY: a static string from the vendored liblzma.
    unsafe { std::ffi::CStr::from_ptr(kz_lzma_version()).to_string_lossy().into_owned() }
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
    if codec & CODEC_LZMA != 0 {
        let ob = n + n / 16 + 4096;
        // SAFETY: out is ob bytes and the encoder reports how many it wrote.
        unsafe {
            let mut out = raw_buf(ob);
            let op = kz_lzma_encode(src.as_ptr(), n, out.as_mut_ptr(), ob);
            if op == 0 { return None; }
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
    let mut out = raw_buf(crate::ckalloc(raw));
    paydecomp_into(src, codec, &mut out, raw).then_some(out)
}

/// paydecomp onto the end of `out`, as the C decodes both frames into one
/// blkbuf: false when the frame does not come back to `raw` bytes.
pub fn paydecomp_into(src: &[u8], codec: i32, out: &mut Vec<u8>, raw: usize) -> bool {
    crate::ckalloc(raw);
    out.reserve(raw);
    let n0 = out.len();
    // SAFETY: out has room for raw bytes past n0; each decoder is told so and
    // checked to have written all of them before the length moves.
    unsafe {
        let dst = out.as_mut_ptr().add(n0);
        let ok = if codec & CODEC_LZMA != 0 {
            kz_lzma_decode(src.as_ptr(), src.len(), dst, raw) != 0
        } else if codec & CODEC_BZIP3 != 0 {
            let mut got = raw;
            bz3_decompress(src.as_ptr(), dst, src.len(), &mut got) == 0 && got == raw
        } else {
            if raw > u32::MAX as usize { return false; }
            let mut ul = raw as c_uint;
            let q = if codec & CODEC_ARITH != 0 {
                arith_uncompress_to(src.as_ptr() as *mut u8, src.len() as c_uint, dst, &mut ul)
            } else {
                rans_uncompress_to_4x16(src.as_ptr() as *mut u8, src.len() as c_uint, dst, &mut ul)
            };
            !q.is_null() && ul as usize == raw
        };
        if ok { out.set_len(n0 + raw); }
        ok
    }
}

/// The CRAM name tokeniser over one chunk. It overwrites every separator
/// with NUL and does not put it back, so it gets a copy, never the caller's
/// buffer: that trap has bitten this project before.
/// In place, as the C's tokc: tok3 overwrites every separator with NUL, and
/// they are put back here, so no copy of the names is made.
pub fn tok_encode(names: &mut [u8]) -> Option<Vec<u8>> {
    if names.is_empty() || names.len() > c_int::MAX as usize { return None; }
    // SAFETY: out_len is written by the callee; the returned buffer is malloc'd
    // by C and freed here after copying.
    unsafe {
        let mut ol: c_int = 0;
        let p = tok3_encode_names(names.as_mut_ptr() as *mut c_char, names.len() as c_int,
                                  9, 0, &mut ol, std::ptr::null_mut());
        for c in names.iter_mut() { if *c == 0 { *c = b'\n'; } }
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

unsafe extern "C" {
    fn dup2(a: i32, b: i32) -> i32;
    fn close(fd: i32) -> i32;
    fn unlink(p: *const c_char) -> i32;
    fn _exit(c: i32) -> !;
    fn signal(sig: c_int, handler: usize) -> usize;
}

/// -o's temporary name and the final one, as the C's OUTTMP and OUTFIN.
/// The signal handler reads the leaked C string, never the mutex.
static OUT: std::sync::Mutex<Option<(std::ffi::CString, String)>> = std::sync::Mutex::new(None);
static OUTSIG: std::sync::atomic::AtomicPtr<c_char> = std::sync::atomic::AtomicPtr::new(std::ptr::null_mut());
/// --out-dir's staged BEDs, for the handler: a leaked array of leaked C strings
static STAGE: std::sync::atomic::AtomicPtr<*mut c_char> = std::sync::atomic::AtomicPtr::new(std::ptr::null_mut());
static STAGEN: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
extern "C" fn outsig(_: c_int) {
    use std::sync::atomic::Ordering::SeqCst;
    let (p, a, n) = (OUTSIG.load(SeqCst), STAGE.load(SeqCst), STAGEN.load(SeqCst));
    // SAFETY: unlink and _exit are async-signal-safe; every pointer is leaked.
    unsafe {
        if !p.is_null() { unlink(p); }
        if !a.is_null() { for i in 0..n { unlink(*a.add(i)); } }
        _exit(1)
    }
}

/// SIGINT, SIGTERM and SIGHUP remove -o's and --out-dir's temporaries, as the C.
fn catch_signals() {
    // SAFETY: outsig only unlinks and exits.
    unsafe { for s in [1, 2, 15] { signal(s, outsig as extern "C" fn(c_int) as usize); } }
}

/// The staged BEDs a signal removes; an empty list once they are swapped in.
pub fn sig_stage(paths: &[std::path::PathBuf]) {
    use std::{os::unix::ffi::OsStrExt, sync::atomic::Ordering::SeqCst};
    let n = STAGEN.swap(0, SeqCst);
    let a = STAGE.swap(std::ptr::null_mut(), SeqCst);
    // SAFETY: a and its strings came from the leaks below, n long.
    if !a.is_null() { unsafe { for c in Box::from_raw(std::ptr::slice_from_raw_parts_mut(a, n)) { drop(std::ffi::CString::from_raw(c)); } } }
    if paths.is_empty() { return; }
    let v: Vec<*mut c_char> = paths.iter().filter_map(|p| std::ffi::CString::new(p.as_os_str().as_bytes()).ok())
        .map(|c| c.into_raw()).collect();
    let n = v.len();
    STAGE.store(Box::leak(v.into_boxed_slice()).as_mut_ptr(), SeqCst);
    STAGEN.store(n, SeqCst);
    catch_signals();
}

/// Frees the handler's copy of -o's temporary name.
fn sig_out_clear() {
    let p = OUTSIG.swap(std::ptr::null_mut(), std::sync::atomic::Ordering::SeqCst);
    // SAFETY: p came from CString::into_raw in redirect.
    if !p.is_null() { unsafe { drop(std::ffi::CString::from_raw(p)); } }
}

/// A panic is a bug, yet leaves no temporaries behind either.
pub fn panic_clean() {
    let hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |i| { crate::db::drop_tmp(); drop_out(); hook(i) }));
}

/// Removes -o's temporary file, on any failure. try_lock: a panic hook may
/// run with the lock held.
pub fn drop_out() {
    sig_out_clear();
    if let Ok(mut o) = OUT.try_lock() { if let Some((t, _)) = o.take() { unsafe { unlink(t.as_ptr()); } } }
}

/// On success: -o's temporary file renamed into place.
pub fn commit_out() {
    use std::io::Write;
    let Some((t, fin)) = OUT.lock().ok().and_then(|mut o| o.take()) else { return };
    sig_out_clear();
    let tp = t.to_string_lossy().into_owned();   // made from a &str
    // NFS reports some errors only at close
    if std::io::stdout().flush().is_err() || unsafe { close(1) } != 0 {
        let _ = std::fs::remove_file(&tp); crate::bail("write failed");
    }
    if std::fs::rename(&tp, &fin).is_err() { let _ = std::fs::remove_file(&tp); crate::bail(&format!("cannot rename into {fin}")) }
}

/// -i and -o: the file takes the place of stdin or stdout. -o to a regular
/// file is written under a temporary name, renamed in by `commit_out`; a
/// device or a FIFO is written in place.
pub fn redirect(path: &str, fd: i32) {
    use std::os::{fd::AsRawFd, unix::fs::{OpenOptionsExt, PermissionsExt}};
    let old = std::fs::symlink_metadata(path).ok();
    let mut f = None;
    if fd == 1 && old.as_ref().is_none_or(|m| m.file_type().is_file()) {
        let t = format!("{path}.tmp.{}", std::process::id());
        if let Ok(h) = std::fs::OpenOptions::new().write(true).create_new(true).mode(0o644).open(&t) {
            if let Some(m) = &old { let _ = h.set_permissions(std::fs::Permissions::from_mode(m.permissions().mode() & 0o7777)); }
            let c = std::ffi::CString::new(t).unwrap_or_else(|_| crate::bail("cannot open the output"));
            OUTSIG.store(c.clone().into_raw(), std::sync::atomic::Ordering::SeqCst);
            *OUT.lock().unwrap() = Some((c, path.to_string()));
            catch_signals();
            f = Some(Ok(h));
        }
    }
    let f = f.unwrap_or_else(|| if fd == 0 { std::fs::File::open(path) }
            else { std::fs::OpenOptions::new().write(true).create(true).truncate(true).mode(0o644).open(path) });
    if !f.is_ok_and(|f| unsafe { dup2(f.as_raw_fd(), fd) } >= 0) { crate::bail(&format!("cannot open {path}")) }
}

// glibc's byte search: on short fields 2.4x fewer instructions than the
// memchr crate's dispatch layers, and the same code the C calls
unsafe extern "C" {
    #[link_name = "memchr"] fn c_memchr(s: *const u8, c: i32, n: usize) -> *const u8;
    #[link_name = "memrchr"] fn c_memrchr(s: *const u8, c: i32, n: usize) -> *const u8;
    #[link_name = "strchrnul"] fn c_strchrnul(s: *const u8, c: i32) -> *const u8;
}
fn at(s: &[u8], p: *const u8) -> Option<usize> { if p.is_null() { None } else { Some(p as usize - s.as_ptr() as usize) } }
/// The first `b` in `s`.
pub fn find(b: u8, s: &[u8]) -> Option<usize> { at(s, unsafe { c_memchr(s.as_ptr(), b as i32, s.len()) }) }
/// The last `b` in `s`.
pub fn rfind(b: u8, s: &[u8]) -> Option<usize> { at(s, unsafe { c_memrchr(s.as_ptr(), b as i32, s.len()) }) }
/// The first newline or NUL in `s`, which must end in a newline.
pub fn line_end(s: &[u8]) -> usize {
    assert_eq!(s.last(), Some(&b'\n'));
    unsafe { c_strchrnul(s.as_ptr(), b'\n' as i32) as usize - s.as_ptr() as usize }
}

// BGZF blocks state their size: libdeflate's whole-buffer inflate and
// PCLMUL CRC32, vendored in c/, as the C's
#[repr(C)] struct LdDec { _p: [u8; 0] }
unsafe extern "C" {
    fn libdeflate_alloc_decompressor() -> *mut LdDec;
    fn libdeflate_free_decompressor(d: *mut LdDec);
    fn libdeflate_deflate_decompress(d: *mut LdDec, i: *const u8, il: usize, o: *mut u8, ol: usize, got: *mut usize) -> i32;
    fn libdeflate_crc32(crc: u32, b: *const u8, n: usize) -> u32;
}
pub struct Inflater(*mut LdDec);
unsafe impl Send for Inflater {}
impl Inflater {
    pub fn new() -> Inflater { let d = unsafe { libdeflate_alloc_decompressor() }; assert!(!d.is_null(), "out of memory"); Inflater(d) }
    /// Inflates raw deflate `i` into all of `o`: true when it fills `o` exactly.
    pub fn inflate(&mut self, i: &[u8], o: &mut [u8]) -> bool {
        let mut got = 0usize;
        (unsafe { libdeflate_deflate_decompress(self.0, i.as_ptr(), i.len(), o.as_mut_ptr(), o.len(), &mut got) }) == 0 && got == o.len()
    }
}
impl Drop for Inflater { fn drop(&mut self) { unsafe { libdeflate_free_decompressor(self.0) } } }
pub fn crc32(b: &[u8]) -> u32 { unsafe { libdeflate_crc32(0, b.as_ptr(), b.len()) } }

unsafe extern "C" { fn lseek(fd: i32, off: i64, whence: i32) -> i64; }
/// Where stdin's file offset is; None for a pipe.
pub fn stdin_at() -> Option<u64> { u64::try_from(unsafe { lseek(0, 0, 1) }).ok() }
