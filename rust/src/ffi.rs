//! The FFI: raw pointers and length checks for the codecs, which are the C ones from ../c.
//! The only other unsafe is unchecked indexing in kst.rs (put_row, putu64_at) and smooth.rs (smooth_one, emit).

use std::os::raw::{c_char, c_int, c_uchar, c_uint, c_void};
use std::os::unix::ffi::{OsStrExt, OsStringExt};

pub const CODEC_ARITH: i32 = 0x10000;
pub const CODEC_BZIP3: i32 = 0x20000;
pub const CODEC_LZMA: i32 = 0x40000;
/// 0x80000|MB: libzpaq level 5 on pieces of MB, -e high only
pub const CODEC_ZPAQ: i32 = 0x80000;

/// Stores no size, so a size check cannot be made
pub fn nosz(c: i32) -> bool {
    let u = c as u32;
    u != 0 && u & (CODEC_BZIP3 | CODEC_ZPAQ) as u32 == 0 && u & 0x10 != 0
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
    fn bz3_compress_cuts(block_size: c_uint, inp: *const c_uchar, out: *mut c_uchar,
                         in_size: usize, out_size: *mut usize, cut: *const c_uint, ncut: c_uint) -> c_int;
    fn bz3_decompress(inp: *const c_uchar, out: *mut c_uchar,
                      in_size: usize, out_size: *mut usize) -> c_int;
    // c/xz/kzlzma.c: LZMA2 with every setting pinned there, for both binaries
    fn kz_lzma_encode(src: *const c_uchar, n: usize, out: *mut c_uchar, ob: usize) -> usize;
    fn kz_lzma_decode(src: *const c_uchar, n: usize, dst: *mut c_uchar, raw: usize) -> c_int;
    fn kz_lzma_version() -> *const c_char;
    // c/zpaq/kzzpaq.cpp: libzpaq on the caller's pieces, for both binaries
    fn kz_zpaq_encode(inp: *const c_uchar, n: usize, cut: *const c_uint, ncut: c_uint, method: *const c_char,
                      out: *mut c_uchar, cap: usize) -> usize;
    fn kz_zpaq_decode(inp: *const c_uchar, n: usize, dst: *mut c_uchar, raw: usize) -> c_int;
    fn bz3_version() -> *const c_char;
    fn free(p: *mut c_void);
}

/// A write past the file-size limit returns EFBIG instead of killing us, as the C; SIGXFSZ is 25 on Linux
pub fn ignore_sigxfsz() {
    // SAFETY: SIG_IGN is 1, and signal takes no argument we own
    unsafe { signal(25, 1); }
}

/// glibc's mmap threshold pinned at the default, so big buffers are mapped and unmapped as the C's; OpenMP's team is -t
pub fn omp_threads(n: usize) {
    unsafe extern "C" { fn omp_set_num_threads(n: c_int); }
    // SAFETY: takes one int and only sets libgomp's default team size
    unsafe { omp_set_num_threads(n as c_int); }
}

pub fn pin_mmap_threshold() {
    unsafe extern "C" { fn mallopt(param: c_int, value: c_int) -> c_int; }
    // SAFETY: M_MMAP_THRESHOLD is -3; mallopt takes two ints
    unsafe { mallopt(-3, 128 << 10); }
}

/// Whether fds 0, 1 and 2 were open when the process started
pub static STDFDS_OPEN: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(true);
#[used]
#[unsafe(link_section = ".init_array")]
static CHECK_STDFDS: extern "C" fn() = check_stdfds;
extern "C" fn check_stdfds() {
    unsafe extern "C" { fn fcntl(fd: c_int, cmd: c_int, ...) -> c_int; }
    // SAFETY: F_GETFD is 1 and takes no argument
    let open = (0..3).all(|fd| unsafe { fcntl(fd, 1) } >= 0);
    STDFDS_OPEN.store(open, std::sync::atomic::Ordering::Relaxed);
}

/// The vendored htscodecs, from its own version.h at build time
pub const HTSCODECS_VERSION: &str = env!("HTSCODECS_VERSION");
pub const LIBDEFLATE_VERSION: &str = env!("LIBDEFLATE_VERSION");

/// libbz3's version string, as the C `--version` prints it
pub fn bzip3_version() -> String {
    // SAFETY: bz3_version returns a static NUL-terminated string
    unsafe { std::ffi::CStr::from_ptr(bz3_version()).to_string_lossy().into_owned() }
}

thread_local! {
    /// Spent buffers this thread may reuse, as the C decodes every block into one blkbuf
    static SPARE: std::cell::RefCell<Vec<Vec<u8>>> = const { std::cell::RefCell::new(Vec::new()) };
}
/// A buffer with room for cap bytes, length 0, a spare when one is big enough; no zero fill, the callee fills it
pub fn raw_buf(cap: usize) -> Vec<u8> {
    SPARE.with(|s| {
        let mut s = s.borrow_mut();
        match s.iter().position(|v| v.capacity() >= cap) {
            Some(i) => { let mut v = s.swap_remove(i); v.clear(); v }
            None => Vec::with_capacity(cap),
        }
    })
}
/// Hand back a buffer no one reads any more; four are kept, the largest
pub fn recycle(v: Vec<u8>) {
    if v.capacity() < 1 << 17 { return; }   // below the mmap threshold: malloc reuses it already
    SPARE.with(|s| {
        let mut s = s.borrow_mut();
        s.push(v);
        if s.len() > 4 { let i = (0..s.len()).min_by_key(|&i| s[i].capacity()).unwrap(); s.swap_remove(i); }
    })
}

/// zstd at `level`, with zstd's own bound
pub fn zstd(src: &[u8], level: i32) -> Vec<u8> {
    // SAFETY: the bound is what zstd asks for, and both slices are live
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

/// The decoded size a zstd frame declares, or None when it does not
pub fn zstd_frame_size(src: &[u8]) -> Option<usize> {
    // SAFETY: reads the frame header of a live slice
    unsafe {
        let n = zstd_sys::ZSTD_getFrameContentSize(src.as_ptr() as *const c_void, src.len());
        if n == u64::MAX || n == u64::MAX - 1 { None } else { Some(n as usize) }
    }
}

pub fn zstd_decompress(src: &[u8], raw: usize) -> Option<Vec<u8>> {
    let mut out = raw_buf(crate::ckalloc(raw));
    zstd_decompress_into(src, &mut out, raw).then_some(out)
}

/// zstd onto the end of `out`, as paydecomp_into
pub fn zstd_decompress_into(src: &[u8], out: &mut Vec<u8>, raw: usize) -> bool {
    crate::ckalloc(raw);
    out.reserve(raw);
    let n0 = out.len();
    // SAFETY: out has room for raw bytes past n0, and zstd is told so
    unsafe {
        let n = zstd_sys::ZSTD_decompress(out.as_mut_ptr().add(n0) as *mut c_void, raw,
                                          src.as_ptr() as *const c_void, src.len());
        if zstd_sys::ZSTD_isError(n) != 0 || n != raw { return false; }
        out.set_len(n0 + raw);
        true
    }
}

pub fn lzma_version() -> String {
    // SAFETY: a static string from the vendored liblzma
    unsafe { std::ffi::CStr::from_ptr(kz_lzma_version()).to_string_lossy().into_owned() }
}

pub fn zstd_version() -> String {
    // SAFETY: a static string from the linked libzstd
    unsafe { std::ffi::CStr::from_ptr(zstd_sys::ZSTD_versionString()).to_string_lossy().into_owned() }
}

/// The C's strerror(errno) words, without Rust's " (os error N)"
pub fn strerror(e: &std::io::Error) -> String {
    unsafe extern "C" { fn strerror(n: c_int) -> *const c_char; }
    let Some(n) = e.raw_os_error() else { return e.to_string() };
    // SAFETY: strerror returns a NUL-terminated string, read before any other call can reuse it
    unsafe { std::ffi::CStr::from_ptr(strerror(n)).to_string_lossy().into_owned() }
}

/// bzip3 with the pieces' sizes given, each at most the block size, summing to the source: the C's bz3_compress_cuts
pub fn bz3_cuts(src: &[u8], codec: i32, cut: &[u32]) -> Option<Vec<u8>> {
    let n = src.len();
    // SAFETY: out is bz3_bound(n), cut holds cut.len() sizes; the C checks them
    unsafe {
        let bb = bz3_bound(n);
        let mut out = raw_buf(bb);
        let mut got = bb;
        let bs = ((codec & 0xffff) as c_uint) << 20;
        if bz3_compress_cuts(bs, src.as_ptr(), out.as_mut_ptr(), n, &mut got, cut.as_ptr(), cut.len() as c_uint) != 0 { return None; }
        out.set_len(got);
        Some(out)
    }
}

/// n bytes in pieces of b, the last shorter, as the C's fixed zpaq pieces
pub fn zpaq_pieces(n: usize, b: usize) -> Vec<u32> { (0..n).step_by(b.max(1)).map(|o| (n - o).min(b) as u32).collect() }

/// libzpaq level 5, or the interleaved model when ilv, on pieces summing to the source: the C's zpaq branch of paycut
pub fn zpaq_cuts(src: &[u8], cut: &[u32], ilv: bool) -> Option<Vec<u8>> {
    let cap = src.len() + src.len() / 8 + 4096 + 1024 * cut.len();
    // SAFETY: out holds cap bytes, cut holds cut.len() sizes; the C checks them
    unsafe {
        let mut out = raw_buf(cap);
        let z = kz_zpaq_encode(src.as_ptr(), src.len(), cut.as_ptr(), cut.len() as c_uint, if ilv { c"ilv" } else { c"5" }.as_ptr(), out.as_mut_ptr(), cap);
        if z == 0 { return None; }
        out.set_len(z);
        Some(out)
    }
}

/// One payload codec, as the C's paycomp; None when the codec refuses, which is "did not win"
pub fn paycomp(src: &[u8], codec: i32) -> Option<Vec<u8>> {
    let n = src.len();
    if codec & CODEC_ZPAQ != 0 {   // pieces of MB, as the C's
        return zpaq_cuts(src, &zpaq_pieces(n, ((codec & 0xffff) as usize) << 20), false);
    }
    if codec & CODEC_LZMA != 0 {
        let ob = n + n / 16 + 4096;
        // SAFETY: out is ob bytes and the encoder reports how many it wrote
        unsafe {
            let mut out = raw_buf(ob);
            let op = kz_lzma_encode(src.as_ptr(), n, out.as_mut_ptr(), ob);
            if op == 0 { return None; }
            out.set_len(op);
            return Some(out);
        }
    }
    if codec & CODEC_BZIP3 != 0 {
        // SAFETY: out is bz3_bound(n), which is what bz3_compress needs
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
        // SAFETY: ol is the buffer's true length and is updated in place
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
        // SAFETY: the bound comes from rans itself
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

/// The matching decompressor; `raw` is the decoded length the file claims
pub fn paydecomp(src: &[u8], codec: i32, raw: usize) -> Option<Vec<u8>> {
    let mut out = raw_buf(crate::ckalloc(raw));
    paydecomp_into(src, codec, &mut out, raw).then_some(out)
}

/// paydecomp onto the end of `out`, as the C decodes both frames into one blkbuf; false unless it comes back to `raw` bytes
pub fn paydecomp_into(src: &[u8], codec: i32, out: &mut Vec<u8>, raw: usize) -> bool {
    crate::ckalloc(raw);
    out.reserve(raw);
    let n0 = out.len();
    // SAFETY: out has room for raw bytes past n0; each decoder is told so and checked to have written all of them before the length moves
    unsafe {
        let dst = out.as_mut_ptr().add(n0);
        let ok = if codec & CODEC_ZPAQ != 0 {
            kz_zpaq_decode(src.as_ptr(), src.len(), dst, raw) == 0
        } else if codec & CODEC_LZMA != 0 {
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

/// The CRAM name tokeniser over one chunk, in place as the C's tokc: tok3 writes NUL on every separator, put back here
pub fn tok_encode(names: &mut [u8]) -> Option<Vec<u8>> {
    if names.is_empty() || names.len() > c_int::MAX as usize { return None; }
    // SAFETY: out_len is written by the callee; the returned buffer is malloc'd by C and freed here after copying
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
    // SAFETY: as above; the callee reports the decoded length
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
    fn fdatasync(fd: i32) -> i32;
    fn raise(sig: c_int) -> c_int;
    fn sigemptyset(set: *mut SigSet) -> c_int;
    fn sigaddset(set: *mut SigSet, sig: c_int) -> c_int;
    fn pthread_sigmask(how: c_int, set: *const SigSet, old: *mut SigSet) -> c_int;
    fn nanosleep(req: *const TimeSpec, rem: *mut TimeSpec) -> c_int;
}
#[repr(C)] struct TimeSpec { s: i64, ns: i64 }
/// glibc's sigset_t: 1024 bits
#[repr(C)] pub struct SigSet([u64; 16]);
const SIG_BLOCK: c_int = 0;
const SIG_SETMASK: c_int = 2;

/// A staged output finished before it is renamed over anything: synced, then closed, both checked, since NFS, a
/// quota or a failing disk may report a lost write only at fsync or close (File's drop discards a close error). The
/// rename that follows is not followed by a directory fsync: the replacement is atomic, not promised to survive a
/// crash. As the C's fdfinish; true on success, fd closed either way
pub fn finish_fd(fd: i32) -> bool {
    // SAFETY: fd is owned by the caller and given up here
    let synced = unsafe { fdatasync(fd) } == 0;
    (unsafe { close(fd) } == 0) && synced
}
/// `finish_fd` on a File
pub fn finish(f: std::fs::File) -> bool { use std::os::fd::IntoRawFd; finish_fd(f.into_raw_fd()) }

/// Set through --out-dir's renames: a handler then only records the signal, raised once the renames are done
static INCOMMIT: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
static PENDSIG: std::sync::atomic::AtomicI32 = std::sync::atomic::AtomicI32::new(0);
/// The commit window opens: SIGINT, SIGTERM and SIGHUP blocked in this thread, deferred by a handler on any other
/// (a rayon worker has them unblocked). Returns the mask to put back, as the C's stage_swap
pub fn commit_begin() -> SigSet {
    let (mut hs, mut was) = (SigSet([0; 16]), SigSet([0; 16]));
    // SAFETY: both sets are owned, the size glibc's sigset_t is
    unsafe {
        sigemptyset(&mut hs);
        for s in [1, 2, 15] { sigaddset(&mut hs, s); }
        pthread_sigmask(SIG_BLOCK, &hs, &mut was);
    }
    INCOMMIT.store(true, std::sync::atomic::Ordering::SeqCst);
    was
}
/// The window closes. A deferred signal is raised, pending until the mask is put back, then the handler runs with
/// nothing staged left; unless `quiet`, when the caller refuses and exits anyway
pub fn commit_end(was: SigSet, quiet: bool) {
    use std::sync::atomic::Ordering::SeqCst;
    INCOMMIT.store(false, SeqCst);
    if quiet { return; }
    let p = PENDSIG.load(SeqCst);
    // SAFETY: raise and pthread_sigmask on a valid signal and the saved mask
    unsafe { if p != 0 { raise(p); } pthread_sigmask(SIG_SETMASK, &was, std::ptr::null_mut()); }
}

/// -o's temporary name, where it is renamed to and the name given, as the C's OUTTMP, OUTDST and OUTFIN; the signal handler reads the leaked C string, never the mutex
static OUT: std::sync::Mutex<Option<(std::ffi::CString, std::path::PathBuf, String)>> = std::sync::Mutex::new(None);
static OUTSIG: std::sync::atomic::AtomicPtr<c_char> = std::sync::atomic::AtomicPtr::new(std::ptr::null_mut());
/// --out-dir's staged BEDs for the handler: a leaked array of leaked C strings
static STAGE: std::sync::atomic::AtomicPtr<*mut c_char> = std::sync::atomic::AtomicPtr::new(std::ptr::null_mut());
static STAGEN: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
/// KZ_TEST_SIG_PAUSE_MS=N, for tests only: the handler sleeps N ms right after it decides, to exit or to defer, so
/// a test can move the main thread on meanwhile (KZ_TEST_SWAP_PAUSE_MS, db.rs, holds it before and inside the commit
/// window). Read once by sig_test_init at startup, before any handler is installed (getenv is not async-signal-safe);
/// unset or 0, nothing sleeps. As the C's SIGPAUSE_MS
static SIGPAUSE_MS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub fn sig_test_init() {
    let ms = std::env::var("KZ_TEST_SIG_PAUSE_MS").ok().and_then(|v| v.parse().ok()).unwrap_or(0);
    SIGPAUSE_MS.store(ms, std::sync::atomic::Ordering::SeqCst);
}
/// nanosleep, async-signal-safe, through EINTR
fn sig_pause() {
    let ms = SIGPAUSE_MS.load(std::sync::atomic::Ordering::SeqCst);
    if ms == 0 { return; }
    let mut t = TimeSpec { s: (ms / 1000) as i64, ns: (ms % 1000) as i64 * 1_000_000 };
    // SAFETY: both pointers are to the one owned timespec, which nanosleep may rewrite with what is left
    while unsafe { nanosleep(&t, &mut t) } != 0 {}
}
extern "C" fn outsig(sig: c_int) {
    use std::sync::atomic::Ordering::SeqCst;
    if INCOMMIT.load(SeqCst) { PENDSIG.store(sig, SeqCst); sig_pause(); return; }
    sig_pause();
    let (p, a, n) = (OUTSIG.load(SeqCst), STAGE.load(SeqCst), STAGEN.load(SeqCst));
    // SAFETY: unlink and _exit are async-signal-safe; every pointer is leaked
    unsafe {
        if !p.is_null() { unlink(p); }
        if !a.is_null() { for i in 0..n { unlink(*a.add(i)); } }
        _exit(1)
    }
}

/// SIGINT, SIGTERM and SIGHUP remove -o's and --out-dir's temporaries, as the C.
fn catch_signals() {
    // SAFETY: outsig only unlinks and exits
    unsafe { for s in [1, 2, 15] { signal(s, outsig as extern "C" fn(c_int) as usize); } }
}

/// The staged BEDs a signal removes; empty once they are swapped in
pub fn sig_stage(paths: &[std::path::PathBuf]) {
    use std::{os::unix::ffi::OsStrExt, sync::atomic::Ordering::SeqCst};
    let n = STAGEN.swap(0, SeqCst);
    let a = STAGE.swap(std::ptr::null_mut(), SeqCst);
    // SAFETY: a and its strings came from the leaks below, n long
    if !a.is_null() { unsafe { for c in Box::from_raw(std::ptr::slice_from_raw_parts_mut(a, n)) { drop(std::ffi::CString::from_raw(c)); } } }
    if paths.is_empty() { return; }
    let v: Vec<*mut c_char> = paths.iter().filter_map(|p| std::ffi::CString::new(p.as_os_str().as_bytes()).ok())
        .map(|c| c.into_raw()).collect();
    let n = v.len();
    STAGE.store(Box::leak(v.into_boxed_slice()).as_mut_ptr(), SeqCst);
    STAGEN.store(n, SeqCst);
    catch_signals();
}

/// Frees the handler's copy of -o's temporary name
fn sig_out_clear() {
    let p = OUTSIG.swap(std::ptr::null_mut(), std::sync::atomic::Ordering::SeqCst);
    // SAFETY: p came from CString::into_raw in redirect
    if !p.is_null() { unsafe { drop(std::ffi::CString::from_raw(p)); } }
}

/// A panic is a bug, yet leaves no temporaries behind either
pub fn panic_clean() {
    let hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |i| {
        if i.payload().is::<crate::Deferred>() { return; }   // a refusal on its way to the walk
        crate::db::drop_tmp(); drop_out(); hook(i)
    }));
}

/// Removes -o's temporary file on any failure; try_lock, since a panic hook may run with the lock held
pub fn drop_out() {
    sig_out_clear();
    if let Ok(mut o) = OUT.try_lock() { if let Some((t, ..)) = o.take() { unsafe { unlink(t.as_ptr()); } } }
}

/// On success: -o's temporary file renamed into place
pub fn commit_out() {
    use std::io::Write;
    let Some((t, dst, fin)) = OUT.lock().ok().and_then(|mut o| o.take()) else { return };
    sig_out_clear();
    let tp = std::path::PathBuf::from(std::ffi::OsStr::from_bytes(t.as_bytes()));
    // NFS reports some errors only at close
    if std::io::stdout().flush().is_err() || !finish_fd(1) {
        let _ = std::fs::remove_file(&tp); crate::bail("write failed");
    }
    if std::fs::rename(&tp, &dst).is_err() { let _ = std::fs::remove_file(&tp); crate::bail(&format!("cannot rename into {fin}")) }
}

/// A new file PREFIX.XXXXXX, as mkstemp: created exclusively, so never through a symlink, and a taken name retried
pub fn mkstemp(prefix: &std::path::Path, mode: u32) -> std::io::Result<(std::fs::File, std::path::PathBuf)> {
    use std::{hash::{BuildHasher, Hasher}, os::unix::fs::OpenOptionsExt};
    const A: &[u8; 62] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
    let mut n = prefix.as_os_str().as_bytes().to_vec(); n.extend_from_slice(b".XXXXXX");
    let at = n.len() - 6;
    for i in 0..1000u64 {
        // RandomState: random keys a process, moved on at every call
        let mut h = std::collections::hash_map::RandomState::new().build_hasher();
        h.write_u64(i); h.write_u128(std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0, |d| d.as_nanos()));
        let mut r = h.finish();
        for c in &mut n[at..] { *c = A[(r % 62) as usize]; r /= 62; }
        let p = std::path::PathBuf::from(std::ffi::OsStr::from_bytes(&n));
        match std::fs::OpenOptions::new().write(true).create_new(true).mode(mode).open(&p) {
            Ok(f) => return Ok((f, p)),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {}
            Err(e) => return Err(e),
        }
    }
    Err(std::io::ErrorKind::AlreadyExists.into())
}

/// -i and -o in place of stdin and stdout. -o to a regular file, or to none yet, is staged in a new random name
/// beside it, renamed in by `commit_out`, so a failure to stage is fatal and never truncates the old file. A symlink
/// is resolved and its target staged and replaced, the link kept; a dangling one is refused. As the C's redirect
pub fn redirect(path: &str, fd: i32) {
    use std::os::{fd::AsRawFd, unix::fs::PermissionsExt};
    fn refuse(p: &str) -> ! { crate::bail(&format!("cannot open {p}")) }
    let mut old = std::fs::symlink_metadata(path).ok();
    let mut dst = std::path::PathBuf::from(path);
    if fd == 1 && old.as_ref().is_some_and(|m| m.file_type().is_symlink()) {
        dst = std::fs::canonicalize(path).unwrap_or_else(|_| refuse(path));
        old = Some(std::fs::metadata(&dst).unwrap_or_else(|_| refuse(path)));
    }
    let f = if fd == 1 && old.as_ref().is_none_or(|m| m.file_type().is_file()) {
        let mut pre = dst.clone().into_os_string(); pre.push(".tmp");
        let (h, t) = mkstemp(pre.as_ref(), 0o644)
            .unwrap_or_else(|_| crate::bail(&format!("cannot create a temporary file beside {path}")));
        if let Some(m) = &old { let _ = h.set_permissions(std::fs::Permissions::from_mode(m.permissions().mode() & 0o7777)); }
        let c = std::ffi::CString::new(t.into_os_string().into_vec()).unwrap_or_else(|_| refuse(path));
        OUTSIG.store(c.clone().into_raw(), std::sync::atomic::Ordering::SeqCst);
        *OUT.lock().unwrap() = Some((c, dst, path.to_string()));
        catch_signals();
        Ok(h)
    } else if fd == 0 { std::fs::File::open(path) } else { std::fs::OpenOptions::new().write(true).open(&dst) };
    if !f.is_ok_and(|f| unsafe { dup2(f.as_raw_fd(), fd) } >= 0) { refuse(path) }
}

// glibc's byte search: fewer instructions than the memchr crate on short fields, the code the C calls
unsafe extern "C" {
    #[link_name = "memchr"] fn c_memchr(s: *const u8, c: i32, n: usize) -> *const u8;
    #[link_name = "memrchr"] fn c_memrchr(s: *const u8, c: i32, n: usize) -> *const u8;
    #[link_name = "strchrnul"] fn c_strchrnul(s: *const u8, c: i32) -> *const u8;
}
fn at(s: &[u8], p: *const u8) -> Option<usize> { if p.is_null() { None } else { Some(p as usize - s.as_ptr() as usize) } }
/// The first `b` in `s`
pub fn find(b: u8, s: &[u8]) -> Option<usize> { at(s, unsafe { c_memchr(s.as_ptr(), b as i32, s.len()) }) }
/// The last `b` in `s`
pub fn rfind(b: u8, s: &[u8]) -> Option<usize> { at(s, unsafe { c_memrchr(s.as_ptr(), b as i32, s.len()) }) }
/// The first newline or NUL in `s`, which must end in a newline
pub fn line_end(s: &[u8]) -> usize {
    assert_eq!(s.last(), Some(&b'\n'));
    unsafe { c_strchrnul(s.as_ptr(), b'\n' as i32) as usize - s.as_ptr() as usize }
}

// BGZF blocks state their size: libdeflate's whole-buffer inflate and PCLMUL CRC32, as the C's
#[repr(C)] struct LdDec { _p: [u8; 0] }
unsafe extern "C" {
    fn libdeflate_alloc_decompressor() -> *mut LdDec;
    fn libdeflate_free_decompressor(d: *mut LdDec);
    fn libdeflate_deflate_decompress(d: *mut LdDec, i: *const u8, il: usize, o: *mut u8, ol: usize, got: *mut usize) -> i32;
    fn libdeflate_crc32(crc: u32, b: *const u8, n: usize) -> u32;
}
pub struct Inflater(*mut LdDec);
unsafe impl Send for Inflater {}
impl Default for Inflater { fn default() -> Self { Inflater::new() } }
impl Inflater {
    pub fn new() -> Inflater { let d = unsafe { libdeflate_alloc_decompressor() }; assert!(!d.is_null(), "out of memory"); Inflater(d) }
    /// Inflates raw deflate `i` into all of `o`: true when it fills `o` exactly
    pub fn inflate(&mut self, i: &[u8], o: &mut [u8]) -> bool {
        let mut got = 0usize;
        (unsafe { libdeflate_deflate_decompress(self.0, i.as_ptr(), i.len(), o.as_mut_ptr(), o.len(), &mut got) }) == 0 && got == o.len()
    }
}
impl Drop for Inflater { fn drop(&mut self) { unsafe { libdeflate_free_decompressor(self.0) } } }
pub fn crc32(b: &[u8]) -> u32 { unsafe { libdeflate_crc32(0, b.as_ptr(), b.len()) } }

unsafe extern "C" { fn lseek(fd: i32, off: i64, whence: i32) -> i64; }
/// Where stdin's file offset is; None for a pipe
pub fn stdin_at() -> Option<u64> { u64::try_from(unsafe { lseek(0, 0, 1) }).ok() }
