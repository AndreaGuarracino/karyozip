//! The library writes the same bytes as the CLI and refuses the same contract violations: the Rust half of tests/libtest.c
use karyozip::{ffi, writer::{self, Options, Writer}, Error, ErrorKind, KzResult};
use std::io::Write;

/// One BED through `Writer::row`, as an annotator would call it

fn build(bed: &[u8], opt: Options) -> Vec<u8> {
    let dir = std::env::temp_dir();
    build_in(bed, opt, &dir)
}
fn build_in(bed: &[u8], opt: Options, dir: &std::path::Path) -> Vec<u8> {
    let mut w = Writer::new(opt, dir).expect("temp file");
    for line in bed.split(|&c| c == b'\n') {
        if line.is_empty() { continue; }
        let f: Vec<&[u8]> = line.splitn(4, |&c| c == b'\t').collect();
        let st: u64 = std::str::from_utf8(f[1]).unwrap().parse().unwrap();
        let en: u64 = std::str::from_utf8(f[2]).unwrap().parse().unwrap();
        w.row(f[0], en - st, f[3], st).expect("row");
    }
    let mut out = Vec::new();
    w.close(&mut out).expect("close");
    out
}

/// The -t of every CLI child: KZ_TEST_THREADS (the Makefile passes each test's share), else 1; the bytes never depend on it
fn test_threads() -> String {
    match std::env::var("KZ_TEST_THREADS").ok().and_then(|v| v.parse::<u32>().ok()) {
        Some(n) if (1..=1024).contains(&n) => n.to_string(),
        _ => "1".into(),
    }
}

fn cli(args: &[&str], bed: &[u8]) -> Vec<u8> {
    let exe = env!("CARGO_BIN_EXE_karyozip");
    let mut c = std::process::Command::new(exe)
        .arg("-t").arg(test_threads())
        .args(args).arg("compress")
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::null())
        .spawn().expect("spawn");
    c.stdin.take().unwrap().write_all(bed).unwrap();
    c.wait_with_output().unwrap().stdout
}

fn fixture() -> Vec<u8> {
    let mut b = Vec::new();
    for s in 0..40 {
        let mut p = 0u64;
        for k in 0..25 {
            let l = 7 + (s * 13 + k * 5) % 400;
            b.extend_from_slice(format!("r{s}\t{p}\t{}\tL{}\n", p + l, k % 6).as_bytes());
            p += l;
        }
    }
    b
}

#[test]
fn same_bytes_as_the_cli() {
    let bed = fixture();
    for (args, opt) in [
        (vec![], Options::default()),
        (vec!["--block", "4"], Options { block: 4, ..Options::default() }),
        (vec!["--block", "all"], Options { block: u64::MAX, ..Options::default() }),
        (vec!["--names", "drop"], Options { drop_names: true, ..Options::default() }),
    ] {
        assert_eq!(build(&bed, opt), cli(&args, &bed), "differs for {args:?}");
    }
}

/// Every contract violation is refused, not written out, in the CLI's words
#[test]
fn refuses_what_the_header_promises() {
    let big = vec![b'x'; 4096];
    type Case = (Box<dyn Fn(&mut Writer) -> KzResult<()>>, &'static str);
    let bad: Vec<Case> = vec![
        (Box::new(|w| w.row(b"r1", 10, b"A", 5)), "a sequence does not start at 0"),
        (Box::new(|w| { w.row(b"r1", 10, b"A", 0)?; w.row(b"r1", 10, b"A", 50) }), "gap or overlap: start is not the previous end"),
        (Box::new(|w| { w.row(b"r1", 10, b"A", 0)?; w.row(b"r1", 10, b"A", 5) }), "gap or overlap: start is not the previous end"),
        (Box::new(|w| { w.row(b"r1", 10, b"A", 0)?; w.row(b"r2", 10, b"A", 0)?;
                        w.row(b"r1", 10, b"A", 10) }), "a sequence does not start at 0"),   // a sequence returns
        (Box::new(|w| w.row(b"r1", 0, b"A", 0)), "end equals start: a run must cover at least one base"),
        (Box::new(move |w| w.row(&big, 10, b"A", 0)), "name too long"),
        (Box::new(|w| w.row(b"", 10, b"A", 0)), "empty sequence name"),
        (Box::new(|w| w.row(b"r1", 10, b"A\nB", 0)), "a newline inside a name or label"),
        (Box::new(|w| w.row(b"r\0x", 10, b"A", 0)), "a NUL byte inside a row"),
        (Box::new(|w| w.row(b"a\tb", 10, b"A", 0)), "a tab inside a sequence name"),
        (Box::new(|w| { w.row(b"r1", u64::MAX, b"A", 0)?;
                        w.row(b"r1", 10, b"A", u64::MAX) }), "a coordinate does not fit in 64 bits"),
    ];
    let dir = std::env::temp_dir();
    for (i, (f, want)) in bad.iter().enumerate() {
        let mut w = Writer::new(Options::default(), &dir).expect("temp file");
        let e = f(&mut w).expect_err(&format!("case {i} was accepted"));
        assert_eq!((e.message(), e.kind()), (*want, ErrorKind::Input), "case {i}");
    }
    // a repeated name is found at the close
    let mut w = Writer::new(Options::default(), &dir).expect("temp file");
    for n in [&b"r1"[..], b"r2", b"r1"] { w.row(n, 10, b"A", 0).expect("row"); }
    let e = w.close(&mut Vec::new()).err().expect("a repeated name was accepted");
    assert_eq!(e.message(), "a sequence name repeats");
}

/// A directory of its own under cargo's scratch space, so a test can count what is left in it
fn scratch(name: &str) -> std::path::PathBuf {
    let d = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("kz-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).expect("scratch dir");
    d
}
/// karyozip.* temporaries left under `d`
fn leftovers(d: &std::path::Path) -> Vec<String> {
    std::fs::read_dir(d).map(|r| r.filter_map(|e| e.ok()).map(|e| e.file_name().to_string_lossy().into_owned())
                                  .filter(|n| n.starts_with("karyozip")).collect()).unwrap_or_default()
}
fn tmp_msg(d: &std::path::Path) -> String { format!("cannot create a temp file under {} (see --tmp-dir)", d.display()) }

fn feed(w: &mut Writer, bed: &[u8]) {
    for line in bed.split(|&c| c == b'\n').filter(|l| !l.is_empty()) {
        let f: Vec<&[u8]> = line.splitn(4, |&c| c == b'\t').collect();
        let st: u64 = std::str::from_utf8(f[1]).unwrap().parse().unwrap();
        let en: u64 = std::str::from_utf8(f[2]).unwrap().parse().unwrap();
        w.row(f[0], en - st, f[3], st).expect("row");
    }
}

/// A temp directory the writer cannot use is an Err, in the CLI's words, and the process goes on
#[test]
fn a_bad_tmp_dir_is_an_error() {
    let bed = fixture();
    let d = scratch("tmpdir");
    // missing
    let gone = d.join("nosuch");
    let e = Writer::new(Options::default(), &gone).err().expect("a missing tmp dir was accepted");
    assert_eq!((e.message(), e.kind()), (tmp_msg(&gone).as_str(), ErrorKind::Io));
    assert_eq!(e.to_string(), tmp_msg(&gone));
    // moved away between new and close (a rename, since NFS keeps an unlinked open file as .nfs*): the close's body file has nowhere to go
    let sub = d.join("moved"); std::fs::create_dir(&sub).unwrap();
    let mut w = Writer::new(Options { block: 4, ..Options::default() }, &sub).expect("temp file");
    feed(&mut w, &bed);
    std::fs::rename(&sub, d.join("moved.away")).unwrap();
    let e = w.close(&mut Vec::new()).err().expect("close without a tmp dir was accepted");
    assert_eq!((e.message(), e.kind()), (tmp_msg(&sub).as_str(), ErrorKind::Io));
    assert!(leftovers(&d.join("moved.away")).is_empty(), "left: {:?}", leftovers(&d.join("moved.away")));
    // read-only (root writes anyway: then only that it works)
    use std::os::unix::fs::PermissionsExt;
    let ro = d.join("ro"); std::fs::create_dir(&ro).unwrap();
    std::fs::set_permissions(&ro, std::fs::Permissions::from_mode(0o555)).unwrap();
    let root = std::fs::File::create(ro.join("probe")).is_ok();
    if !root {
        let e = Writer::new(Options::default(), &ro).err().expect("a read-only tmp dir was accepted");
        assert_eq!(e.message(), tmp_msg(&ro));
        // writable at new, read-only at close
        std::fs::set_permissions(&ro, std::fs::Permissions::from_mode(0o755)).unwrap();
        let mut w = Writer::new(Options::default(), &ro).expect("temp file");
        feed(&mut w, &bed);
        std::fs::set_permissions(&ro, std::fs::Permissions::from_mode(0o555)).unwrap();
        let e = w.close(&mut Vec::new()).err().expect("close in a read-only dir was accepted");
        assert_eq!(e.message(), tmp_msg(&ro));
        assert!(leftovers(&ro).is_empty(), "left in the tmp dir: {:?}", leftovers(&ro));
    }
    std::fs::set_permissions(&ro, std::fs::Permissions::from_mode(0o755)).unwrap();
    // still alive, and the next writer in the same directory gives the CLI's bytes
    assert_eq!(build_in(&bed, Options::default(), &d), cli(&[], &bed));
    assert!(leftovers(&d).is_empty(), "left in the tmp dir: {:?}", leftovers(&d));
    let _ = std::fs::remove_dir_all(&d);
}

/// An output that fails after `room` bytes, as a full disk does
struct Full { room: usize, got: Vec<u8> }
impl Write for Full {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if self.got.len() + b.len() > self.room {
            return Err(std::io::Error::new(std::io::ErrorKind::StorageFull, "no space left"));
        }
        self.got.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> { Ok(()) }
}

/// A failing output is "write failed", as the CLI says it; no temp file stays, and the process goes on
#[test]
fn a_failing_output_is_an_error() {
    let bed = fixture();
    let d = scratch("full");
    let whole = build_in(&bed, Options { block: 4, ..Options::default() }, &d);
    // nothing, inside the header, inside the tables, inside the body, one byte short
    for room in [0, 8, 100, whole.len() / 2, whole.len() - 1] {
        let mut w = Writer::new(Options { block: 4, ..Options::default() }, &d).expect("temp file");
        feed(&mut w, &bed);
        let mut out = Full { room, got: Vec::new() };
        let e = w.close(&mut out).err().unwrap_or_else(|| panic!("an output full at {room} B was accepted"));
        assert_eq!((e.message(), e.kind()), ("write failed", ErrorKind::Io), "full at {room} B");
        assert!(out.got.len() <= room);
        assert_eq!(out.got[..], whole[..out.got.len()], "what was written is the file's start");
        assert!(leftovers(&d).is_empty(), "left in the tmp dir: {:?}", leftovers(&d));
    }
    let mut out = Full { room: whole.len(), got: Vec::new() };
    let mut w = Writer::new(Options { block: 4, ..Options::default() }, &d).expect("temp file");
    feed(&mut w, &bed);
    w.close(&mut out).expect("an output with room");
    assert_eq!(out.got, whole);
    // a Writer dropped unclosed, rows and all, leaves nothing either
    let mut w = Writer::new(Options::default(), &d).expect("temp file");
    feed(&mut w, &bed);
    drop(w);
    assert!(leftovers(&d).is_empty(), "left in the tmp dir: {:?}", leftovers(&d));
    let _ = std::fs::remove_dir_all(&d);
}

/// Writers one after another in one process, and two at once on two threads, each the CLI's bytes
#[test]
fn writers_in_sequence_and_at_once() {
    let bed = fixture();
    let want = cli(&["--block", "4"], &bed);
    let opt = || Options { block: 4, ..Options::default() };
    for _ in 0..3 { assert_eq!(build(&bed, opt()), want, "a later writer in the same process"); }
    let (a, b) = std::thread::scope(|s| {
        let a = s.spawn(|| build(&bed, opt()));
        let b = s.spawn(|| build(&bed, opt()));
        (a.join().unwrap(), b.join().unwrap())
    });
    assert_eq!((a == want, b == want), (true, true), "two writers at once");
}

/// The error type: a message, a kind, and std's Error
#[test]
fn the_error_type() {
    let e = Error::input("name too long");
    assert_eq!((e.message(), e.kind(), e.to_string()), ("name too long", ErrorKind::Input, "name too long".to_string()));
    let b: Box<dyn std::error::Error + Send + Sync> = Box::new(e.clone());
    assert_eq!(b.to_string(), "name too long");
}

/// One file of tests/fixtures.txt, where "@@ NAME" opens each
fn fixfile(name: &str) -> Vec<u8> {
    let (all, mut out, mut on) = (std::fs::read_to_string(format!("{}/../tests/fixtures.txt", env!("CARGO_MANIFEST_DIR"))).unwrap(), String::new(), false);
    for l in all.lines() { if let Some(n) = l.strip_prefix("@@ ") { on = n == name } else if on { out += l; out.push('\n') } }
    out.into_bytes()
}

/// Every codec candidate, and the name tokeniser, round-trips its input
#[test]
fn every_codec_round_trips() {
    for f in ["db.rall.bed", "db.chromosome.bed", "ls.bed"] {
        let mut src = fixfile(f);
        let mut n = 0;
        for c in writer::candlist(true) {
            if let Some(v) = ffi::paycomp(&src, c) {
                n += 1;
                assert_eq!(ffi::paydecomp(&v, c, src.len()).as_deref(), Some(&src[..]), "codec 0x{c:x} on {f}");
            }
        }
        assert!(n > 0, "no codec took {f}");
        if let Some(v) = ffi::tok_encode(&mut src) {
            assert_eq!(ffi::tok_decode(&v).as_deref(), Some(&src[..]), "tokeniser on {f}");
        }
    }
}
