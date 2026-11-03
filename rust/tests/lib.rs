//! The library writes the same bytes as the CLI and refuses the same contract violations: the Rust half of tests/libtest.c
use karyozip::{ffi, writer::{self, Options, Writer}};
use std::io::Write;

/// One BED through `Writer::row`, as an annotator would call it

fn build(bed: &[u8], opt: Options) -> Vec<u8> {
    let dir = std::env::temp_dir();
    let mut w = Writer::new(opt, &dir).expect("temp file");
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

/// Every contract violation is refused, not written out
#[test]
fn refuses_what_the_header_promises() {
    let big = vec![b'x'; 4096];
    let bad: Vec<Box<dyn Fn(&mut Writer) -> Result<(), String>>> = vec![
        Box::new(|w| w.row(b"r1", 10, b"A", 5)),                        // not at 0
        Box::new(|w| { w.row(b"r1", 10, b"A", 0)?; w.row(b"r1", 10, b"A", 50) }), // gap
        Box::new(|w| { w.row(b"r1", 10, b"A", 0)?; w.row(b"r1", 10, b"A", 5) }),  // overlap
        Box::new(|w| { w.row(b"r1", 10, b"A", 0)?; w.row(b"r2", 10, b"A", 0)?;
                       w.row(b"r1", 10, b"A", 10) }),                   // a sequence returns
        Box::new(|w| w.row(b"r1", 0, b"A", 0)),                         // zero-length run
        Box::new(move |w| w.row(&big, 10, b"A", 0)),                    // name too long
        Box::new(|w| w.row(b"", 10, b"A", 0)),                          // empty name
        Box::new(|w| w.row(b"r1", 10, b"A\nB", 0)),                     // newline label
        Box::new(|w| w.row(b"r\0x", 10, b"A", 0)),                      // NUL in the name
        Box::new(|w| w.row(b"a\tb", 10, b"A", 0)),                      // tab in the name
        Box::new(|w| { w.row(b"r1", u64::MAX, b"A", 0)?;
                       w.row(b"r1", 10, b"A", u64::MAX) }),             // coordinate overflow
    ];
    let dir = std::env::temp_dir();
    for (i, f) in bad.iter().enumerate() {
        let mut w = Writer::new(Options::default(), &dir).expect("temp file");
        assert!(f(&mut w).is_err(), "case {i} was accepted");
    }
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
