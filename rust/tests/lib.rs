//! The library must write the same bytes as the CLI, and refuse the same
//! contract violations. This is the Rust half of c/tests/libtest.c.
use karyozip::writer::{Options, Writer};
use std::io::Write;

/// One BED through `Writer::row`, which is what an annotator would call
/// instead of writing a BED at all.
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

fn cli(args: &[&str], bed: &[u8]) -> Vec<u8> {
    let exe = env!("CARGO_BIN_EXE_karyozip");
    let mut c = std::process::Command::new(exe)
        .args(args).arg("encode")
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
        (vec!["--level", "3"], Options { level: 3, levelset: true, ..Options::default() }),
        (vec!["--names", "keep"], Options { keep_names: true, ..Options::default() }),
    ] {
        assert_eq!(build(&bed, opt), cli(&args, &bed), "differs for {args:?}");
    }
}

/// Every contract violation is refused, not written out.
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
