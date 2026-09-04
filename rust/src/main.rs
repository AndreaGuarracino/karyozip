mod ffi;

use std::io::{Read, Write};

/// Milestone check: every codec called from Rust must give the bytes the C
/// gives. `codecs FILE` prints one line per candidate.
fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() >= 2 && a[1] == "--version" {
        println!("karyozip {} (rust), format KST00001", env!("CARGO_PKG_VERSION"));
        println!("  zstd       {}", ffi::zstd_version());
        println!("  libbz3     {}", ffi::bzip3_version());
        return;
    }
    if a.len() == 3 && a[1] == "codecs" {
        let mut src = Vec::new();
        std::fs::File::open(&a[2]).unwrap().read_to_end(&mut src).unwrap();
        let mut out = std::io::stdout().lock();
        for lvl in [1, 3, 17, 19, 22] {
            let z = ffi::zstd(&src, lvl);
            writeln!(out, "zstd:{lvl}\t{}\t{:016x}", z.len(), fnv(&z)).unwrap();
        }
        // the eight default candidates, then the six bzip3 block sizes
        for c in [0x100, 1, 65, 521, ffi::CODEC_ARITH, ffi::CODEC_ARITH | 1,
                  ffi::CODEC_BZIP3 | 16, ffi::CODEC_BZIP3 | 64] {
            match ffi::paycomp(&src, c) {
                Some(v) => {
                    let back = ffi::paydecomp(&v, c, src.len());
                    let ok = back.as_deref() == Some(&src[..]);
                    writeln!(out, "codec:0x{c:x}\t{}\t{:016x}\troundtrip={}", v.len(), fnv(&v), ok).unwrap();
                }
                None => writeln!(out, "codec:0x{c:x}\trefused").unwrap(),
            }
        }
        match ffi::tok_encode(&src) {
            Some(v) => {
                let back = ffi::tok_decode(&v);
                writeln!(out, "tok\t{}\t{:016x}\troundtrip={}", v.len(), fnv(&v),
                         back.as_deref() == Some(&src[..])).unwrap();
            }
            None => writeln!(out, "tok\trefused").unwrap(),
        }
        return;
    }
    eprintln!("usage: karyozip --version | codecs FILE");
    std::process::exit(2);
}

fn fnv(b: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &x in b { h ^= x as u64; h = h.wrapping_mul(0x100_0000_01b3); }
    h
}
