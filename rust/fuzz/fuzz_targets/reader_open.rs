//! A .kst from the fuzzer's bytes: open it, then decode every row. Reaches the header checks,
//! block_dec, load_packed and the chunk tables.
#![no_main]
use libfuzzer_sys::fuzz_target;
use std::io::Write;
#[path = "common.rs"] mod common;

fuzz_target!(|data: &[u8]| {
    // the reader takes a path: one temp file a process, rewritten each run
    let path = std::env::temp_dir().join(format!("karyozip-fuzz-{}.kst", std::process::id()));
    std::fs::File::create(&path).and_then(|mut f| f.write_all(data)).expect("temp file");
    let p = path.to_str().unwrap().to_owned();
    common::refusals_ok(move || {
        let mut r = karyozip::reader::Reader::open_with(&p, &karyozip::reader::ReadOpts::default());
        r.decode(&mut std::io::sink());
    });
});
