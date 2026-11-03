//! Every text input goes through karyozip::unzip: plain, BGZF blocks inflated side by side, or any gzip.
#![no_main]
use libfuzzer_sys::fuzz_target;
use std::io::Read;
#[path = "common.rs"] mod common;

fuzz_target!(|data: &[u8]| {
    let d = data.to_vec();
    common::refusals_ok(move || {
        let mut out = Vec::new();
        let _ = karyozip::unzip(std::io::Cursor::new(d)).read_to_end(&mut out);
    });
});
