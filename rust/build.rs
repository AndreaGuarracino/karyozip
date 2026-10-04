// The codecs are the C ones, from ../c: htscodecs because no Rust port
// exists, libdeflate so both inflate BGZF with one code, and bzip3 because our copy carries two fixes upstream still
// lacks, one of which silently loses half the data on an aligned input.
// See c/bzip3/PROVENANCE.
fn main() {
    let c = std::path::Path::new("../c");
    println!("cargo:rerun-if-changed=../c");
    // The two vendored versions come from the vendored trees, never a copy.
    println!("cargo:rustc-env=HTSCODECS_VERSION={}", grab(&c.join("htscodecs/version.h"), "HTSCODECS_VERSION_TEXT"));
    println!("cargo:rustc-env=LIBDEFLATE_VERSION={}", grab(&c.join("libdeflate/libdeflate.h"), "LIBDEFLATE_VERSION_STRING"));
    let bz3 = grab(&c.join("Makefile"), "BZIP3_VERSION ?=");
    // The usage text is the C's own, so the two cannot drift apart.
    let out = std::path::PathBuf::from(std::env::var("OUT_DIR").unwrap());
    std::fs::write(out.join("usage.txt"), usage(&c.join("karyozip.c"))).unwrap();
    cc::Build::new()
        .files([
            c.join("htscodecs/tokenise_name3.c"),
            c.join("htscodecs/rANS_static4x16pr.c"),
            c.join("htscodecs/rANS_static32x16pr.c"),
            c.join("htscodecs/arith_dynamic.c"),
            c.join("htscodecs/rle.c"),
            c.join("htscodecs/pack.c"),
            c.join("htscodecs/utils.c"),
            c.join("bzip3/libbz3.c"),
            c.join("libdeflate/lib/deflate_decompress.c"),
            c.join("libdeflate/lib/crc32.c"),
            c.join("libdeflate/lib/utils.c"),
            c.join("libdeflate/lib/x86/cpu_features.c"),
            c.join("libdeflate/lib/arm/cpu_features.c"),
        ])
        .include(c)
        .define("VERSION", format!("{bz3:?}").as_str())
        .opt_level(3)
        .flag("-fopenmp")                        // libbz3's blocks, side by side
        .flag("-fno-semantic-interposition")     // cc builds -fPIC: let calls inline, as the C's -fPIE
        .warnings(false)
        .compile("kzcodecs");
    println!("cargo:rustc-link-lib=gomp");
}

/// The quoted or trailing word on the one line naming `key`.
fn grab(path: &std::path::Path, key: &str) -> String {
    let t = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let line = t.lines().find(|l| l.contains(key))
        .unwrap_or_else(|| panic!("{}: no {key}", path.display()));
    match line.split('"').nth(1) {
        Some(v) => v.to_string(),
        None => line.split_whitespace().last().unwrap().to_string(),
    }
}

/// The C's USAGE string literal, unescaped back into plain text.
fn usage(path: &std::path::Path) -> String {
    let t = std::fs::read_to_string(path).unwrap();
    let i = t.find("static const char *const USAGE =").expect("no USAGE in karyozip.c");
    let j = t[i..].find(";\n").expect("unterminated USAGE") + i;
    let mut s = String::new();
    for part in t[i..j].split('"').skip(1).step_by(2) {
        s.push_str(&part.replace("\\n", "\n"));
    }
    assert!(s.starts_with("usage:\n"), "USAGE did not parse");
    s
}
