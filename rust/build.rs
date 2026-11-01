// The codecs are the C ones from ../c: no Rust htscodecs, one BGZF inflate for both, and bzip3 with the fixes in c/bzip3/PROVENANCE
fn main() {
    let c = std::path::Path::new("../c");
    println!("cargo:rerun-if-changed=../c");
    // The two vendored versions come from the vendored trees, never a copy
    println!("cargo:rustc-env=HTSCODECS_VERSION={}", grab(&c.join("htscodecs/version.h"), "HTSCODECS_VERSION_TEXT"));
    println!("cargo:rustc-env=LIBDEFLATE_VERSION={}", grab(&c.join("libdeflate/libdeflate.h"), "LIBDEFLATE_VERSION_STRING"));
    let bz3 = grab(&c.join("Makefile"), "BZIP3_VERSION ?=");
    // The usage text is the C's own, so the two cannot drift apart
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
    // liblzma 5.8.4 on its own: its config.h, common.h and version.h share names with htscodecs' and bzip3's; the C Makefile's file list
    let mk = std::fs::read_to_string(c.join("Makefile")).unwrap();
    let key = "$(addprefix xz/src/liblzma/,";
    let i = mk.find(key).expect("no liblzma list in c/Makefile") + key.len();
    let mut xz = cc::Build::new();
    for f in mk[i..i + mk[i..].find(')').unwrap()].split_whitespace().filter(|w| *w != "\\") {
        xz.file(c.join("xz/src/liblzma").join(f));
    }
    xz.file(c.join("xz/kzlzma.c"));
    xz.define("HAVE_CONFIG_H", None).include(c.join("xz")).include(c.join("xz/src/common"));
    for d in ["api", "common", "check", "lz", "lzma", "rangecoder", "delta", "simple"] {
        xz.include(c.join("xz/src/liblzma").join(d));
    }
    xz.opt_level(3).flag("-std=gnu11").flag("-fno-semantic-interposition").warnings(false).compile("kzxz");
    // libzpaq 7.15 without OpenMP, since divsufsort spawns its own threads; JIT on x86-64 only, as the C's

    let mut zp = cc::Build::new();
    zp.cpp(true).file(c.join("zpaq/libzpaq.cpp")).define("unix", None).opt_level(3).warnings(false);
    if std::env::var("CARGO_CFG_TARGET_ARCH").as_deref() != Ok("x86_64") { zp.define("NOJIT", None); }
    zp.compile("kzlibzpaq");
    let mut zw = cc::Build::new();
    zw.cpp(true).file(c.join("zpaq/kzzpaq.cpp")).include(c.join("zpaq")).define("unix", None).opt_level(3).flag("-fopenmp").warnings(false);
    if std::env::var("CARGO_CFG_TARGET_ARCH").as_deref() != Ok("x86_64") { zw.define("NOJIT", None); }
    zw.compile("kzzpaq");
    println!("cargo:rustc-link-lib=gomp");
}

/// The quoted or trailing word on the one line naming `key`
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

/// The C's USAGE string literal, unescaped back into plain text
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
