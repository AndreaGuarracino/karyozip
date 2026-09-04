// The codecs are the C ones, from ../c: htscodecs because no Rust port
// exists, and bzip3 because our copy carries two fixes upstream still
// lacks, one of which silently loses half the data on an aligned input.
// See c/bzip3/PROVENANCE.
fn main() {
    let c = std::path::Path::new("../c");
    println!("cargo:rerun-if-changed=../c");
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
        ])
        .include(c)
        .define("VERSION", "\"1.5.3\"")
        .opt_level(3)
        .warnings(false)
        .compile("kzcodecs");
}
