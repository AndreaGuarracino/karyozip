# karyozip, in Rust

A second implementation of the same format, to find out whether Rust is the
better host for this tool. It must produce the bytes `c/` produces;
`tests/diff.sh` is the gate.

## Where the codecs come from

None of them is Rust. There is no Rust port of htscodecs, and the `bzip3`
crate lacks the two fixes `c/bzip3/PROVENANCE` describes, one of which
silently loses half the data on an input that is an exact multiple of the
block size. So `build.rs` compiles the same C sources `c/` compiles, and
`src/ffi.rs` wraps them.

zstd comes from `zstd-sys`, pinned to `=2.0.16`, which bundles exactly the
1.5.7 the C links. zstd does not promise identical output across its own
releases, so that pin is what makes byte-identity testable.

## unsafe

`src/ffi.rs` holds the codec calls: every raw pointer and length check for them lives there, and callers see slices and `Option<Vec<u8>>`.
Outside it, four hot loops write unchecked, each behind an assertion and each a measured win: the digit writer in `kst.rs`, the row writers in `reader.rs` and `smooth.rs`, and `smooth_one`.
