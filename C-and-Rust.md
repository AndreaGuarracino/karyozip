# karyozip in C and in Rust

The [README](README.md) describes the Rust binary. The C in `c/` is the reference implementation: both write the same bytes.

## Build the C

```
make            # builds c/build/karyozip
make rust       # builds rust/target/release/karyozip, the binary the tests use
```

The C needs zstd, zlib, an OpenMP compiler and a C++ compiler for libzpaq. liblzma is vendored in `c/xz`, libzpaq in `c/zpaq`.

## C and Rust

`c/` and `rust/` write the same bytes, print the same `--help` and refuse the same inputs with the same words; `rust/build.rs` takes the usage text from `c/karyozip.c`.
`make test` (65 s) runs the shared suite in `tests/`, sanitizers and fuzzing on both and compares them compress by compress.
The C finds zstd, zlib and libgomp at the path in its rpath; the Rust bundles zstd. Both compile the same liblzma from `c/xz`.
The Rust compiles the C codecs too, through `rust/src/ffi.rs`: there is no Rust htscodecs, and the `bzip3` crate lacks the two fixes in `c/bzip3/PROVENANCE`. `zstd-sys` is pinned to `=2.0.16`, the zstd 1.5.7 the C links, because zstd does not promise the same bytes across releases.
Outside `ffi.rs`, four hot loops write unchecked, each behind an assertion: the digit writer in `kst.rs`, the row writers in `reader.rs` and `smooth.rs`, and `smooth_one`.

## The C library

`make` also builds `c/build/libkaryozip.a`, six calls in `c/karyozip.h` (`kt_open`, `kt_block`, `kt_level`, `kt_row`, `kt_close`, `kt_version`), for an annotator that writes `.kst` directly.
Rows must arrive grouped by sequence, each sequence once, starting at 0, gapless; anything else is refused with the CLI's message.

```c
kt_writer *w = kt_open(f);
kt_row(w, name, nlen, run_len, label, llen, start);   /* once per run */
kt_close(w);
```

## Vendored code

htscodecs, bzip3, libdeflate, liblzma and libzpaq under `c/` are pinned; each `PROVENANCE` names the release and every local patch.
The vendored bzip3 also takes each block's size from the caller (`bz3_compress_cuts`), so blocks end at read starts; its frame already recorded each block's size, so the decoder is unchanged.
The vendored bzip3 runs its independent blocks in parallel, writing the same bytes; `c/bzip3/lastblock_check.c` guards the last-block fix that bzip3 1.5.4 took from this tree.
