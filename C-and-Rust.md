# karyozip in C and in Rust

The [README](README.md) describes the Rust binary. The C in `c/` is the reference implementation: both write the same bytes.

## Build the C

```
make            # builds c/build/karyozip
make rust       # builds rust/target/release/karyozip, the binary the tests use
```

The C needs zstd, zlib, an OpenMP compiler and a C++ compiler for libzpaq. liblzma is vendored in `c/xz`, libzpaq in `c/zpaq`.
`CFLAGS`, `CXXFLAGS`, `CPPFLAGS` and `LDFLAGS` from the environment, as a conda env sets them, are honoured; the flags the build needs are added to them. The build stops unless `zstd.h` is zstd 1.5.7, the one the Rust bundles (`ZSTD_CFLAGS` and `ZSTD_LIBS` choose another zstd).

## C and Rust

`c/` and `rust/` write the same bytes, print the same `--help` and refuse the same inputs with the same words; the usage text is `c/usage.txt`, which the C embeds and the Rust includes, and the liblzma file list is `c/xz/sources.txt`, which `c/Makefile` and `rust/build.rs` both read.
`make test` runs the shared suite in `tests/` on both, fuzzes both, and compares them compress by compress. It takes a few minutes, depending on the machine.
The C half fuzzes an ASan/UBSan build. The Rust half checks for leaks under valgrind, not a sanitizer, so `make test` needs valgrind installed.
The C finds zstd, zlib and libgomp at the path in its rpath; the Rust bundles zstd. Both compile the same liblzma from `c/xz`.
The Rust compiles the C codecs too, through `rust/src/ffi.rs`: there is no Rust htscodecs, and the `bzip3` crate lacks the local changes in `c/bzip3/PROVENANCE` (`bz3_compress_cuts` among them). `zstd-sys` is pinned to `=2.1.2`, which bundles the zstd 1.5.7 the C links, because zstd does not promise the same bytes across releases.
Outside `ffi.rs`, four hot loops write unchecked, each behind an assertion: the digit writer in `kst.rs`, the row writers in `reader.rs` and `smooth.rs`, and `smooth_one`.

## The C library

`make` also builds `c/build/libkaryozip.a`, six calls in `c/karyozip.h` (`kt_open`, `kt_block`, `kt_level`, `kt_row`, `kt_close`, `kt_version`), for an annotator that writes `.kst` directly.
Rows must arrive grouped by sequence, each sequence once, starting at 0, gapless; anything else is refused with the CLI's message.
`kt_block` sets the rows per block, 0 or by default one block (no random access; the CLI's `-e low` uses 262144); `kt_level` the zstd level, default 17. Both set process-wide state: call them after `kt_open` and before the first `kt_row`.
The library is not thread-safe: one writer at a time per process, from one thread; `kt_open` resets the encoder, so a second writer writes what a fresh process would.
Errors are not returned: each prints `karyozip: ...` to stderr and calls `exit(1)`, as the CLI does.

```c
kt_writer *w = kt_open(f);
kt_row(w, name, nlen, run_len, label, llen, start);   /* once per run */
kt_close(w);
```

## Vendored code

htscodecs, bzip3, libdeflate, liblzma and libzpaq under `c/` are pinned; each `PROVENANCE` names the release and every local patch.
The vendored bzip3 also takes each block's size from the caller (`bz3_compress_cuts`), so blocks end at read starts; its frame already recorded each block's size, so the decoder is unchanged.
The vendored bzip3 runs its independent blocks in parallel, writing the same bytes; `c/bzip3/lastblock_check.c` guards the last-block fix that bzip3 1.5.4 took from this tree.
