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

## Tests

`make test` runs `tests/run.sh` on the C (with `c/build/karyozip-san`, the ASan/UBSan build, for the leak and fuzz stages), `tests/run.sh` on the Rust, and `tests/diff.sh`. `make test-rust-debug` runs `cargo test` in a debug build, with overflow checks on. The suites are Linux/GNU only and need python3, zstd, gzip, GNU coreutils (`realpath -m`, `stat -c`, `sha256sum`, `timeout`), GNU grep (`-P`) and awk; each script names what is missing and stops.
The Rust half has no sanitizer binary, so its leak stage runs the Rust under **valgrind**. Without valgrind that stage prints `leaks: skipped: valgrind not installed` and the rest runs; set `REQUIRE_VALGRIND=1` to make it a failure, as CI should.
Every check prints one line, failures start with `FAIL: `, and each script ends with `N passed, M failed, K skipped`. Environment knobs: `CASE=REGEX` (diff.sh and fuzz.sh run only the matching sections; run.sh runs everything but reports only matching checks), `TIMEOUT_SCALE` (multiplies fixed timeouts, default 1, 4 under an ASan binary) and `FUZZ_SEED` (fuzz.sh prints the seed it used; set it to repeat a run's mutations), and three for the CPU budget:
- `TEST_JOBS` is the **total** budget, shared and never multiplied; default `SLURM_CPUS_PER_TASK`, else `nproc`. `make test` runs min(3, `TEST_JOBS`) suites at once and gives each floor(`TEST_JOBS` / suites), at least 1, as its own `TEST_JOBS`; `TEST_JOBS=1` runs the three one after another. Inside a suite it bounds the cases run at once (run.sh's concurrent readers, fuzz.sh's `xargs -P`) and `cargo test`'s `-j` and `--test-threads`.
- `KZ_TEST_THREADS` is the `-t` every karyozip run gets, default the suite's `TEST_JOBS`, so no binary falls back to every allocated CPU inside a parallel suite; N runs at once split it (`KZ_TEST_THREADS / N` each), the Rust's own tests get `RAYON_NUM_THREADS` = `KZ_TEST_THREADS / --test-threads`, and the C library's test, which has no `-t`, runs under `taskset` on `KZ_TEST_THREADS` CPUs. The bytes never depend on it. `-t` bounds the codec threads, not every thread: compress also runs a fixed input reader and a names thread beside them, so a `-t 1` compress can have 2 or 3 runnable threads. The stress cases keep the thread counts they test.
- `TEST_LIGHT=1` is the serial, lightweight mode: it skips the stress cases and prints one `skipped (TEST_LIGHT): CASE` line for each, counted as skipped, not passed. The stress cases are the multi-thread equivalence checks (smoothing at 2, 4 and 8 threads against 1; `--db --aln` and `--read-blocks 7 --smooth` at `-t 4`; decode ahead at `-t 4` against `-t 1`, with its 120000-read input), `-t 8` on 200000 one-row sequences, `-t 1024` accepted, `OMP_NUM_THREADS` ignored at the default thread count and at `--threads 2`, the concurrent readers (up to 8), the full body temp file at `-t 2` (run.sh and diff.sh), and the BGZF 4 GiB and 300000-block inputs under a 1.5 GB address limit. They use the thread counts they test, not the budget, so run them in a bounded Slurm job, not on a laptop: `TEST_JOBS=1 TEST_LIGHT=1 make test` is the run for a small machine.
`ci/github-actions.yml` is a ready GitHub Actions workflow; copy it to `.github/workflows/` to enable it.

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
