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

`c/` and `rust/` write the same bytes, print the same `--help` and refuse the same inputs with the same words, `compress -d` from per-set `.kst` inputs included (the C walks each input as `dump` does, setting and restoring the reader's per-track globals around every row, as `reffill_r` does for `--ref`; the Rust pulls rows through `Reader::next_row`); the usage text is `c/usage.txt`, which the C embeds and the Rust includes, and the liblzma file list is `c/xz/sources.txt`, which `c/Makefile` and `rust/build.rs` both read.
`make test` runs the shared suite in `tests/` on both, fuzzes both, and compares them compress by compress. It takes a few minutes, depending on the machine.
The C half fuzzes an ASan/UBSan build. The Rust half checks for leaks under valgrind, not a sanitizer, so `make test` needs valgrind installed.
The C finds zstd, zlib and libgomp at the path in its rpath; the Rust bundles zstd. Both compile the same liblzma from `c/xz`.
The Rust compiles the C codecs too, through `rust/src/ffi.rs`: there is no Rust htscodecs, and the `bzip3` crate lacks the local changes in `c/bzip3/PROVENANCE` (`bz3_compress_cuts` among them). `zstd-sys` is pinned to `=2.1.2`, which bundles the zstd 1.5.7 the C links, because zstd does not promise the same bytes across releases.
Outside `ffi.rs`, four hot loops index or write unchecked, each behind a reserve or an assertion: in `kst.rs` the row writer `put_row` and the digit writer `putu64_at`, in `smooth.rs` `smooth_one` and the row writer `emit`. `reader.rs` has no `unsafe`.

## Tests

`make test` runs `tests/run.sh` on the C (with `c/build/karyozip-san`, the ASan/UBSan build, for the leak and fuzz stages), `tests/run.sh` on the Rust, and `tests/diff.sh`. `make test-rust-debug` runs `cargo test` in a debug build, with overflow checks on. The suites are Linux/GNU only and need python3, zstd, gzip, GNU coreutils (`realpath -m`, `stat -c`, `sha256sum`, `timeout`), GNU grep (`-P`) and awk; each script names what is missing and stops.
The Rust half has no sanitizer binary, so its leak stage runs the Rust under **valgrind**. Without valgrind that stage prints `leaks: skipped: valgrind not installed` and the rest runs; set `REQUIRE_VALGRIND=1` to make it a failure, as CI should.
Every check prints one line, failures start with `FAIL: `, and each script ends with `N passed, M failed, K skipped`. Environment knobs: `CASE=REGEX` (diff.sh and fuzz.sh run only the matching sections; run.sh runs everything but reports only matching checks), `TIMEOUT_SCALE` (multiplies fixed timeouts, default 1, 4 under an ASan binary) and `FUZZ_SEED` (fuzz.sh prints the seed it used; set it to repeat a run's mutations), and three for the CPU budget:
- `TEST_JOBS` is the **total** budget, shared and never multiplied; default `SLURM_CPUS_PER_TASK`, else `nproc`. `make test` runs min(3, `TEST_JOBS`) suites at once and gives each floor(`TEST_JOBS` / suites), at least 1, as its own `TEST_JOBS`; `TEST_JOBS=1` runs the three one after another. Inside a suite it bounds the cases run at once (run.sh's concurrent readers, fuzz.sh's `xargs -P`) and `cargo test`'s `-j` and `--test-threads`.
- `KZ_TEST_THREADS` is the `-t` every karyozip run gets, default the suite's `TEST_JOBS`, so no binary falls back to every allocated CPU inside a parallel suite; N runs at once split it (`KZ_TEST_THREADS / N` each), the Rust's own tests get `RAYON_NUM_THREADS` = `KZ_TEST_THREADS / --test-threads` (and pass that share as `-t` to the CLI they spawn, `rust/tests/lib.rs`), and the C library's test, which has no `-t`, runs under `taskset` on `KZ_TEST_THREADS` CPUs. The bytes never depend on it. `-t` bounds the codec threads, not every thread: compress also runs a fixed input reader and a names thread beside them, so a `-t 1` compress can have 2 or 3 runnable threads. The stress cases keep the thread counts they test.
- `TEST_LIGHT=1` is the serial, lightweight mode; `TEST_JOBS=1 TEST_LIGHT=1 make test` is the run for a small machine. It skips the stress cases and prints one `skipped (TEST_LIGHT): CASE` line for each, counted as skipped, not passed. The stress cases are of two kinds. Thread counts: the multi-thread equivalence checks (smoothing at 2, 4 and 8 threads against 1; `--db --aln` and `--read-blocks 7 --smooth` at `-t 4`; decode ahead at `-t 4` against `-t 1`, with its 120000-read input), `-t 8` on 200000 one-row sequences, `-t 1024` accepted, `OMP_NUM_THREADS` ignored at the default thread count and at `--threads 2`, diff.sh's default count (no `-t`: min(the affinity mask's CPUs, 8) in both) and its `-t4 compress` and `--threads 4` argument forms, the concurrent readers (up to 8), the full body temp file at `-t 2` (run.sh and diff.sh), and the BGZF 4 GiB and 300000-block inputs under a 1.5 GB address limit. Size: every generated input above 20000 rows (run.sh: the 80000-row long temp dir, the 100000-row file-size limit, the 60000-row label cache, the 400000-row full staging filesystem; diff.sh: the 300000-row many-BGZF-blocks file, the 400000-row damaged decompress, the 65537-row forged chunk tables; fuzz.sh: the forged-index stage and its 140000-row file), `-e high` on large inputs (the six contigs of 20000 to 40000 rows and the 60000 short reads, in run.sh and diff.sh), and featureID 2^26-1, whose `--db` tables take about 1 GB for one row. fuzz.sh also runs 10 cases per generator instead of 60. What runs: every other check. The largest generated input kept is 20000 rows (and one row with a 2.1 MB label); `-e high` runs only on the fixtures (at most 175 rows) and diff.sh's 200-read tail inputs. In every mode, every karyozip child gets `-t KZ_TEST_THREADS` (1 under `TEST_JOBS=1`) unless the case tests a thread count, the Rust integration tests' CLI child included (diff.sh's argument-form cases put it first, where a later `-t` still wins; the forms that end before any thread starts, no arguments, `--help` and `--version`, run with none, and run.sh's light default-count check runs the bare binary under `taskset` on one CPU), and two outputs are compared by running their karyozip one after the other into files (`seqcmp` in tests/lib.sh), never through `cmd | cmp - <(cmd2)`. Measured in two 1-CPU, 8 GB Slurm jobs, sampling every 0.2 s: 154 to 170 s, at most 3 runnable threads (one `-t 1` compress and its reader and names threads), one karyozip at a time (a second process appeared only as LeakSanitizer's tracer of the same `karyozip-san`, in the leak stage), and at most 302 MB MaxRSS (Slurm; the largest process is an `-e high` compress of a fixture, about 250 MB); before these limits the same run took 278 s with 4 runnable threads, 2 karyozip processes and 1.06 GB summed RSS. With 8 CPUs allocated and `KZ_TEST_THREADS=1` the budget still holds, since no case falls back to the default count: 357 s, one karyozip at a time, and at most 4 threads in any karyozip process (a `-t 1` compress: main, one codec thread, reader, names), sampled every 0.2 s and, for diff.sh's short-lived argument-form cases, counted with `strace -f` (up to 13 threads each before they got the budget). The stress cases use the thread counts and sizes they test, not the budget, so run them in a bounded Slurm job, not on a laptop.
Four variables are read by the release binaries, C and Rust, for the signal tests only: `KZ_TEST_SIG_PAUSE_MS=N` makes a SIGINT/SIGTERM/SIGHUP handler sleep N ms right after it decides to exit or to defer (read once at startup); `KZ_TEST_SWAP_PAUSE_MS=N` makes an output commit (`--out-dir`'s set of BEDs, `-o`'s rename) sleep N ms before its window opens and after each rename, so run.sh can signal just before and inside the window; `KZ_TEST_STAGE_PAUSE_MS=N` makes `--out-dir` sleep N ms after it creates each staged BED and before it records the name for the signal handler, the window in which a signal once leaked that file; and `KZ_TEST_STAGE_ENTRY_PAUSE_MS=N` makes `--out-dir`'s staging sleep N ms once it has told the handlers it is staging and before it looks for a handler that has claimed the exit, where the two once waited for each other forever; as that sleep starts it writes `karyozip: test: staging entry paused` to stderr in one write(2), so run.sh signals only once the process is inside it. They only add sleeps at those fixed points (and that one line on stderr), so they change timing, never the bytes written; unset (or 0), nothing sleeps. They are not for use outside the test suite.
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

## The Rust library

`rust/` is also a library crate, `karyozip`, for a program that writes `.kst` itself (a Cargo git dependency on the `rust/` subdirectory). The write path returns errors; it never exits the process or panics on anything a caller can cause:

```rust
use karyozip::{writer::{Options, Writer, Stats}, Error, ErrorKind, KzResult};
pub fn Writer::new(opt: Options, tmpdir: &Path) -> KzResult<Writer>;
pub fn Writer::row(&mut self, name: &[u8], run_len: u64, label: &[u8], start: u64) -> KzResult<()>;
pub fn Writer::push(&mut self, name: &[u8], run_len: u64, label: &[u8], start: u64, same: bool) -> KzResult<()>;  // row without its byte checks
pub fn Writer::close(self, out: &mut impl Write) -> KzResult<Stats>;
pub fn Writer::set_aln(&mut self, a: aln::Aln);  pub fn Writer::set_feathash(&mut self, h: u64, cart: bool);
pub fn db::Feats::load(path: &str) -> KzResult<Feats>;   pub fn db::Feats::load_hier(&mut self, path: &str) -> KzResult<()>;
pub fn db::encode(files: &[String], feats: &Feats, w: &mut Writer) -> KzResult<bool>;
pub fn aln::Aln::{load_ref, load_paf, load_aln, load_ordinals}(&mut self, ...) -> KzResult<()>;
```

`KzResult<T>` is `Result<T, karyozip::Error>`; an `Error` has `message()` (the C's words, which the CLI prints after `karyozip: `, then exits 1), `kind()` (`Input`, `Io` or `Internal`) and implements `std::error::Error`.
A failing `out` is `write failed`, a temp directory that cannot take the staging files `cannot create a temp file under DIR (see --tmp-dir)`, at `new` or at `close`. The staging files are unlinked when made, so neither a failed close nor a dropped `Writer` leaves a `karyozip.*` file, and the names threads are joined before either returns.
Writers may run one after another or at once on several threads; they share rayon's global pool. Process-wide state left: the OpenMP team size, which a close sets and restores (the speed of a writer beside it changes, never its bytes), and the allocation cap, which turns on once any file is read, as in the C. The writer prints its decisions (`codec:`, `layout:`, `names:`) on stderr.
The reader path still exits on a refused or corrupt file, as the CLI does: `reader::Reader`, smoothing, query, index, `decompress --db`, `aln-pack`, and a `.kst` read as `--ref` or as a `-d` input. `rust/tests/lib.rs` checks the refusals' words, a missing, read-only or vanished temp directory, an output that fails after N bytes, writers in sequence and on two threads, and `db::encode` over per-set `.kst` files.

The Rust CLI restores SIGPIPE's default, as the C runs: `karyozip ... | head` ends quietly (killed by SIGPIPE, 141 in the shell, or 0 when its last write landed first), never with a panic or `write failed`.

## Vendored code

htscodecs, bzip3, libdeflate, liblzma and libzpaq under `c/` are pinned; each `PROVENANCE` names the release and every local patch.
The vendored bzip3 also takes each block's size from the caller (`bz3_compress_cuts`), so blocks end at read starts; its frame already recorded each block's size, so the decoder is unchanged.
The vendored bzip3 runs its independent blocks in parallel, writing the same bytes; `c/bzip3/lastblock_check.c` guards the last-block fix that bzip3 1.5.4 took from this tree.
