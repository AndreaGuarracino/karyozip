# Changelog

## 0.3.0 (2026-11-03)

The first tagged release. File format `KST00001` ([FORMAT.md](FORMAT.md)); files written by 0.3.0 read with this release, and readers now refuse header flags they do not know instead of misreading them.

### Fixed
- `-o`, `aln-pack`, `--out-dir` and the `.kzi` index never truncate or replace a good output on failure: outputs are staged under exclusive random names, synced and closed with checks, then renamed in.
- `--out-dir` replaces its set of BEDs all or nothing with respect to SIGINT/SIGTERM/SIGHUP and failed renames; backups have unique names, and a run warns about `*.bed.old.*` leftovers of an earlier crash.
- BGZF input: the extra field and each block's ISIZE are validated before any allocation (a forged header could read past the input, or ask for many GB).
- A ten-byte varint that overflows 64 bits is refused (it was silently truncated).
- The Rust compressor no longer hangs when writing its temporary body fails.
- The C library resets all encoder state in `kt_open`; `karyozip.h` documents the real contract (one block by default, not thread-safe, exits on error).
- C: the lazy `.kst` reference fill is race-free; host byte order no longer reaches the file; 32-bit counters on very large inputs; the allocation cap applies on every thread; the default thread count follows the CPU affinity mask (Slurm allocations).
- Rust: `query --smooth` without `-r` is refused, as in the C; a packed alignment's magic is read before the whole file; a closed stdout pipe ends quietly, as the C.
- PAF and packed-alignment coordinates that leave 64 bits are refused (`PAF coordinates overflow 64 bits`, `corrupt packed alignment: a coordinate overflows`); they used to wrap or, in a Rust debug build, panic.
- A signal while `--out-dir` creates its staged BEDs no longer leaves one behind.
- A `.kst` given to `-d` through a pipe is refused with the regular-file message whatever its ninth byte; the C name table refuses 2^32 or more sequences, as the Rust's does.
- A second SIGINT, SIGTERM or SIGHUP arriving while the first one's handler runs no longer hangs the process: the handlers are installed with `sigaction`, each masking the other two.
- Rust: a header-only `features.tsv` no longer panics (a row labelled `novel` is accepted, as in the C).

### Added
- The Rust library's write path (`Writer`, `db::encode`, the `features.tsv` and alignment loaders, including a packed alignment's tokenised names) returns `Result<_, karyozip::Error>` instead of exiting.
- `compress -d` takes each feature set as its BED or its plain `.kst`, mixed as you like, for the same bytes. BED inputs may also come through a pipe, a FIFO or `<(...)`; a `.kst` input must be a regular file. A `.kst` is known by its magic and a first line no BED row could be, so a BED whose first sequence name begins with `KST00001` is a BED whatever its form.
- FORMAT.md: a specification of `KST00001`, the `.kzi` index and `KSTAL1`.
- Tests: a shared harness with pass/fail counts, `CASE`, a total `TEST_JOBS` budget, `TEST_LIGHT=1` for small machines, seeded fuzzing, Rust unit tests, run by `make test`, and a CI workflow template in `ci/`.

### Changed
- The C builds inside an activated conda environment; the usage text and the liblzma source list are shared files, no longer scraped from the C; the build checks for zstd 1.5.7.
- `-o` together with `--out-dir`, and an empty `--out-dir`, are refused.

## 0.2.0 (untagged)

- Commands `compress`, `decompress`, `query`, `index` (the `.kzi` name index) and `aln-pack` (the `KSTAL1` packed alignment).
- Compression modes: standard, db-informed (`-d`), alignment-informed (`-a`/`--ref`), and the two combined. `--smooth` rebuilds KaryoScope's smoothed track from the presmoothed one.
- Effort levels `low`, `medium` and `high`. The codecs are zstd, htscodecs rANS and arith, bzip3, LZMA2 and zpaq. Names are stored as zstd text, integer deltas, CRAM-tokenised or field deltas.
- A C implementation in `c/` and a Rust one in `rust/` write the same bytes. The C also builds the library `libkaryozip.a`.
- Format: the magic `KST00001` stayed the same while flags and codecs were added, so a file from an earlier commit is not guaranteed to read. See the Compatibility section of FORMAT.md.
