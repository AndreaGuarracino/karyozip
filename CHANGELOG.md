# Changelog

No release has been tagged yet: the version `0.2.0` in `rust/Cargo.toml` and `c/karyozip.c` points at no tag.

## Unreleased

The current state is version 0.2.0, file format `KST00001` ([FORMAT.md](FORMAT.md)).

- Commands `compress`, `decompress`, `query`, `index` (the `.kzi` name index) and `aln-pack` (the `KSTAL1` packed alignment).
- Compression modes: standard, db-informed (`-d`), alignment-informed (`-a`/`--ref`), and the two combined. `--smooth` rebuilds KaryoScope's smoothed track from the presmoothed one.
- Effort levels `low`, `medium` and `high`. The codecs are zstd, htscodecs rANS and arith, bzip3, LZMA2 and zpaq. Names are stored as zstd text, integer deltas, CRAM-tokenised or field deltas.
- A C implementation in `c/` and a Rust one in `rust/` write the same bytes. The C also builds the library `libkaryozip.a`.
- Format: the magic `KST00001` stayed the same while flags and codecs were added, so a file from an earlier commit is not guaranteed to read. See the Compatibility section of FORMAT.md.
