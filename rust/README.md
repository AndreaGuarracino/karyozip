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

`src/ffi.rs` is the only `unsafe` in the crate.
Every raw pointer and every length check lives there; callers see slices and
`Option<Vec<u8>>`. Beyond that, `unsafe` needs a measured win recorded in
`notes/log.md`, not an argument.

## State

Milestone 1 done: the codec layer. Five zstd levels, rANS order 0 and 1,
rANS with RLE and with pack, arith order 0 and 1, and bzip3 at two block
sizes all give byte-identical output from Rust and from C on a 23.8 MB
slice of a real track, and all eight round-trip.

Next: the writer, then the reader, then the CLI.
