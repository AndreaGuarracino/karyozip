# C against Rust

Two implementations of the same format, in `c/` and `rust/`. They are
equivalent, not similar: `tests/diff.sh` asserts that every encode is
byte-identical, every decode gives the input back, and every refusal uses
the same words. Run it with `make diff`.

```
codec versions: zstd=1.5.7 libbz3=1.5.3
codecs: every candidate round-trips
encode: 342 forms byte-identical
decode: 95 forms give the input back
cross-decode: the C reads the Rust's files
region: 110 queries identical
index: byte-identical trailer
sorted decode: 12 forms identical
reader flags: 16 forms identical
smooth: 27 forms identical
--db encode / decode: identical, restores all six inputs
--aln and aln-pack: 11 forms identical
refusals: same message from both
```

## Size

|                        |     C |  Rust |
|------------------------|------:|------:|
| lines, excluding vendored C | 4,207 | 3,170 |
| files                  |     2 |     8 |
| `unsafe` blocks        |     — |    14 |
| stripped binary        | 358 KB | 1,896 KB |

The Rust is 25% shorter. Some of that is real, the parallel arrays and the
manual memory management that `Vec` and ownership remove; some is that the
C carries its own hash tables and buffer type where the Rust uses the
standard library.

All 14 `unsafe` blocks are in `src/ffi.rs`, the one module that calls the
vendored C codecs. Nothing else in the crate is unsafe.

## Speed

Interleaved medians of five, on a realistic 40M-row track built by tiling a
real ONT chromosome track, 2.3 GB of BED in, 22.3 MB of `.kst` out.

|          |    C |  Rust |
|----------|-----:|------:|
| encode   | 4.23 s, 273 MB | 5.32 s, 405 MB |
| decode   | 0.74 s, 12 MB | 0.98 s, 17 MB |
| build from clean | 14.2 s | 30.4 s |

The C is 26% faster on encode and 32% on decode. Both numbers were much
worse before profiling: the Rust encode started at 55% slower because its
parser was serial, and decode at 2.5x because `to_string()` allocated a
`String` per coordinate, which was a quarter of decode's time.

## Deployment

This is where they differ most, and it is not close.

```
C:    libzstd.so.1 => /home/guarracino/miniforge3/lib/libzstd.so.1
      libgomp.so.1 => /home/guarracino/miniforge3/lib/libgomp.so.1
      libm, libc, libpthread, libdl
Rust: libgcc_s, libm, libc
```

The C binary needs libzstd and libgomp from a conda environment, at the
path burned into its rpath. Move it to a machine without that path and it
does not start. The Rust binary needs only what every glibc system has,
because `zstd-sys` bundles zstd and rayon needs no runtime.

`cargo build --target x86_64-unknown-linux-musl` would give one static
binary with nothing to install. It was not tested here: `zstd-sys` runs
bindgen, which needs the musl headers, and they are not on this machine.

## Which bugs the compiler would have caught

Of the seven defects audit round 32 found in the C:

- **argv read past the end** after the flag loop compacted it. Rust
  indexes a `Vec`, so this is a panic, not a segfault, and the parser is
  written against a positional list that cannot run off the end.
- **`refend_at`'s cursor was a plain `static`**, shared across threads.
  In Rust that state is a field on the reader, and sharing it across
  threads would not compile.
- **`paycomp` computed a 32-bit bound that wrapped**. Rust's release
  profile wraps too, so this one is not caught; the port computes it in
  64 bits because the C now does.
- **`die()` from inside an OpenMP region**. Rust returns `Result` from the
  parse workers, so there is no equivalent.
- The other three, `region` closing stdout, `strtoul` without an endptr,
  and a test that passed for the wrong reason, are design mistakes that
  neither language prevents.

## What the port cost

Three real bugs, all mine, all caught by the equivalence gate rather than
by reading:

- decode walked the block from its head for every row, O(n squared), the
  same bug the C fixed when smooth went 31.5 s to 3.5 s
- the tokeniser returns names NUL-separated, not newline-separated, which
  only showed at `--level 3` where the tokeniser wins the name stream
- no block cache, so a track whose sequences share a block decompressed it
  once per sequence: smooth was 7.6x slower until that was added

The gate found each of them within a minute of the code being written.
