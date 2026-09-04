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
cli: 69 forms give the same words and code
```

## Size

|                        |     C |  Rust |
|------------------------|------:|------:|
| lines, excluding vendored C | 4,218 | 3,637 |
| files                  |     2 |    10 |
| `unsafe` blocks        |     — |    14 |
| stripped binary        | 358 KB | 2,457 KB |

The Rust is 14% shorter. Some of that is real, the parallel arrays and the
manual memory management that `Vec` and ownership remove; some is that the
C carries its own hash tables and buffer type where the Rust uses the
standard library. The command line is `clap`, which costs 589 KB of binary
and buys back nothing in lines, because the C's exact refusal wording still
has to be reproduced by hand around it. What it does buy is that the flag
table, the scope rules and the parser cannot drift apart.

The usage text is not copied. `rust/build.rs` lifts the `USAGE` literal out
of `c/karyozip.c` at build time, so `--help` is the same bytes from both by
construction. The htscodecs and bzip3 version strings come the same way,
from `c/htscodecs/version.h` and `c/Makefile`.

All 14 `unsafe` blocks are in `src/ffi.rs`, the one module that calls the
vendored C codecs. Nothing else in the crate is unsafe.

## Speed

Fastest of 11 interleaved runs, on 4,000,000 rows, 99,354,996 B of BED in,
7,586,839 B of `.kst` out, byte-identical from both.

|          |    C |  Rust |
|----------|-----:|------:|
| encode   | 0.31 s, 55 MB | 0.35 s, 73 MB |
| decode   | 0.39 s, 9 MB | 0.40 s, 12 MB |
| build from clean | 14.2 s | 30.4 s |

Decode is a tie. Encode is 13% apart, and the phase breakdown says where:
the Rust searches codecs faster (0.120 s against 0.123 s at equal work,
after both were fixed) and appends rows slower (0.083 s against 0.067 s).

Both were much slower before profiling. The search compressed two bzip3
candidates that the seekable-track rule then refused, and both ran their
zstd baselines serially ahead of the parallel sweep; fixing those took the
C from 0.63 s to 0.31 s and the Rust from 0.77 s to 0.35 s, with no change
to a single stored byte. On the Rust side the parser also built its row
array twice, and the label table hashed one-byte keys with SipHash.

The memory difference is glibc, not the program. Rust peak scales with
threads (43 MB at `-t 1`, below the C's 53 MB; 69 MB at `-t 8`) because
every codec call allocates its own `Vec`, and each rayon worker keeps its
own malloc arena. `MALLOC_ARENA_MAX=1` brings it to 57 MB and costs 13% of
the encode. That trade is left alone: the tool sets no environment
variable, and the memory is bounded either way.

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
