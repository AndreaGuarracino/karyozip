# C against Rust

Two implementations of the same format, in `c/` and `rust/`. They are
equivalent, not similar: `tests/diff.sh` asserts that every encode is
byte-identical, every decode gives the input back, and every refusal uses
the same words. Run it with `make diff`.

```
codec versions: zstd=1.5.7 libbz3=1.5.3
codecs: every candidate round-trips
encode: 342 forms identical
decode: 95 forms give the input back
cross-decode: the C reads the Rust's files
region: 110 forms identical
index: byte-identical trailer / index cross-read: 1 forms identical
sorted decode: 12 forms identical
reader flags: 16 forms identical
smooth: 24 forms identical
--db encode: identical / --db decode: identical, restores all six inputs / --db --hier: 2 forms identical
--aln and aln-pack: 11 forms identical
refusals: same message from both
gzip input: 10 forms match the plain track
cli: 74 forms give the same words and code
DIFF_PASS
```

## Size

|                        |     C |  Rust |
|------------------------|------:|------:|
| lines, excluding vendored C | 4,634 | 4,740 |
| files                  |     2 |    11 |
| `unsafe` blocks        |     — |    14 |
| stripped binary        | 382 KB | 2,269 KB |

The two are the same length within 1%. Some of that is real, the parallel arrays and the
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
vendored C codecs and checks the standard fds before the runtime touches
them. Nothing else in the crate is unsafe.

## Speed

### A real track: ONT long reads, 2.20 GB of `.bed.gz`

The KaryoScope `chromosome` track of an ONT_frag sample against CHM13 v2:
293,614,022 rows, 9,387,314 reads, 1,121 blocks, 17.17 GB of BED, read
straight from gzip. All 13 encode forms byte-identical between C and Rust
(default, `--block 4096`, `--block all`, `--level 3`, `--level 22`, `--names
keep`, `--names drop`, `--tail 0`, `--tail 1`, `--tok-chunk 250000`,
`--sort`, `--sort=256M`, `--slow`), round trip exact, identical at 1, 4 and
16 threads, `index`, `region` and `decode --raw` identical.

Fastest of 3 interleaved, 8 threads, on the cluster login node:

|        | C | Rust |
|--------|--:|-----:|
| encode | 151.6 s, 16.3 GB | 157.9 s, 17.1 GB |
| decode |  20.4 s, 1.54 GB |  28.2 s, 1.24 GB |

The Rust encode was 2.3x slower before this track was measured: it ranked
the three names codecs over all 164 MB of read names where the C ranks them
on a 200,000-name sample, and it compressed the tokeniser chunks serially.
Both fixed, its names phase is now the faster one (92 s against 107 s) and
its parser too (42 s against 55 s); `parse:append`, the serial row path, is
where it still loses (14.0 s against 8.2 s). The decode number above
predates two fixes found afterwards, the tokeniser chunks decoded serially
and one `Vec` per sequence name; on a 5M-name track those took the Rust
decode from 3.12 s to 1.76 s against the C's 1.52 s. The big track is to be
re-measured.

Compression: 17.17 GB of BED to 0.552 GB of `.kst`, 31.1x over the plain
BED and 3.98x over gzip, 1.880 bytes per row. Sequence names are 30% of
the file (`--names drop` gives 0.387 GB).

### Small files

Fastest of 9 to 11 interleaved runs. The small case is 4,000,000 rows,
99,354,996 B of BED in; the large one 40,000,000 rows over 658 sequences,
1,390,855,987 B in, 27,940,756 B of `.kst` out.

|          | C, 4M | Rust, 4M | C, 40M | Rust, 40M |
|----------|------:|---------:|-------:|----------:|
| encode   | 0.31 s, 55 MB | 0.35 s, 73 MB | 1.78 s, 43 MB | 2.02 s, 78 MB |
| decode   | 0.39 s, 9 MB | 0.40 s, 12 MB | 1.33 s, 8 MB | 1.26 s, 10 MB |
| `--slow` | | | 32.8 s, 103 MB | 36.7 s, 120 MB |
| build from clean | 14.2 s | 30.4 s | | |

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
