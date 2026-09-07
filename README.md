# karyozip

Reads and writes the **KaryoScope Track (KST)** format: lossless, compact
storage and O(1) random access for k-mer annotation tracks. Files are
`.kst`, and there is only ever one of them: `karyozip index F` appends the
index to the track itself.

Two implementations, `c/` and `rust/`, which produce the same bytes.
`make diff` proves it. `COMPARISON.md` has the numbers.

```
make            # both binaries
make test       # the C suite: ~256 checks, LeakSanitizer, five fuzz stages
make diff       # the C and the Rust must agree on the small fixtures
make matrix     # and on two real HPRC assembly cells, every set, strategy, flag and reader
make rust-test  # the C's own suite and fuzzer, run against the Rust
```

`tests/bench.sh CELLDIR PAF FEATURES HIER REFDIR` prints the benchmark
table, one TSV row per set, smoothing, strategy, seekability and
implementation: sizes, ratios, round trip, encode/decode/smooth time and
memory, and query latency. `tests/mkfix.py` cuts a real cell down to a
fixture on the cluster.

## Commands

Every text input may be plain or gzipped, and is inflated on the fly: the
BED on stdin, the `--db` BEDs, `--ref`, the PAF, `--aln-names`, `--names F`
and the hierarchy. A bgzipped file works too, being many gzip members.

```
karyozip encode [flags] < in.bed        > out.kst
karyozip encode [flags] < in.bed.gz     > out.kst
karyozip decode F | --db F.tsv [--hier H] [--out DIR] [--set A,B] F
karyozip region F SEQ | SEQ POS | SEQ START END [HALO]
karyozip smooth F HIER FS [SEQ [START END]]
karyozip index  F                       appends an index to the track
karyozip aln-pack P.paf REF.bed OUT     packs a PAF to only what --aln reads
```

## Flags

| flag | effect |
|---|---|
| `--block N`, `--block all` | rows per block, default 262144; `all` is one block, no random access, smallest file |
| `--level N` | zstd level, 1 to 22, default 17 |
| `--slow` | every codec over every block: smaller, slower both ways |
| `--sort`, `--sort=MEM` | order by label tuple, then lengths; spill past `MEM` |
| `--names keep` | store names even when they are 0,1,2,... |
| `--names drop`, `--names drop:F` | store none, or write them beside the track |
| `--db features.tsv` | one BED per feature set, or one merged featureID track; every reader takes it, and region then prints one label per set |
| `--hier hierarchy.tsv` | with `--db`: a label tuple outside features.tsv, as a smoothed track carries, gets 2^26 plus its rank in the sorted per-set vocabulary (the features.tsv column, the hierarchy nodes of that set, `novel`), the same id in every file. A presmoothed track is unchanged; a track that used one such id needs the same two files to read, and smooth takes its HIER as the `--hier` |
| `--aln P.paf --ref R` | do not store what the reference already gives; every reader takes the same pair. `R` is the reference BED, or that BED encoded plainly as a `.kst`, which a query reads one contig of instead of parsing the whole text |
| `--aln-names F` | ordinal N is line N+1 of `F` |
| `--names F`, `--names ordinal` | on the readers |
| `--raw`, `--out DIR` | on decode |
| `--set A,B` | decode and region under `--db`: which feature sets, default all. decode without `--out` and region print the picked labels as columns, merged where they all agree, so one set gives that set's own track |
| `--max-gap N`, `--v1` | on smooth; `SEQ START END` smooths the whole sequence and prints the overlap, so it equals the full smooth |
| `--threads N`, `-t N` | 1 to 1024. Default `OMP_NUM_THREADS`, else cores capped at 8 |
| `--tmp-dir DIR`, `--profile`, `--help`, `--version` | any command |

Every flag is read before the command runs, so it works on either side of
the command name. Two environment variables are read and neither changes
the stored bytes: `TMPDIR` when `--tmp-dir` is not given, and
`OMP_NUM_THREADS` when `--threads` is not given, so a SLURM allocation is
honoured.

`--sort` reorders records, so a decode of a sorted track returns the rows
in sorted order. Compare such a decode against a sorted copy of the input.

## Library

An annotator can write a `.kst` directly instead of emitting a BED for
karyozip to re-parse. `make` builds `c/build/libkaryozip.a`; the API is six
calls in `c/karyozip.h`.

```c
kt_writer *w = kt_open(f);
kt_row(w, name, nlen, run_len, label, llen, start);   /* once per run */
kt_close(w);
```

Rows must arrive as a track is already laid out: grouped by sequence, each
sequence appearing once, starting at 0, gapless. Those three properties
were checked on 264,914,318 real rows across all four data classes with
zero violations, and anything else is refused with the message the CLI
gives.

## Dependencies

zstd, and the vendored htscodecs and bzip3 under `c/`. Both are pinned:
`c/htscodecs/PROVENANCE` and `c/bzip3/PROVENANCE` name the release and
every local patch. The bzip3 patches matter: upstream 1.5.3 loses half the
data on an input that is an exact multiple of the block size, and
`c/bzip3/lastblock_check.c` reproduces it.
