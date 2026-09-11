# karyozip

Compresses KaryoScope annotation tracks losslessly, and reads them back without decompressing the whole file.
A `.bed.gz` track becomes one `.kst` file, several times smaller; `decode` returns the original BED byte for byte.

## Install

```
make            # builds c/build/karyozip and rust/target/release/karyozip
```

Either binary works; they write and read the same files.
The C needs zstd, liblzma, zlib and an OpenMP compiler; the Rust needs cargo, liblzma and libgomp.

## Quick start

```
karyozip encode < sample.chromosome.presmoothed.bed.gz > sample.kst
karyozip decode sample.kst > sample.bed                  # the original BED, byte for byte
karyozip region sample.kst chr1_hap1 1000000 1200000     # the runs overlapping a range
karyozip smooth sample.kst hierarchy.tsv chromosome      # KaryoScope's smoothing, from the presmoothed track
```

`encode` reads plain or gzipped BED from stdin, so it can sit behind `karyoscope annotate` in a pipe.
Every other text input may be plain or gzipped too.

## Commands

```
karyozip encode [flags] < in.bed[.gz] > out.kst
karyozip decode F                               the whole track as BED
karyozip region F SEQ                           one sequence
karyozip region F SEQ POS                       the run covering POS
karyozip region F SEQ START END [HALO]          the runs overlapping START..END, plus HALO runs each side
karyozip smooth F HIER FS [SEQ [START END]]     KaryoScope's hierarchy smoothing of feature set FS
karyozip index F                                adds a name index, so region finds a read fast on tracks with many reads
karyozip aln-pack P.paf REF.bed|REF.kst OUT     cuts a PAF to what --aln needs
```

## Size or speed

| you want | use |
|---|---|
| the smallest file, to archive | `encode`, the default |
| a smaller file still, whatever the time | `encode --slow` |
| a quicker encode, a larger file | `encode --fast` |
| fast `region` queries on a big track | `encode --queryable`, 3.1% larger in our test |

The default looks at the first 16 MB of a track to choose how to store it, `--fast` at 4 MB, `--slow` at all of it.
On our four biggest test tracks, up to 84,741,909 rows, `--slow` took 4.5 to 6.6 times as long as the default and came out 0% to 3.1% smaller.
`--queryable` combines with any of the three.

## Feature sets, smoothing and alignments

A KaryoScope run writes one BED per feature set, presmoothed and smoothed.

- **Smoothed tracks need not be stored.** `karyozip smooth` rebuilds a smoothed track from the presmoothed one, row for row as KaryoScope writes it. The exception is the v1 `region` and `acrocentric` sets, whose hierarchy changed after they were written: store those.
- **All feature sets in one file.** `--db features.tsv` merges one sample's per-set BEDs into a single, smaller track. Smoothed sets merged this way also need `--hier hierarchy.tsv`, on encode and on every read.

```
karyozip encode --db features.tsv s.acrocentric.bed.gz s.chromosome.bed.gz ... > s.all.kst
karyozip decode --db features.tsv --out DIR s.all.kst       # every set's BED back
karyozip decode --db features.tsv --set gene s.all.kst      # one set
```

- **Against a reference.** With the reads aligned to a reference (a PAF) and the reference's own track, `--aln P.paf --ref ref.bed` stores only what differs from the reference. Every read of such a file takes the same `--aln` and `--ref`.

## Read names

Read names are stored in the file by default.

| flag | effect |
|---|---|
| `--names drop` | store none; the file is smaller and a read is known by its number |
| `--names drop:F` | write the names to `F` beside the track |
| `--names F` | on decode, region and smooth: take names from `F`, a `drop:F` file, a FASTA, a FASTQ or a list, one a line |
| `--names ordinal` | on decode: accept a track stored with `--names drop` and print the numbers |
| `--names keep` | store names even when they are 0, 1, 2, ... |

`--sort` groups records by their labels, which pays on short reads stored with `--names drop`; a decode then returns the rows in that order.

## All flags

| flag | on | effect |
|---|---|---|
| `--fast`, `--slow`, `--queryable` | encode | see Size or speed |
| `--block N`, `--block all` | encode | rows per block by hand; `all`, the default, is one block |
| `--level N` | encode | zstd level, 1 to 22, default 17 |
| `--sort`, `--sort=MEM` | encode | group records by labels; spill to disk past `MEM` |
| `--names ...` | encode, readers | see Read names |
| `--db features.tsv`, `--hier H`, `--set A,B` | encode, readers | see Feature sets |
| `--aln P.paf --ref R`, `--aln-names F` | encode, readers | see Feature sets; `--aln-names`: read N is line N+1 of `F` |
| `--out DIR`, `--raw` | decode | one BED per set under `--db`; `--raw` prints what is stored, for debugging |
| `--max-gap N`, `--v1` | smooth | KaryoScope's smoothing options |
| `--threads N`, `-t N` | all | default `OMP_NUM_THREADS`, else the cores, up to 8 |
| `--tmp-dir DIR` | all | where temporary files go, default `TMPDIR`, else /tmp |
| `--profile`, `--help`, `--version` | all | time per phase; usage; library versions |
| `--tail 0\|1`, `--tok-chunk N`, `--id-codec N`, `--len-codec N` | encode | pin a choice the encoder otherwise makes itself |

Flags work on either side of the command name.

## Supplementary: technical notes

**Two implementations.** `c/` and `rust/` write the same bytes, print the same `--help` and refuse the same inputs with the same words; `rust/build.rs` takes the usage text from `c/karyozip.c`.
`make test` (64.5 s) runs the test suite, sanitizers and fuzzing on both and compares them encode by encode.
The C finds zstd, liblzma, zlib and libgomp at the path in its rpath; the Rust bundles zstd.

**How encode chooses.** The three modes are one rule with two settings: `--fast` tries 8 of the 51 codecs on 4 MB, the default all 51 plus zstd at 17 and 22 on 16 MB, `--slow` the same on the whole file.
The smallest wins, with no margin, and must read back what it was ranked on.
A one-block file also chooses its layout on the same sample: lengths and labels apart, run counts heading each read's lengths, lengths grouped by label, or labels interleaved with lengths; every layout and codec is ranked in one parallel pool.
Names are ranked the same way among zstd, integer deltas, field deltas and the CRAM name tokeniser.
Encode prints what it chose on stderr, on its `layout:` and `names:` lines.
`--queryable` asks for 262144-row blocks and a block index and drops the codecs slow to decode: bzip3, arith and LZMA.

**Name index.** `karyozip index` appends a 32-bit fingerprint and a bit-packed ordinal per name, 2,044,644 B for 300,000 reads; a fingerprint match is confirmed against the decoded name.

**Benchmark.** `tests/bench.sbatch` runs `tests/bench.sh` on the cluster at three sizes over the same rows, `small`, `wide` and `big`, checking every row's round trip, C against Rust bytes and smoothing against KaryoScope's, and recording the layout, codecs and phase times each encode chose.

**Library.** `make` also builds `c/build/libkaryozip.a`, six calls in `c/karyozip.h` (`kt_open`, `kt_block`, `kt_level`, `kt_row`, `kt_close`, `kt_version`), for an annotator that writes `.kst` directly.
Rows must arrive grouped by sequence, each sequence once, starting at 0, gapless; anything else is refused with the CLI's message.

```c
kt_writer *w = kt_open(f);
kt_row(w, name, nlen, run_len, label, llen, start);   /* once per run */
kt_close(w);
```

**Vendored code.** htscodecs and bzip3 under `c/` are pinned; `c/htscodecs/PROVENANCE` and `c/bzip3/PROVENANCE` name the release and every local patch.
Upstream bzip3 1.5.3 loses half the data on an input that is an exact multiple of the block size, which `c/bzip3/lastblock_check.c` reproduces; the vendored copy also runs its independent blocks in parallel, writing the same bytes.
