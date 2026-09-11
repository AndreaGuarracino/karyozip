# karyozip

Reads and writes the **KaryoScope Track (KST)** format: lossless, compact storage and O(1) random access for k-mer annotation tracks.

Two implementations, `c/` and `rust/`, which produce the same bytes and print the same `--help`.
`make test` proves it; `COMPARISON.md` has the numbers.

```
make            # both binaries
make test       # 64.5 s: the suite, sanitizers and fuzz on both binaries, C against Rust
```

`tests/bench.sbatch` runs the benchmark on the cluster at three sizes over the same rows, `small`, `wide` and `big`, every row checked for its round trip, for C against Rust bytes and for smoothing against KaryoScope's own, and recording the layout and codecs each encode chose.
`tests/bench.sh CELLDIR PAF FEATURES HIER REFDIR [C] [RUST] [KIND:PART[:SEEK[:IMPL]]]` prints its table, one TSV row per feature set, smoothing, storage form, seekability, mode and implementation.

## Commands

Every text input may be plain, gzipped or bgzipped, and is inflated on the fly: the BED on stdin, the `--db` BEDs, `--ref`, the PAF, `--aln-names`, `--names F` and the hierarchy.

```
karyozip encode [flags] < in.bed[.gz]   > out.kst
karyozip decode F | --db F.tsv [--hier H] [--out DIR] [--set A,B] F
karyozip region F SEQ | SEQ POS | SEQ START END [HALO]
karyozip smooth F HIER FS [SEQ [START END]]
karyozip index  F                       appends a name index to the track
karyozip aln-pack P.paf REF.bed|REF.kst OUT   packs a PAF to only what --aln reads
```

## Encode modes

With no flag the file is one block, the smallest, with no random access.
The three search modes are one algorithm with two settings, which codecs and how much data decides:

| mode | codecs | decides on |
|---|---|---|
| `--fast` | 8 of the 51 codecs, plus zstd | the first 4 MB of each stream |
| default | every codec, plus zstd at the level and at 22 | the first 16 MB |
| `--slow` | every codec, plus zstd at the level and at 22 | the whole file |

In all three the smallest wins, with no margin, and must read back what it was ranked on.
A one-block file also picks its layout on the same sample: lengths and labels apart, run counts heading each read's lengths, lengths grouped by label, or labels interleaved with lengths.
Read names get the same rule: zstd, integer deltas, field deltas or the CRAM name tokeniser, ranked on the mode's sample.
Encode says on stderr what it chose, on its `layout:` and `names:` lines; `--profile` adds the time each phase took.

`--queryable` asks for 262144-row blocks and a block index, and keeps only codecs that decode fast: no bzip3, arith or LZMA.
It composes with the three modes.

## Flags

| flag | effect |
|---|---|
| `--queryable` | 262144-row blocks, a block index and fast-decoding codecs only |
| `--fast`, `--slow` | the search modes above; not both |
| `--block N`, `--block all` | rows per block by hand; `all`, the default, is one block |
| `--level N` | zstd level, 1 to 22, default 17; given, zstd-22 is no longer tried |
| `--sort`, `--sort=MEM` | order records by label tuple, then lengths; spill to `--tmp-dir` past `MEM` |
| `--names keep` | store names even when they are 0,1,2,... |
| `--names drop`, `--names drop:F` | store none, or write them to `F` beside the track |
| `--db features.tsv` | one BED per feature set, or one merged featureID track; every reader takes it, and region then prints one label per set |
| `--hier hierarchy.tsv` | with `--db`: a label tuple outside features.tsv, as a smoothed track carries, gets 2^26 plus its rank in the sorted per-set vocabulary (the features.tsv column, the hierarchy nodes of that set, `novel`), the same id in every file; a track that used one needs the same two files to read, and smooth takes its HIER as the `--hier` |
| `--aln P.paf --ref R` | do not store what the reference already gives; every reader takes the same pair. `R` is the reference BED, or that BED encoded plainly as a `.kst`, which a query reads one contig of instead of parsing the whole text |
| `--aln-names F` | with `--aln`, on encode and decode: ordinal N is line N+1 of `F` |
| `--tail 0\|1`, `--tok-chunk N`, `--id-codec N`, `--len-codec N` | pin what the search would choose: the last-run trick, the name tokeniser's chunk, the label and length codecs |
| `--names F`, `--names ordinal` | on the readers: names from a `drop:F` file, a FASTA, a FASTQ or a list, one a line; or accept a track whose names were dropped |
| `--raw`, `--out DIR` | on decode: what is stored, reversing no transform; one BED per feature set under `--db` |
| `--set A,B` | decode and region under `--db`: which feature sets, default all. decode without `--out` and region print the picked labels as columns, merged where they all agree, so one set gives that set's own track |
| `--max-gap N`, `--v1` | on smooth; `SEQ START END` smooths the whole sequence and prints the overlap, so it equals the full smooth |
| `--threads N`, `-t N` | 1 to 1024. Default `OMP_NUM_THREADS`, else cores capped at 8 |
| `--tmp-dir DIR`, `--profile`, `--help`, `--version` | any command |

Every flag is read before the command runs, so it works on either side of the command name.
Two environment variables are read and neither changes the stored bytes: `TMPDIR` when `--tmp-dir` is not given, and `OMP_NUM_THREADS` when `--threads` is not given, so a SLURM allocation is honoured.

`--sort` reorders records, so a decode of a sorted track returns the rows in sorted order.
Compare such a decode against a sorted copy of the input.

The name index `karyozip index` appends holds a 32-bit fingerprint and a bit-packed ordinal a name, 2,044,644 bytes for 300,000 reads, so `region` finds a read by name without decoding every name.

## Library

An annotator can write a `.kst` directly instead of emitting a BED for karyozip to re-parse.
`make` builds `c/build/libkaryozip.a`; the API is six calls in `c/karyozip.h`: `kt_open`, `kt_block`, `kt_level`, `kt_row`, `kt_close`, `kt_version`.

```c
kt_writer *w = kt_open(f);
kt_row(w, name, nlen, run_len, label, llen, start);   /* once per run */
kt_close(w);
```

Rows must arrive as a track is already laid out: grouped by sequence, each sequence appearing once, starting at 0, gapless.
Those three properties were checked on 264,914,318 real rows across all four data classes with zero violations, and anything else is refused with the message the CLI gives.

## Dependencies

zstd and liblzma, and the vendored htscodecs and bzip3 under `c/`.
Both vendored libraries are pinned: `c/htscodecs/PROVENANCE` and `c/bzip3/PROVENANCE` name the release and every local patch.
The bzip3 patches matter: upstream 1.5.3 loses half the data on an input that is an exact multiple of the block size, and `c/bzip3/lastblock_check.c` reproduces it.
The vendored `libbz3.c` also encodes and decodes its independent blocks in parallel with OpenMP, writing the same bytes as the serial loop.
