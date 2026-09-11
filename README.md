# karyozip

Compresses KaryoScope annotation tracks losslessly, and reads them back without decompressing the whole file.

## Install

```
make            # builds c/build/karyozip and rust/target/release/karyozip
```

The C needs zstd, liblzma, zlib and an OpenMP compiler; the Rust needs cargo, liblzma and libgomp.

## Quick start

```
karyozip encode -i sample.chromosome.presmoothed.bed.gz -o sample.kst
karyozip decode -i sample.kst -o sample.bed                               # the original BED, byte for byte
karyozip region -i sample.kst -s chr1_hap1 --start 1000000 --end 1200000  # the runs overlapping a range
karyozip smooth -i sample.kst --hier hierarchy.tsv -f chromosome          # KaryoScope's smoothing, from the presmoothed track
```

Every text input may be plain or gzipped. Without `-i`, `encode` reads stdin, so it can sit behind `karyoscope annotate` in a pipe; without `-o`, every command writes stdout.


## Commands

```
karyozip encode -i IN.bed[.gz] -o OUT.kst [flags]
karyozip decode -i F                                  the whole track as BED
karyozip region -i F -s SEQ                           one sequence
karyozip region -i F -s SEQ -p POS                    the run covering POS
karyozip region -i F -s SEQ --start N --end N [--halo N]   the runs overlapping a range, plus N runs each side
karyozip smooth -i F --hier HIER -f FS [-s SEQ [--start N --end N]]   KaryoScope's smoothing of feature set FS
karyozip index  -i F                                  adds a name index, so region finds a read fast on tracks with many reads
karyozip aln-pack -a P.paf -r REF.bed|REF.kst -o OUT  cuts a PAF to what -a needs
```

## Size or speed

| you want | use |
|---|---|
| a small file, the usual choice | `encode` |
| the smallest file, however long it takes | `encode --slow` |
| a quick encode | `encode --fast` |
| fast `region` queries on a big track | `encode --queryable`, with any of the above |

## Feature sets, smoothing and alignments

A KaryoScope run writes one BED per feature set, presmoothed and smoothed.

- **Smoothed tracks need not be stored.** `karyozip smooth` rebuilds a smoothed track from the presmoothed one, row for row as KaryoScope writes it. The exception is the v1 `region` and `acrocentric` sets, whose hierarchy changed after they were written: store those.
- **All feature sets in one file.** `-d features.tsv` merges one sample's per-set BEDs into a single, smaller track. Smoothed sets merged this way also need `--hier hierarchy.tsv`, on encode and on every read.

```
karyozip encode -d features.tsv -i s.acrocentric.bed.gz -i s.chromosome.bed.gz ... -o s.all.kst
karyozip decode -d features.tsv -i s.all.kst --out-dir DIR      # every set's BED back
karyozip decode -d features.tsv -i s.all.kst --set gene         # one set
```

- **Against a reference.** With the reads aligned to a reference (a PAF) and the reference's own track, `-a P.paf -r ref.bed` stores only what differs from the reference. Every read of such a file takes the same `-a` and `-r`.

## Read names

Read names are stored in the file by default.

| flag | effect |
|---|---|
| `-n drop` | store none; the file is smaller and a read is known by its number |
| `-n drop:F` | write the names to `F` beside the track |
| `-n F` | on decode, region and smooth: take names from `F`, a `drop:F` file, a FASTA, a FASTQ or a list, one a line |
| `-n ordinal` | on decode: accept a track stored with `-n drop` and print the numbers |
| `-n keep` | store names even when they are 0, 1, 2, ... |

`--sort` groups records by their labels, which pays on short reads stored with `-n drop`; a decode then returns the rows in that order.

## All flags

| flag | on | effect |
|---|---|---|
| `-i`, `--input FILE` | all | encode: the BED, one `-i` per set under `-d`, else stdin; the others: the track |
| `-o`, `--output FILE` | all but index | where the output goes, else stdout; aln-pack: the packed file |
| `-s`, `--seq SEQ` | region, smooth | the sequence |
| `-p`, `--pos N` | region | the run covering position N |
| `--start N`, `--end N` | region, smooth | a range; region prints the runs overlapping it, smooth smooths the whole sequence and prints the overlap |
| `--halo N` | region | N more runs each side of the range |
| `-f`, `--feature-set FS` | smooth | the feature set to smooth |
| `--hier FILE` | smooth; encode and readers with `-d` | the KaryoScope hierarchy.tsv |
| `--fast`, `--slow`, `--queryable` | encode | see Size or speed |
| `-b`, `--block N`, `--block all` | encode | rows per block by hand; `all`, the default, is one block |
| `-l`, `--level N` | encode | zstd level, 1 to 22, default 17 |
| `--sort`, `--sort=MEM` | encode | group records by labels; spill to disk past `MEM` |
| `-n`, `--names ...` | encode, readers | see Read names |
| `-d`, `--db features.tsv`, `--set A,B` | encode, readers | see Feature sets |
| `-a`, `--aln P.paf`, `-r`, `--ref R`, `--aln-names F` | encode, readers, aln-pack | see Feature sets; `--aln-names`: read N is line N+1 of `F` |
| `--out-dir DIR`, `--raw` | decode | one BED per set under `-d`; `--raw` prints what is stored, for debugging |
| `--max-gap N`, `--v1` | smooth | KaryoScope's smoothing options |
| `-t`, `--threads N` | all | default `OMP_NUM_THREADS`, else the cores, up to 8 |
| `--tmp-dir DIR` | all | where temporary files go, default `TMPDIR`, else /tmp |
| `--profile`, `-h`, `--help`, `--version` | all | time per phase; usage; library versions |
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
