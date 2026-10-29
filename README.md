# karyozip

Compresses KaryoScope annotation tracks losslessly, and reads them back without decompressing the whole file.

## Install

```
git clone https://github.com/AndreaGuarracino/karyozip && cd karyozip
cargo install --path rust       # puts karyozip in ~/.cargo/bin
```

It needs cargo, a C compiler and libgomp: the codecs under `c/` are compiled into the binary.

## Quick start

```
karyozip compress   -i sample.chromosome.presmoothed.bed.gz -o sample.kst
karyozip decompress -i sample.kst -o sample.bed                           # the original BED, byte for byte
karyozip query      -i sample.kst -r chr1_hap1:1000000-1200000            # the runs overlapping a range
karyozip decompress -i sample.kst --smooth --hier hierarchy.tsv           # KaryoScope's smoothed track, from the presmoothed one
```

Every text input may be plain or gzipped. Without `-i`, `compress` reads stdin, so it can sit behind `karyoscope annotate` in a pipe; without `-o`, decompress and query write stdout.

## Commands

```
karyozip compress   -i IN.bed[.gz] -o OUT.kst [flags]
karyozip decompress -i F                       the whole track as BED
karyozip query      -i F -r NAME[:START-END]   one sequence, or the runs overlapping a range of it
karyozip index      -i F                       writes F.kzi beside the track: query finds a read by name without decoding the names
```

`START-END` is BED: 0-based, END excluded, so `-r chr1:150-151` is the run covering position 150.
A name may contain `:`; only a final `:START-END` of digits is read as a range.
`index` never changes the track; `query` uses `F.kzi` when it is there and matches the track, and never builds it.
The index costs 7.5 B a read: on short reads it is larger than the track (223334174 B track, 481554594 B index), on long reads and assemblies 0 to 2% of it.

## Size or speed

| you want | use |
|---|---|
| a small file and a `query` that decodes one block, the usual choice | `compress`, or `-e medium`, then `index` for many reads |
| the smallest file; a `query` decodes the whole file | `compress -e high` |
| the quickest compress and `query`, a larger file | `compress -e low`, then `index` for many reads |

## Smoothing, feature sets and alignments

A KaryoScope run writes one BED per feature set, presmoothed and smoothed.

- **Smoothed tracks need not be stored.** `--smooth --hier hierarchy.tsv` on decompress or query rebuilds the smoothed track from the presmoothed one, row for row as KaryoScope writes it. The set is the one whose tree in hierarchy.tsv holds every label of the track; when two do, or for a `-d` track, `--set SET` names it. Tracks written by KaryoScope 1.x, whose smoothed labels carry no `_specific` suffix, need `--v1`; those its first release (2025-12-10) smoothed also need `--chunk 20000`, its chunks of 20000 runs per sequence.
- **All feature sets in one file.** `-d features.tsv` merges one sample's per-set BEDs into a single, smaller track. Smoothed sets merged this way also need `--hier hierarchy.tsv`, on compress and on every read.

```
karyozip compress   -d features.tsv -i s.acrocentric.bed.gz -i s.chromosome.bed.gz ... -o s.all.kst
karyozip decompress -d features.tsv -i s.all.kst --out-dir DIR      # every set's BED back
karyozip decompress -d features.tsv -i s.all.kst --set gene         # one set
karyozip decompress -d features.tsv -i s.all.kst --smooth --hier hierarchy.tsv --out-dir DIR   # every set smoothed, one BED each
```

- **Against a reference.** With the reads aligned to a reference (a PAF) and the reference's own track, `-a P.paf --ref ref.bed` stores only what differs from the reference. Every read of such a file takes the same `-a` and `--ref`. Give `--ref` as a `.kst` track (`compress -e low`) for queries: a point query on the whole CHM13 reference takes 0.10 s against 13.0 s from the BED.

## Read names

Read names are stored in the file by default.

| flag | effect |
|---|---|
| `-n drop` | store none; the file is smaller, and every read of it needs `-n F` |
| `-n F` | on decompress and query: take names from `F`, a FASTA, a FASTQ or a list, one a line |

## All flags

| flag | on | effect |
|---|---|---|
| `-i`, `--input FILE` | all but aln-pack | compress: the BED, one `-i` per set under `-d`, else stdin; the others: the track |
| `-o`, `--output FILE` | all but index | where the output goes, else stdout; aln-pack: the packed file, required |
| `-r`, `--region NAME[:START-END]` | query | a sequence, or a range of it |
| `--smooth` | decompress, query | KaryoScope's smoothed track instead of the stored one; with query the whole sequence is smoothed and the range printed |
| `--hier FILE` | `--smooth`; compress and readers with `-d` | the KaryoScope hierarchy.tsv |
| `--v1`, `--chunk N` | `--smooth` | KaryoScope 1.x output; its chunked smoothing before 2026-01-08 |
| `-e`, `--effort low\|medium\|high` | compress | see Size or speed; default medium |
| `-n`, `--names ...` | compress, readers | see Read names |
| `-d`, `--db features.tsv` | compress, readers | see Smoothing, feature sets and alignments |
| `--set A,B` | decompress, query | the sets of a `-d` track to print, default all; with `--smooth`, the one set smoothed |
| `--out-dir DIR` | decompress | one BED per set under `-d` |
| `-a`, `--aln P.paf`, `--ref R` | compress, readers, aln-pack | see Smoothing, feature sets and alignments |
| `-t`, `--threads N` | all | default the cores, up to 8 |
| `--tmp-dir DIR` | all | where temporary files go, default /tmp |
| `-h`, `--help`, `--version` | all | usage; library versions |

Flags work on either side of the command name. The readers are decompress and query.

## Advanced

| command or flag | on | effect |
|---|---|---|
| `aln-pack -a P.paf --ref R -o OUT` | | cuts a PAF to what `-a` reads, for many tracks against one alignment |
| `--aln-names F` | compress, readers | with `-a`: read N is line N+1 of `F` |
| `--max-gap N` | `--smooth` | KaryoScope's gap bound, 1000 by default |

## How it works

**How compress chooses.** The three efforts are one rule with two settings: `low` tries 8 of the 51 codecs plus zstd at 17 on 4 MB, `medium` all 51 plus zstd at 17 and 22 on 16 MB, `high` the same on the whole file.
The smallest wins, with no margin, and must read back what it was ranked on.
A one-block file also chooses its layout on the same sample: lengths and labels apart, run counts heading each read's lengths, lengths grouped by label, or labels interleaved with lengths; every layout and codec is ranked in one parallel pool.
Names are ranked the same way among zstd, integer deltas, field deltas and the CRAM name tokeniser.
Compress prints what it chose on stderr, on its `layout:` and `names:` lines.
`-e low` also writes 262144-row blocks and a block index, which a query decodes one at a time, and so drops the codecs slow to decode, bzip3, arith and LZMA, on a track of more than one block.
`-e medium` writes blocks of whole reads, each closed at the first read that starts past 4194304 rows; every block takes the layout the first block chose, so a query decodes one block. `-e high` writes one block.
Decompress decodes the next blocks on the other threads while it writes the rows of one, when the blocks hold 65536 rows or more on average.

**Name index.** `karyozip index` writes to `F.kzi` a 32-bit fingerprint and a bit-packed ordinal per name, 2,044,644 B for 300,000 reads; one fingerprint match is answered without decoding the names, so a name absent from the track passes for a present one about 4 times in 2^32 lookups; two or more matches are confirmed against the decoded names.

## C and Rust

A C implementation in `c/` writes the same bytes. How to build it, its library, and how the two are kept equal: [C-and-Rust.md](C-and-Rust.md).
