# karyozip

Compresses KaryoScope annotation tracks losslessly, and reads them back without decompressing the whole file.

[KaryoScope](https://github.com/barthel-lab/KaryoScope) annotates assemblies and sequencing reads without alignment. It gives every base of every read or contig a label drawn from one or more hierarchical feature sets (chromosome of origin, satellite, repeat class, and so on), and writes one BED track per feature set. Those tracks are what karyozip stores.

## Install

```
git clone https://github.com/AndreaGuarracino/karyozip && cd karyozip
cargo install --path rust       # puts karyozip in ~/.cargo/bin
```

It needs cargo, a C and a C++ compiler and libgomp: the codecs under `c/` are compiled into the binary.
It builds on Linux with glibc and GCC (or a compiler with `-fopenmp` and libgomp), and from a clone only: `rust/build.rs` reads `c/`, so the crate is not published.

## Input

A track is a 4-column BED: sequence name, start, end, label.

```
m84039_230928_213653_s3/1000/ccs	0	18701	chr1
m84039_230928_213653_s3/1000/ccs	18701	27109	chr13
m84039_230928_213653_s3/1007/ccs	0	14778	chrX
```

Every row must meet these rules, and compress refuses the input with a message otherwise:
- Each sequence's rows come together, and a sequence never comes back once left.
- Each sequence starts at 0 and is gapless: every row starts where the one before ended.
- Every run has a positive length.
- Names are unique, under 512 bytes, and contain no tab, newline or NUL. Labels contain no newline or NUL.
- Under `-d`, two touching runs of one sequence must not carry the same label (featureID): merge them first.

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
The index costs 7.5 B a read: on short reads it can be larger than the track, on long reads and assemblies it is 2% of it or less.

## Compression modes

| mode | input | flags, on compress and on every read |
|---|---|---|
| standard, the default | one BED, or stdin | none |
| db-informed | one BED or plain `.kst` per feature set of a sample, in `features.tsv` column order | `-d features.tsv` |
| alignment-informed | one BED, with its reads aligned to a reference | `-a P.paf --ref ref.bed` |
| db + alignment | as db-informed | `-d features.tsv -a P.paf --ref ref.ids.bed` |

Use db-informed when you have every feature set of a sample: one file replaces the per-set files and is smaller than their sum.
Each set may come as its BED or as its plain `.kst` (what `compress` wrote from that BED, in any mode), known by its magic `KST00001` and a first line no BED row could be (a BED whose first name begins with `KST00001` is still a BED), and the two may be mixed: the `.kst` inputs are decoded in step, one row at a time, and the merged file is byte-identical to the one their BEDs give.
A `.kst` input must be plain: not written with `-d`, `-a` or `-n drop`. A missing set, a wrong order or a refused input is refused in the same words as for BEDs.
Standard stays the default because db-informed cannot take one BED or stdin, and its file reads back only with the same `features.tsv`.

An alignment-informed file stores only what differs from the reference. The PAF and the reference stay outside the file, and on short reads they cost more than one file saves.
`ref.bed` is the reference's track of the same feature set. `ref.ids.bed` is the reference's db-informed track, written by `karyozip decompress --raw -i ref.all.kst -o ref.ids.bed`.
For queries, give `--ref` as a plain `.kst` track (`compress -e low`): it is read one contig at a time, not whole.

```
karyozip compress   -d features.tsv -i s.acrocentric.bed.gz -i s.chromosome.bed.gz ... -o s.all.kst
karyozip compress   -d features.tsv -i s.acrocentric.kst -i s.chromosome.kst ... -o s.all.kst   # the same file, from per-set tracks
karyozip decompress -d features.tsv -i s.all.kst --out-dir DIR      # every set's BED back
karyozip decompress -d features.tsv -i s.all.kst --set gene         # one set
karyozip decompress -d features.tsv -i s.all.kst --smooth --hier hierarchy.tsv --out-dir DIR   # every set smoothed, one BED each
```

Smoothed sets under `-d` also need `--hier hierarchy.tsv`, on compress and on every read.

## Size or speed

| you want | use |
|---|---|
| a small file and a `query` that decodes one block, the usual choice | `compress`, or `-e medium`, then `index` for many reads |
| the smallest file; a `query` decodes the whole file | `compress -e high` |
| the quickest compress and `query`, a larger file | `compress -e low`, then `index` for many reads |

## Smoothed tracks

Smoothed tracks need not be stored. `--smooth --hier hierarchy.tsv` on decompress or query rebuilds the smoothed track from the presmoothed one, row for row as KaryoScope writes it. The set is the one whose tree in hierarchy.tsv holds every label of the track; when two do, or for a `-d` track, `--set SET` names it. Tracks of KaryoScope 1.x need `--v1`, and those smoothed before 2026-01-08 also `--chunk 20000`.

## Read names

Read names are stored in the file by default.

| flag | effect |
|---|---|
| `-n drop` | store none; the file is smaller, and every read of it needs `-n F` |
| `-n F` | on decompress and query: take names from `F`, a FASTA, a FASTQ or a list, one a line |

## All flags

| flag | on | effect |
|---|---|---|
| `-i`, `--input FILE` | all but aln-pack | compress: the BED, else stdin; under `-d` one `-i` per set, a BED or a plain `.kst`; the others: the track |
| `-o`, `--output FILE` | all but index | where the output goes, else stdout; aln-pack: the packed file, required |
| `-r`, `--region NAME[:START-END]` | query | a sequence, or a range of it |
| `--smooth` | decompress, query | KaryoScope's smoothed track instead of the stored one; with query the whole sequence is smoothed and the range printed |
| `--hier FILE` | `--smooth`; compress and readers with `-d` | the KaryoScope hierarchy.tsv |
| `--v1`, `--chunk N` | `--smooth` | KaryoScope 1.x output; its chunked smoothing before 2026-01-08 |
| `-e`, `--effort low\|medium\|high` | compress | see Size or speed; default medium |
| `-n`, `--names ...` | compress, readers | see Read names |
| `-d`, `--db features.tsv` | compress, readers | see Compression modes |
| `--set A,B` | decompress, query | the sets of a `-d` track to print, default all; with `--smooth`, the one set smoothed |
| `--out-dir DIR` | decompress | one BED per set under `-d`; the BEDs are replaced all or nothing on a SIGINT, SIGTERM, SIGHUP or failed rename (every new, or every old or absent as before), though an uncatchable kill (SIGKILL, OOM, power loss) mid-commit may leave the old ones as SET.bed.old.XXXXXX beside them |
| `-a`, `--aln P.paf`, `--ref R` | compress, readers, aln-pack | see Compression modes |
| `-t`, `--threads N` | all | default the cores, up to 8 |
| `--tmp-dir DIR` | all | where temporary files go, default /tmp |
| `-h`, `--help`, `--version` | all | usage; library versions |

Flags work on either side of the command name. The readers are decompress and query.
No environment variable changes the bytes karyozip writes. Four are read, for the test suite only: `KZ_TEST_SIG_PAUSE_MS`, `KZ_TEST_SWAP_PAUSE_MS`, `KZ_TEST_STAGE_PAUSE_MS` and `KZ_TEST_STAGE_ENTRY_PAUSE_MS` add sleeps at fixed points of a signal, an output commit or `--out-dir`'s staging, and do nothing when unset ([C-and-Rust.md](C-and-Rust.md#tests)).

## Advanced

| command or flag | on | effect |
|---|---|---|
| `aln-pack -a P.paf --ref R -o OUT` | | cuts a PAF to what `-a` reads, for many tracks against one alignment |
| `--aln-names F` | compress, readers | with `-a`: read N is line N+1 of `F` |
| `--max-gap N` | `--smooth` | KaryoScope's gap bound, 1000 by default |
| `--raw` | decompress | prints what is stored, without the companion files: featureIDs of a `-d` track, `=` labels and 0 lengths of an `-a` track, ordinals for dropped names |
| `-b`, `--block N\|all` | compress | blocks of exactly N rows (a block may split a sequence); `all` is one block. Overrides the block size `-e` picks |
| `--read-blocks N` | compress | blocks of whole sequences, each closing at the first sequence start at or past N rows. The later of `-b` and `--read-blocks` wins |
| `--profile` | all | a debug flag: per-phase times and peak memory on stderr |

## How it works

Compress tries its codecs on a sample, and under `medium` and `high` four layouts of lengths and labels too, keeps the smallest that reads back, and prints its choice on stderr (`layout:` and `names:` lines): `low` tries the fast codecs on 4 MB, `medium` all of them on 16 MB, `high` all of them on the whole file, plus bzip3 at many block sizes and zpaq.
`-e low` writes blocks of 262144 rows with a block index, `-e medium` blocks of whole reads from 4194304 rows, `-e high` one block; a query decodes one block, so under `high` the whole file.
`index` stores a 32-bit fingerprint and an ordinal per name: one matching fingerprint is answered without decoding the names, so a name absent from the track passes for a present one about 4 times in 2^32 lookups.

## C and Rust

A C implementation in `c/` writes the same bytes. How to build it, its library, and how the two are kept equal: [C-and-Rust.md](C-and-Rust.md).
The Rust library's writer (`karyozip::writer::Writer`) returns every error as a `Result` with the CLI's words and never exits the process, and writers may run at once on several threads; the Rust reader, and the C library `libkaryozip.a`, still print a message and exit on an error, and the C library is not thread-safe.

## File format

[FORMAT.md](FORMAT.md) specifies the `.kst` track, the `.kzi` index and the `aln-pack` file. The magic `KST00001` has not changed while flags and codecs were added. 0.3.0 is the first tagged release: what it writes, it reads, and its readers refuse a header flag they do not know by name. A file from an untagged commit before 0.3.0 is not guaranteed to read.
