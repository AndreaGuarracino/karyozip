# karyozip

Reads and writes the **KaryoScope Track (KST)** binary format: lossless, compact storage and O(1) random access for k-mer annotation tracks. Files are `.kst`, and there is only ever one of them: `karyozip index F` appends the index to the track itself.

Compresses a BED-like track (`seq  start  end  label`, sorted, gapless, one block of rows per sequence) and gives it back byte-identical. Build with `make`. Every knob is a flag, and every flag is read before the command runs, so it works on either side of the command name. `karyozip --help` prints the usage. Two environment variables are read, and neither changes the bytes karyozip stores: `TMPDIR`, where temp files go when `--tmp-dir` is not given, and `OMP_NUM_THREADS`, the thread count when `--threads` is not given, so a SLURM allocation is honoured. Nothing else changes what karyozip does.

Two compression modes: none (quick-decoding codecs) and `--slow` (smaller, slower both ways). Names live in one of two places: in the track, or outside it, in which case decode needs `--names FILE`.

## encode

```
karyozip encode [flags] < in.bed > out.kst        # or: producer | karyozip encode > out.kst
```

| flag | effect |
|---|---|
| `--block N` | rows per block, default 262144. Small blocks, fast queries. |
| `--block all` | one block, no random access, smallest file. |
| `--level N` | zstd level, 1 to 22, default 17. |
| `--slow` | every codec over every block: smaller file, slower encode and decode. |
| `--names keep` | store names even when they are 0,1,2,... |
| `--names drop` | store no names. Decode needs `--names FILE`, a FASTA, FASTQ or one name per line. |
| `--names drop:F` | store no names, write them to `F`. Decode needs `--names F` itself, nothing else. |
| `--sort` | order records by label tuple, then lengths, then input rank. |
| `--sort=MEM` | the same, spilling past `MEM` bytes instead of holding every record. `1M` and up; spills under `--tmp-dir`. |
| `--db features.tsv` | inputs are one BED per feature set, or one merged featureID track. |
| `--aln P.paf --ref R` | labels and lengths the reference already gives are not stored. Keep the packed alignment beside the track: it is part of the stored size. |
| `--aln-names F` | for `--aln` on a track that names rows by ordinal: ordinal N is line N+1 of `F`. |
| `--tail 0\|1` | force the tail predictor off or on. |
| `--tok-chunk N` | names per tokeniser chunk. |
| `--id-codec N`, `--len-codec N` | force one codec, for tests. |

## decode

```
karyozip decode out.kst > back.bed
karyozip decode --db features.tsv --out DIR out.kst      # one BED per set
karyozip decode --aln P.paf --ref R out.kst              # same PAF and reference as encode
karyozip decode --names F out.kst                        # drop:F wrote F; a drop track takes a FASTA, FASTQ or list
karyozip decode --aln P.paf --ref R --aln-names F out.kst # ordinal-named rows: F maps ordinal to query name
karyozip decode --names ordinal out.kst                  # names were dropped, print ordinals
karyozip decode --raw out.kst                            # what is stored, no transform reversed
```

## read

```
karyozip region F SEQ                    the whole sequence
karyozip region F SEQ POS                the run covering POS
karyozip region F SEQ START END          rows overlapping the range
karyozip region F SEQ START END HALO      plus HALO runs on each side
karyozip smooth F HIER FS                promote noisy runs, whole track
karyozip smooth F HIER FS SEQ            one sequence
karyozip smooth F HIER FS SEQ START END  one range
karyozip index  F                        appends an index to F, faster lookups
karyozip aln-pack P.paf REF.bed OUT      packs a PAF to only what --aln reads
```

`HIER` is the database `hierarchy.tsv`, three columns, `feature_set  child  parent`, with a header row; `FS` is one `feature_set` in it.

`--names FILE` and `--names ordinal` work on every reader, and `--aln-names F` goes with `--aln`, `smooth` also takes `--max-gap N` (default 1000), `--halo N` (default 200) and `--v1` (drop the `_specific` suffix).

## any command

| flag | effect |
|---|---|
| `--threads N`, `-t N` | threads, 1 to 1024. Default `OMP_NUM_THREADS`, else all cores. Output is byte-identical at every count. |
| `--profile` | phase times and peak memory on stderr. |
| `--tmp-dir DIR` | where the block temp file goes, else `TMPDIR`, else `/tmp`. |

## combinations

Encode flags combine freely except: `--names keep` with `--names drop`, and `--aln` with `--names drop` and no name file.
Decode must be given the same `--db` table, the same PAF and reference, and the same name source as encode; a mismatch is refused, never guessed.
`--sort` reorders records, so a decode of a sorted track returns the same rows in the sorted order, never the input order. Compare such a decode against a sorted copy of the input.
A flag used on the wrong command is refused: `encode --names ordinal` and `decode --tail 1` both fail with a message.
An unknown flag, a flag with no value, and a second track given to `decode` are all refused.

## library

An annotator can write a `.kst` directly instead of emitting a BED for karyozip to re-parse. `make` builds `build/libkaryozip.a` beside the binary; the API is six calls in `karyozip.h`.

```c
#include "karyozip.h"
kt_writer *w = kt_open(f);
kt_block(w, 0);                    /* optional: 0 is one block */
kt_row(w, name, nlen, run_len, label, llen, start);   /* once per run */
kt_close(w);
```

Rows must arrive as a track is already laid out: grouped by sequence, each sequence appearing once, starting at 0, gapless. Those three properties were checked on 264,914,318 real rows across all four data classes with zero violations, and anything else is refused with the message the CLI gives. One writer at a time.

A track built this way is byte-identical to `karyozip encode` on the equivalent BED; `tests/libtest.c` asserts it and the suite runs it.

## tests

```
make test                                    # about 230 checks, LeakSanitizer, five fuzz stages
bash ../scripts/quicktest.sh                 # the above plus five real cells, 80 s
bash ../scripts/realtest.sh BINARY           # one real cell, 21 checks
```
