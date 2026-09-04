/* karyozip: write a KaryoScope Track (.kst) directly, without a BED in
   between. The pipeline today emits a BED and karyozip re-parses it; on one
   3 Gbp HiFi sample that intermediate is 3.65 GB written and read back, and
   the parse is about 40% of encode.

   Rows must arrive the way a track is already laid out, which was checked on
   264,914,318 real rows across assemblies, HiFi, ONT and Illumina with zero
   violations:
     - runs grouped by sequence, and a sequence never reappears once left
     - every sequence starts at 0
     - the tiling is gapless: each start is the previous end

   Anything else is refused with the same message the CLI gives. One writer
   at a time: the encoder keeps global state.

       kt_writer *w = kt_open(f);
       kt_block(w, 0);                      // optional: 0 = one block
       for each run:  kt_row(w, name, nlen, len, label, llen, start);
       kt_close(w);                         // writes the file, frees w   */
#ifndef KARYOZIP_H
#define KARYOZIP_H

#include <stdio.h>
#include <stddef.h>
#include <stdint.h>

typedef struct ktw kt_writer;

/* out must be open for writing. Nothing is written until kt_close. */
kt_writer *kt_open(FILE *out);

/* Rows per block, before the first row. 0 means one block: smallest file,
   no random access. The default is 262144. */
void kt_block(kt_writer *w, uint64_t rows);

/* zstd level 1 to 22, before the first row. The default is 17. */
void kt_level(kt_writer *w, int level);

/* One run. start is the previous run's end, or 0 for a new sequence. */
void kt_row(kt_writer *w, const char *name, size_t nlen,
            uint64_t run_len, const char *label, size_t llen, uint64_t start);

/* Finishes the file and frees w. Any error exits with a message, as the
   CLI does. */
void kt_close(kt_writer *w);

/* The versions the stored bytes depend on. */
const char *kt_version(void);

#endif
