/* karyozip: write a KaryoScope Track (.kst) directly, without a BED in between */
/* rows grouped by sequence, a sequence never back once left, each from 0, gapless; else refused with the CLI's message; one writer at a time */
/* w = kt_open(f); kt_block(w, 0) optional, 0 one block; kt_row(w, name, nlen, len, label, llen, start) a run; kt_close(w) writes and frees */
#ifndef KARYOZIP_H
#define KARYOZIP_H

#include <stdio.h>
#include <stddef.h>
#include <stdint.h>

typedef struct ktw kt_writer;

/* out must be open for writing; nothing is written until kt_close */
kt_writer *kt_open(FILE *out);

/* Rows per block, before the first row; 0 is one block, no random access; default 262144 */
void kt_block(kt_writer *w, uint64_t rows);

/* zstd level 1 to 22, before the first row; default 17 */
void kt_level(kt_writer *w, int level);

/* One run: start is the previous run's end, or 0 for a new sequence */
void kt_row(kt_writer *w, const char *name, size_t nlen,
            uint64_t run_len, const char *label, size_t llen, uint64_t start);

/* Finishes the file and frees w; an error exits with a message, as the CLI
 */
void kt_close(kt_writer *w);

/* The versions the stored bytes depend on */
const char *kt_version(void);

#endif
