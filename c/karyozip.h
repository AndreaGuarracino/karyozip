/* karyozip: write a KaryoScope Track (.kst) directly, without a BED in between */
/* rows grouped by sequence, a sequence never back once left, each from 0, gapless; else refused with the CLI's message */
/* w = kt_open(f); kt_block, kt_level optional; kt_row(w, name, nlen, len, label, llen, start) a run; kt_close(w) writes and frees */
/* Not thread-safe: the writer keeps process-wide state, so one writer at a time per process, from one thread; kt_open starts the next
   from nothing. Every error, a refused row included, prints "karyozip: ..." to stderr and calls exit(1); nothing returns an error */
#ifndef KARYOZIP_H
#define KARYOZIP_H

#include <stdio.h>
#include <stddef.h>
#include <stdint.h>

typedef struct ktw kt_writer;

/* out must be open for writing; nothing is written until kt_close */
kt_writer *kt_open(FILE *out);

/* Rows per block; 0 is one block, no random access; default one block (the CLI's -e low sets 262144).
   Process-wide, as kt_level: call after kt_open and before the first kt_row; w is not consulted */
void kt_block(kt_writer *w, uint64_t rows);

/* zstd level 1 to 22, after kt_open and before the first kt_row; default 17 */
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
