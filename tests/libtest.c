/* The library and the CLI write the same file: a BED fed row by row through kt_row, as an annotator would */
#include "karyozip.h"
#include <stdlib.h>
#include <string.h>

/* Contract violations, one per case, each refused with the CLI's message
 */
static int abuse(int c) {
    FILE *f = fopen("/dev/null", "wb");
    kt_writer *w = kt_open(f);
    char big[4096]; memset(big, 'x', sizeof big);
    switch (c) {
    case 1: kt_row(w, "r1", 2, 10, "A", 1, 5); break;             /* not at 0 */
    case 2: kt_row(w, "r1", 2, 10, "A", 1, 0);
            kt_row(w, "r1", 2, 10, "A", 1, 50); break;            /* a gap */
    case 3: kt_row(w, "r1", 2, 10, "A", 1, 0);
            kt_row(w, "r1", 2, 10, "A", 1, 5); break;             /* an overlap */
    case 4: kt_row(w, "r1", 2, 10, "A", 1, 0);
            kt_row(w, "r2", 2, 10, "A", 1, 0);
            kt_row(w, "r1", 2, 10, "A", 1, 10); break;            /* a sequence returns */
    case 5: kt_row(w, "r1", 2, 0, "A", 1, 0); break;              /* zero-length run */
    case 6: kt_row(w, big, sizeof big, 10, "A", 1, 0); break;     /* name too long */
    case 7: kt_row(w, "", 0, 10, "A", 1, 0); break;               /* empty name */
    case 8: kt_row(w, "r1", 2, 10, "A\nB", 3, 0); break;          /* newline label */
    case 9: kt_row(w, "r\0x", 3, 10, "A", 1, 0); break;           /* NUL in the name */
    case 10: kt_row(w, "a\tb", 3, 10, "A", 1, 0); break;          /* tab in the name */
    case 11: kt_row(w, "r1", 2, ~0ULL, "A", 1, 0);
             kt_row(w, "r1", 2, 10, "A", 1, ~0ULL); break;        /* coordinate overflow */
    case 12: kt_level(w, 99); break;                              /* level out of range */
    default: return 2;
    }
    kt_close(w);
    return 0;
}

/* A BED through one writer; rows 0 leaves the block size alone */
static int bedfile(const char *inp, const char *outp, uint64_t rows) {
    FILE *in = fopen(inp, "r"), *out = fopen(outp, "wb");
    if (!in || !out) { perror("open"); return 2; }
    kt_writer *w = kt_open(out);
    if (rows) kt_block(w, rows);
    char *line = NULL; size_t cap = 0; ssize_t n;
    while ((n = getline(&line, &cap, in)) > 0) {
        if (line[n - 1] == '\n') line[--n] = 0;
        char *t1 = strchr(line, '\t'); if (!t1) return 3;
        char *t2 = strchr(t1 + 1, '\t'); if (!t2) return 3;
        char *t3 = strchr(t2 + 1, '\t'); if (!t3) return 3;
        uint64_t st = strtoull(t1 + 1, NULL, 10), en = strtoull(t2 + 1, NULL, 10);
        kt_row(w, line, (size_t)(t1 - line), en - st,
               t3 + 1, strlen(t3 + 1), st);
    }
    kt_close(w);
    free(line); fclose(in);
    if (fclose(out)) { perror("close"); return 2; }
    return 0;
}

int main(int argc, char **argv) {
    if (argc == 3 && !strcmp(argv[1], "-abuse")) return abuse(atoi(argv[2]));
    if (argc == 4 && !strcmp(argv[1], "-twice")) {      /* two writers in one process */
        FILE *a = fopen(argv[2], "wb"), *b = fopen(argv[3], "wb");
        kt_writer *w = kt_open(a); kt_row(w, "s1", 2, 10, "AAA", 3, 0); kt_row(w, "s1", 2, 10, "BBB", 3, 10); kt_close(w); fclose(a);
        w = kt_open(b); kt_row(w, "t1", 2, 7, "CCC", 3, 0); kt_row(w, "t1", 2, 5, "AAA", 3, 7); kt_close(w); fclose(b);
        return 0;
    }
    if (argc == 7 && !strcmp(argv[1], "-seq")) {       /* two BEDs, the second at BLOCK rows a block, in one process */
        int r = bedfile(argv[2], argv[3], 0);
        return r ? r : bedfile(argv[4], argv[5], strtoull(argv[6], NULL, 10));
    }
    if (argc < 3) { fprintf(stderr, "usage: libtest IN.bed OUT.kst [BLOCK]\n"
                                    "       libtest -seq A.bed A.kst B.bed B.kst BLOCK\n"
                                    "       libtest -abuse N\n"); return 2; }
    return bedfile(argv[1], argv[2], argc > 3 ? strtoull(argv[3], NULL, 10) : 0);
}
