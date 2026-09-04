/* Proves the library and the CLI produce the same file. Reads a BED and
   feeds it row by row through kt_row, which is what an annotator would do
   instead of writing a BED at all. */
#include "karyozip.h"
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: libtest IN.bed OUT.kst [BLOCK]\n"); return 2; }
    FILE *in = fopen(argv[1], "r"), *out = fopen(argv[2], "wb");
    if (!in || !out) { perror("open"); return 2; }
    kt_writer *w = kt_open(out);
    if (argc > 3) kt_block(w, strtoull(argv[3], NULL, 10));
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
