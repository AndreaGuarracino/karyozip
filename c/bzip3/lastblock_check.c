#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libbz3.h"
int main(void) {
    size_t bs = 65 * 1024;              /* the minimum block size */
    size_t n  = 2 * bs;                 /* an exact multiple: two full blocks */
    unsigned char *in = malloc(n), *cm = malloc(bz3_bound(n)), *back = malloc(n + 64);
    for (size_t i = 0; i < n; i++) in[i] = (unsigned char)(i * 7 + (i >> 9));
    size_t cn = bz3_bound(n);
    int rc = bz3_compress((unsigned int)bs, in, cm, n, &cn);
    printf("  compress rc=%d  in=%zu  out=%zu\n", rc, n, cn);
    if (rc) return 1;
    size_t bn = n + 64;
    int rd = bz3_decompress(cm, back, cn, &bn);
    printf("  decompress rc=%d  got=%zu  want=%zu  bytes %s\n",
           rd, bn, n, (rd == 0 && bn == n && !memcmp(in, back, n)) ? "MATCH" : "LOST");
    return 0;
}
