#include "src/liblzma/api/lzma.h"
#include <string.h>
#include "kzlzma.h"

/* every field pinned: a preset that moved between liblzma releases would
   change the bytes */
static void filt(lzma_filter *f, lzma_options_lzma *o) {
    memset(o, 0, sizeof *o);
    o->dict_size = 64u << 20; o->lc = 3; o->lp = 0; o->pb = 0;
    o->mode = LZMA_MODE_NORMAL; o->nice_len = 273;
    o->mf = LZMA_MF_BT4; o->depth = 0;
    f[0].id = LZMA_FILTER_LZMA2; f[0].options = o;
    f[1].id = LZMA_VLI_UNKNOWN;  f[1].options = NULL;
}

/* bytes written, 0 when out is too small */
size_t kz_lzma_encode(const uint8_t *src, size_t n, uint8_t *out, size_t ob) {
    lzma_filter f[2]; lzma_options_lzma o; filt(f, &o);
    size_t op = 0;
    return lzma_raw_buffer_encode(f, NULL, src, n, out, &op, ob) == LZMA_OK ? op : 0;
}

const char *kz_lzma_version(void) { return lzma_version_string(); }

/* 1 when exactly raw bytes came back */
int kz_lzma_decode(const uint8_t *src, size_t n, uint8_t *dst, size_t raw) {
    lzma_filter f[2]; lzma_options_lzma o; filt(f, &o);
    size_t ip = 0, op = 0;
    return lzma_raw_buffer_decode(f, NULL, src, &ip, n, dst, &op, raw) == LZMA_OK && op == raw;
}
