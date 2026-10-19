/* LZMA2 raw buffers with every setting pinned, shared by the C and the Rust */
#include <stddef.h>
#include <stdint.h>
size_t kz_lzma_encode(const uint8_t *src, size_t n, uint8_t *out, size_t ob);
int kz_lzma_decode(const uint8_t *src, size_t n, uint8_t *dst, size_t raw);
const char *kz_lzma_version(void);
