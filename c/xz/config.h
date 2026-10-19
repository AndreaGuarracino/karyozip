/* liblzma 5.8.4, raw LZMA2 only, for karyozip: hand-written, not configure's */
#define HAVE_ENCODERS 1
#define HAVE_DECODERS 1
#define HAVE_ENCODER_LZMA1 1
#define HAVE_DECODER_LZMA1 1
#define HAVE_ENCODER_LZMA2 1
#define HAVE_DECODER_LZMA2 1
#define HAVE_MF_BT4 1
#define HAVE_CHECK_CRC32 1
#define HAVE_STDBOOL_H 1
#define HAVE__BOOL 1
#define HAVE_STDINT_H 1
#define HAVE_INTTYPES_H 1
#define HAVE_UINTPTR_T 1
#define HAVE_VISIBILITY 0
#define HAVE___BUILTIN_BSWAPXX 1
#define HAVE___BUILTIN_ASSUME_ALIGNED 1
#define ASSUME_RAM 128
#define PACKAGE_NAME "XZ Utils"
#define PACKAGE_URL "https://tukaani.org/xz/"
#if defined(__SIZEOF_SIZE_T__)
#  define SIZEOF_SIZE_T __SIZEOF_SIZE_T__
#endif
#if defined(__x86_64__) || defined(__i386__) || defined(__aarch64__)
#  define TUKLIB_FAST_UNALIGNED_ACCESS 1
#endif
#if defined(__SSE2__)
#  define HAVE_IMMINTRIN_H 1
#  define HAVE__MM_MOVEMASK_EPI8 1
#endif
