# One target, one binary. The tokeniser and the BWT codec are not optional
# and no longer sit behind a define. Everything vendored lives beside this file.
CC      ?= gcc
CFLAGS  ?= -O3 -std=gnu11 -Wall -Wextra -fopenmp
# VERSION is bzip3's, which libbz3.c needs; the KST format version is in
# the magic. They are not the same number and never were.
BZIP3_VERSION ?= 1.5.4
CFLAGS  += -DVERSION='"$(BZIP3_VERSION)"' -I.
ZSTD_CFLAGS := $(shell pkg-config --cflags libzstd 2>/dev/null)
ZSTD_LIBS   := $(shell pkg-config --libs libzstd 2>/dev/null || echo -lzstd)
ifneq ($(CONDA_PREFIX),)
  CFLAGS   += -I$(CONDA_PREFIX)/include
  LDFLAGS  += -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib
endif

SRC := karyozip.c \
       htscodecs/tokenise_name3.c htscodecs/rANS_static4x16pr.c \
       htscodecs/rANS_static32x16pr.c htscodecs/arith_dynamic.c \
       htscodecs/rle.c htscodecs/pack.c htscodecs/utils.c \
       bzip3/libbz3.c \
       libdeflate/lib/deflate_decompress.c libdeflate/lib/crc32.c libdeflate/lib/utils.c \
       libdeflate/lib/x86/cpu_features.c libdeflate/lib/arm/cpu_features.c

# liblzma 5.8.4, raw LZMA2 only, own flags: its config.h, common.h and
# version.h share names with htscodecs' and bzip3's
XZ_SRC := xz/kzlzma.c $(addprefix xz/src/liblzma/, common/common.c \
       common/filter_buffer_decoder.c common/filter_buffer_encoder.c common/filter_common.c \
       common/filter_decoder.c common/filter_encoder.c check/crc32_fast.c \
       lz/lz_decoder.c lz/lz_encoder.c lz/lz_encoder_mf.c lzma/fastpos_table.c \
       lzma/lzma2_decoder.c lzma/lzma2_encoder.c lzma/lzma_decoder.c lzma/lzma_encoder.c \
       lzma/lzma_encoder_optimum_fast.c lzma/lzma_encoder_optimum_normal.c \
       lzma/lzma_encoder_presets.c rangecoder/price_table.c)
XZ_FLAGS := -DHAVE_CONFIG_H -Ixz -Ixz/src/common \
       $(addprefix -Ixz/src/liblzma/, api common check lz lzma rangecoder delta simple)
XZ_HDR := $(shell find xz -name '*.h')
XZ_OBJ := $(XZ_SRC:%.c=build/xz/%.o)
build/xz/%.o: %.c $(XZ_HDR)
	mkdir -p $(@D)
	$(CC) -O3 -std=gnu11 -w $(XZ_FLAGS) -c $< -o $@
XZ_SAN := $(XZ_SRC:%.c=build/xz-san/%.o)
build/xz-san/%.o: %.c $(XZ_HDR)
	mkdir -p $(@D)
	$(CC) -O1 -g -fsanitize=address,undefined -std=gnu11 -w $(XZ_FLAGS) -c $< -o $@

# libtest is in all on purpose: a stale one silently compares the wrong
# binary and every library check "differs".
all: build/karyozip build/libkaryozip.a build/libtest

HDR := karyozip.h $(wildcard htscodecs/*.h bzip3/*.h libdeflate/*.h libdeflate/lib/*.h libdeflate/lib/*/*.h)
build/karyozip: $(SRC) $(HDR) $(XZ_OBJ)
	mkdir -p build
	$(CC) $(CFLAGS) $(ZSTD_CFLAGS) -o $@ $(SRC) $(XZ_OBJ) $(LDFLAGS) $(ZSTD_LIBS) -lm -lpthread -lz

# the same sources, minus main: an annotator can write .kst without a BED.
# Objects keep their paths: two vendored trees both have a utils.c.
OBJ := $(SRC:%.c=build/lib/%.o)
build/lib/%.o: %.c $(HDR)
	mkdir -p $(@D)
	$(CC) $(CFLAGS) $(ZSTD_CFLAGS) -Wno-unused-function -Wno-unused-variable -DKARYOZIP_LIB -c $< -o $@
build/libkaryozip.a: $(OBJ) $(XZ_OBJ)
	rm -f $@ && $(AR) rcs $@ $(OBJ) $(XZ_OBJ)

build/libtest: build/libkaryozip.a ../tests/libtest.c
	$(CC) $(CFLAGS) $(ZSTD_CFLAGS) -o $@ ../tests/libtest.c build/libkaryozip.a $(LDFLAGS) $(ZSTD_LIBS) -lm -lpthread -lz

# the memory-checked twin fuzz.sh runs: one build, kept until a
# source changes, where each script used to compile its own every run
build/karyozip-san: $(SRC) $(HDR) $(XZ_SAN)
	mkdir -p build
	$(CC) -O1 -g -fsanitize=address,undefined -std=gnu11 -fopenmp -DVERSION='"$(BZIP3_VERSION)"' -I. \
	  $(ZSTD_CFLAGS) -o $@ $(SRC) $(XZ_SAN) $(LDFLAGS) $(ZSTD_LIBS) -lm -lpthread -lz

test: build/karyozip build/libtest build/karyozip-san
	bash ../tests/run.sh build/karyozip build/libtest build/karyozip-san

clean:
	rm -rf build

.PHONY: all test clean
