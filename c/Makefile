# One target, one binary. The tokeniser and the BWT codec are not optional
# and no longer sit behind a define. Everything vendored lives beside this file.
CC      ?= gcc
CFLAGS  ?= -O3 -std=gnu11 -Wall -Wextra -fopenmp
# VERSION is bzip3's, which libbz3.c needs; the KST format version is in
# the magic. They are not the same number and never were.
BZIP3_VERSION ?= 1.5.3
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

# libtest is in all on purpose: a stale one silently compares the wrong
# binary and every library check "differs".
all: build/karyozip build/libkaryozip.a build/libtest

HDR := karyozip.h $(wildcard htscodecs/*.h bzip3/*.h libdeflate/*.h libdeflate/lib/*.h libdeflate/lib/*/*.h)
build/karyozip: $(SRC) $(HDR)
	mkdir -p build
	$(CC) $(CFLAGS) $(ZSTD_CFLAGS) -o $@ $(SRC) $(LDFLAGS) $(ZSTD_LIBS) -llzma -lm -lpthread -lz

# the same sources, minus main: an annotator can write .kst without a BED.
# Objects keep their paths: two vendored trees both have a utils.c.
OBJ := $(SRC:%.c=build/lib/%.o)
build/lib/%.o: %.c $(HDR)
	mkdir -p $(@D)
	$(CC) $(CFLAGS) $(ZSTD_CFLAGS) -Wno-unused-function -Wno-unused-variable -DKARYOZIP_LIB -c $< -o $@
build/libkaryozip.a: $(OBJ)
	rm -f $@ && $(AR) rcs $@ $(OBJ)

build/libtest: build/libkaryozip.a ../tests/libtest.c
	$(CC) $(CFLAGS) $(ZSTD_CFLAGS) -o $@ ../tests/libtest.c build/libkaryozip.a $(LDFLAGS) $(ZSTD_LIBS) -llzma -lm -lpthread -lz

# the memory-checked twin leaks.sh and fuzz.sh run: one build, kept until a
# source changes, where each script used to compile its own every run
build/karyozip-san: $(SRC) $(HDR)
	mkdir -p build
	$(CC) -O1 -g -fsanitize=address,undefined -std=gnu11 -fopenmp -DVERSION='"$(BZIP3_VERSION)"' -I. \
	  $(ZSTD_CFLAGS) -o $@ $(SRC) $(LDFLAGS) $(ZSTD_LIBS) -llzma -lm -lpthread -lz

test: build/karyozip build/libtest build/karyozip-san
	KT=$(CURDIR)/build/karyozip KARYOZIP_SAN=$(CURDIR)/build/karyozip-san bash ../tests/run.sh

clean:
	rm -rf build

.PHONY: all test clean
