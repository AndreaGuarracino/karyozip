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
       bzip3/libbz3.c

# libtest is in all on purpose: a stale one silently compares the wrong
# binary and every library check "differs".
all: build/karyozip build/libkaryozip.a build/libtest

build/karyozip: $(SRC)
	mkdir -p build
	$(CC) $(CFLAGS) $(ZSTD_CFLAGS) -o $@ $(SRC) $(LDFLAGS) $(ZSTD_LIBS) -lm -lpthread

# the same sources, minus main: an annotator can write .kst without a BED
build/libkaryozip.a: $(SRC) karyozip.h
	mkdir -p build/lib
	$(CC) $(CFLAGS) $(ZSTD_CFLAGS) -DKARYOZIP_LIB -c $(SRC) && mv *.o build/lib/
	$(AR) rcs $@ build/lib/*.o

build/libtest: build/libkaryozip.a tests/libtest.c
	$(CC) $(CFLAGS) $(ZSTD_CFLAGS) -o $@ tests/libtest.c build/libkaryozip.a $(LDFLAGS) $(ZSTD_LIBS) -lm -lpthread

test: build/karyozip build/libtest
	KT=$(CURDIR)/build/karyozip tests/run.sh

clean:
	rm -rf build

.PHONY: all test clean
