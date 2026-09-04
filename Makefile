# One target, one binary. The tokeniser and the BWT codec are not optional:
# karyozip.c calls paydecomp unconditionally, so a build without -DKARYOZIP_TOK
# does not link. Everything vendored lives beside this file.
CC      ?= gcc
CFLAGS  ?= -O3 -std=gnu11 -Wall -Wextra -fopenmp
# VERSION is bzip3's, which libbz3.c needs; the KST format version is in
# the magic. They are not the same number and never were.
BZIP3_VERSION ?= 1.5.1
CFLAGS  += -DKARYOZIP_TOK -DVERSION='"$(BZIP3_VERSION)"' -I.
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

all: build/karyozip

build/karyozip: $(SRC)
	mkdir -p build
	$(CC) $(CFLAGS) $(ZSTD_CFLAGS) -o $@ $(SRC) $(LDFLAGS) $(ZSTD_LIBS) -lm -lpthread

test: build/karyozip
	KT=$(CURDIR)/build/karyozip tests/run.sh

clean:
	rm -rf build

.PHONY: all test clean
