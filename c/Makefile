# one binary, every codec, vendored beside
CC      ?= gcc
CFLAGS  ?= -O3 -std=gnu11 -Wall -Wextra -fopenmp
# VERSION is bzip3's, not the format's
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

# liblzma apart: its headers clash with htscodecs'
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

# libtest in all: a stale one misleads
all: build/karyozip build/libkaryozip.a build/libtest

# libzpaq: no OpenMP, JIT on x86-64 only
CXX     ?= g++
# or a libstdc++.a without g++
CXXLIB  ?= -lstdc++
ZPQ_DEF := -Dunix $(if $(filter x86_64,$(shell uname -m)),,-DNOJIT)
ZPQ_OBJ := build/zpaq/libzpaq.o build/zpaq/kzzpaq.o
build/zpaq/libzpaq.o: zpaq/libzpaq.cpp zpaq/libzpaq.h
	mkdir -p $(@D)
	$(CXX) -O3 -w $(ZPQ_DEF) -c $< -o $@
build/zpaq/kzzpaq.o: zpaq/kzzpaq.cpp zpaq/libzpaq.h
	mkdir -p $(@D)
	$(CXX) -O3 -Wall -Wextra -fopenmp $(ZPQ_DEF) -Izpaq -c $< -o $@

HDR := karyozip.h $(wildcard htscodecs/*.h bzip3/*.h libdeflate/*.h libdeflate/lib/*.h libdeflate/lib/*/*.h)
build/karyozip: $(SRC) $(HDR) $(XZ_OBJ) $(ZPQ_OBJ)
	mkdir -p build
	$(CC) $(CFLAGS) $(ZSTD_CFLAGS) -o $@ $(SRC) $(XZ_OBJ) $(ZPQ_OBJ) $(LDFLAGS) $(ZSTD_LIBS) -lm -lpthread -lz $(CXXLIB)

# the library: same sources, no main
OBJ := $(SRC:%.c=build/lib/%.o)
build/lib/%.o: %.c $(HDR)
	mkdir -p $(@D)
	$(CC) $(CFLAGS) $(ZSTD_CFLAGS) -Wno-unused-function -Wno-unused-variable -DKARYOZIP_LIB -c $< -o $@
build/libkaryozip.a: $(OBJ) $(XZ_OBJ) $(ZPQ_OBJ)
	rm -f $@ && $(AR) rcs $@ $(OBJ) $(XZ_OBJ) $(ZPQ_OBJ)

build/libtest: build/libkaryozip.a ../tests/libtest.c
	$(CC) $(CFLAGS) $(ZSTD_CFLAGS) -o $@ ../tests/libtest.c build/libkaryozip.a $(LDFLAGS) $(ZSTD_LIBS) -lm -lpthread -lz $(CXXLIB)

# the sanitizer build fuzz.sh runs
build/karyozip-san: $(SRC) $(HDR) $(XZ_SAN) $(ZPQ_OBJ)
	mkdir -p build
	$(CC) -O1 -g -fsanitize=address,undefined -std=gnu11 -fopenmp -DVERSION='"$(BZIP3_VERSION)"' -I. \
	  $(ZSTD_CFLAGS) -o $@ $(SRC) $(XZ_SAN) $(ZPQ_OBJ) $(LDFLAGS) $(ZSTD_LIBS) -lm -lpthread -lz $(CXXLIB)

test: build/karyozip build/libtest build/karyozip-san
	bash ../tests/run.sh build/karyozip build/libtest build/karyozip-san

clean:
	rm -rf build

.PHONY: all test clean
