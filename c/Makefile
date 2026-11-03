# one binary, every codec, vendored beside
CC      ?= gcc
# CFLAGS, CXXFLAGS, CPPFLAGS and LDFLAGS are the user's (a conda env sets all
# four): optimisation, tuning, include and library paths, on every rule. What
# the build needs is in the KZ_* below, always added after them.
CFLAGS   ?= -O3
CXXFLAGS ?= -O3
WARN     ?= -Wall -Wextra
# VERSION is bzip3's, not the format's
BZIP3_VERSION ?= 1.5.4
KZ_STD    := -std=gnu11
KZ_CFLAGS := $(KZ_STD) -fopenmp -DVERSION='"$(BZIP3_VERSION)"' -I.
SAN       := -O1 -g -fsanitize=address,undefined -U_FORTIFY_SOURCE
KZ_LIBS    = $(ZSTD_LIBS) -lm -lpthread -lz $(CXXLIB)

# libtest in all: a stale one misleads
all: build/karyozip build/libkaryozip.a build/libtest

# zstd: pkg-config when it knows libzstd, else plain -lzstd; ZSTD_CFLAGS and
# ZSTD_LIBS on the command line override both
PKG_CONFIG ?= pkg-config
ifeq ($(shell command -v $(PKG_CONFIG) 2>/dev/null),)
  $(warning $(PKG_CONFIG) not found: zstd from the compiler's default paths, set ZSTD_CFLAGS and ZSTD_LIBS to choose one)
  ZSTD_CFLAGS ?=
  ZSTD_LIBS   ?= -lzstd
else
  ifeq ($(origin ZSTD_CFLAGS),undefined)
    ZSTD_CFLAGS := $(shell $(PKG_CONFIG) --cflags libzstd 2>/dev/null)
  endif
  ifeq ($(origin ZSTD_LIBS),undefined)
    ZSTD_LIBS := $(shell $(PKG_CONFIG) --libs libzstd 2>/dev/null || echo -lzstd)
  endif
endif
ifneq ($(CONDA_PREFIX),)
  KZ_CFLAGS  += -I$(CONDA_PREFIX)/include
  KZ_LDFLAGS += -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib
endif
# zstd does not promise the same bytes across releases, and the Rust pins
# zstd-sys 2.1.2, which bundles zstd 1.5.7: check the zstd.h the C builds
# with. ZSTD_WANT= skips the check.
ZSTD_WANT ?= 1.5.7
zstd-check:
ifneq ($(ZSTD_WANT),)
	@v=$$(printf '#include <zstd.h>\nZSTD_VERSION_STRING\n' | \
	  $(CC) $(CPPFLAGS) $(CFLAGS) $(KZ_CFLAGS) $(ZSTD_CFLAGS) -E -P -x c - 2>/dev/null | tail -n 1); \
	if [ "$$v" != '"$(ZSTD_WANT)"' ]; then \
	  echo "karyozip: zstd.h is $${v:-not found}, not \"$(ZSTD_WANT)\": the C would not write the bytes the Rust" >&2; \
	  echo "  writes with its zstd $(ZSTD_WANT); point ZSTD_CFLAGS and ZSTD_LIBS at zstd $(ZSTD_WANT)" >&2; \
	  exit 1; fi
endif

SRC := karyozip.c \
       htscodecs/tokenise_name3.c htscodecs/rANS_static4x16pr.c \
       htscodecs/rANS_static32x16pr.c htscodecs/arith_dynamic.c \
       htscodecs/rle.c htscodecs/pack.c htscodecs/utils.c \
       bzip3/libbz3.c \
       libdeflate/lib/deflate_decompress.c libdeflate/lib/crc32.c libdeflate/lib/utils.c \
       libdeflate/lib/x86/cpu_features.c libdeflate/lib/arm/cpu_features.c

# the --help text, the Rust's too, as one C string literal without the last newline
build/usage.inc: usage.txt
	mkdir -p $(@D)
	sed -e 's/\\/\\\\/g' -e 's/"/\\"/g' -e 's/^/"/' -e 's/$$/\\n"/' -e '$$s/\\n"$$/"/' $< > $@

# liblzma apart: its headers clash with htscodecs'. The vendored files are
# listed in xz/sources.txt, which rust/build.rs reads too; kzlzma.c is ours,
# so it keeps its warnings
XZ_SRC := $(shell sed 's/\#.*//' xz/sources.txt)
XZ_FLAGS := $(KZ_STD) -DHAVE_CONFIG_H -Ixz -Ixz/src/common \
       $(addprefix -Ixz/src/liblzma/, api common check lz lzma rangecoder delta simple)
XZ_HDR := $(shell find xz -name '*.h')
XZ_OBJ := $(XZ_SRC:%.c=build/xz/%.o) build/xz/xz/kzlzma.o
build/xz/%.o: %.c $(XZ_HDR)
	mkdir -p $(@D)
	$(CC) $(CPPFLAGS) $(CFLAGS) -w $(XZ_FLAGS) -c $< -o $@
build/xz/xz/kzlzma.o: xz/kzlzma.c $(XZ_HDR)
	mkdir -p $(@D)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARN) $(XZ_FLAGS) -c $< -o $@
XZ_SAN := $(XZ_SRC:%.c=build/xz-san/%.o) build/xz-san/xz/kzlzma.o
build/xz-san/%.o: %.c $(XZ_HDR)
	mkdir -p $(@D)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(SAN) -w $(XZ_FLAGS) -c $< -o $@
build/xz-san/xz/kzlzma.o: xz/kzlzma.c $(XZ_HDR)
	mkdir -p $(@D)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(SAN) $(WARN) $(XZ_FLAGS) -c $< -o $@

# libzpaq: no OpenMP, JIT on x86-64 targets only: the compiler's target, not the host
CXX     ?= g++
# or a libstdc++.a without g++
CXXLIB  ?= -lstdc++
ARCH    ?= $(shell $(CC) -dumpmachine)
ZPQ_DEF := -Dunix $(if $(filter x86_64-%,$(ARCH)),,-DNOJIT)
ZPQ_OBJ := build/zpaq/libzpaq.o build/zpaq/kzzpaq.o
build/zpaq/libzpaq.o: zpaq/libzpaq.cpp zpaq/libzpaq.h
	mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -w $(ZPQ_DEF) -c $< -o $@
build/zpaq/kzzpaq.o: zpaq/kzzpaq.cpp zpaq/libzpaq.h
	mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(WARN) -fopenmp $(ZPQ_DEF) -Izpaq -c $< -o $@

HDR := karyozip.h xz/kzlzma.h $(wildcard htscodecs/*.h bzip3/*.h libdeflate/*.h libdeflate/lib/*.h libdeflate/lib/*/*.h)
build/karyozip: $(SRC) $(HDR) build/usage.inc $(XZ_OBJ) $(ZPQ_OBJ) | zstd-check
	mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARN) $(KZ_CFLAGS) $(ZSTD_CFLAGS) -o $@ $(SRC) $(XZ_OBJ) $(ZPQ_OBJ) $(LDFLAGS) $(KZ_LDFLAGS) $(KZ_LIBS)

# the library: same sources, no main
OBJ := $(SRC:%.c=build/lib/%.o)
build/lib/%.o: %.c $(HDR) build/usage.inc | zstd-check
	mkdir -p $(@D)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARN) $(KZ_CFLAGS) $(ZSTD_CFLAGS) -Wno-unused-function -Wno-unused-variable -DKARYOZIP_LIB -c $< -o $@
build/libkaryozip.a: $(OBJ) $(XZ_OBJ) $(ZPQ_OBJ)
	rm -f $@ && $(AR) rcs $@ $(OBJ) $(XZ_OBJ) $(ZPQ_OBJ)

build/libtest: build/libkaryozip.a ../tests/libtest.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARN) $(KZ_CFLAGS) $(ZSTD_CFLAGS) -o $@ ../tests/libtest.c build/libkaryozip.a $(LDFLAGS) $(KZ_LDFLAGS) $(KZ_LIBS)

# the sanitizer build fuzz.sh runs
build/karyozip-san: $(SRC) $(HDR) build/usage.inc $(XZ_SAN) $(ZPQ_OBJ) | zstd-check
	mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) $(SAN) $(KZ_CFLAGS) \
	  $(ZSTD_CFLAGS) -o $@ $(SRC) $(XZ_SAN) $(ZPQ_OBJ) $(LDFLAGS) $(KZ_LDFLAGS) $(KZ_LIBS)

test: build/karyozip build/libtest build/karyozip-san
	bash ../tests/run.sh build/karyozip build/libtest build/karyozip-san

clean:
	rm -rf build

.PHONY: all test clean zstd-check
