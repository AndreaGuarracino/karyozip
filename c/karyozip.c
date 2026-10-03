/* karyozip: the KaryoScope Track (KST) format */
#define _GNU_SOURCE
#include <stdio.h>
#include <sys/resource.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>
#include <zlib.h>
#include <time.h>
#include <zstd.h>
#include <lzma.h>
#include "karyozip.h"
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>
#include <fcntl.h>
#include <malloc.h>
#include <signal.h>
#include <unistd.h>
#include <pthread.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static uint64_t BLOCK = UINT64_MAX;      /* one block: the smallest file */
static int BLOCKSET = 0;                 /* --block was given */
/* 0 means use every core; --threads overrides */
static int NTHREADS = 0;
static int LEVEL = 17;
static int LEVELSET = 0;              /* kt_level was called: do not override */
static int LEVEL_L = -1, LEVEL_I = -1;   /* per stream, chosen by the search */

static const char MAGIC[8] = "KST00001";   /* KaryoScope Track, format version 1 */
#define IDXMAGIC 0x313058444954534bULL   /* "KSTIDX01", the index trailer */
#define KARYOZIP_VERSION "0.2.0"

static int PROF = 0;
static double PH[10];
static const char *PHN[10] = {"parse", "compress", "names", "index",
                             "smooth", "parse:lines", "parse:append", "write", "codec search", "layout search"};
static inline double nowsec(void) {
    struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}
/* the high-water mark is read from /proc at most every 10 ms: a span per
   sequence over 3,000,000 sequences made a 1.7 s smooth take 104 s */
#define PH_ADD(i, t0) do { if (PROF) { double _n = nowsec(); PH[i] += _n - (t0); \
    if (_n - PHLAST > 0.01) { PHRSS[i] = rss_now(); PHLAST = _n; } } } while (0)
static size_t PHRSS[10]; static double PHLAST;
static double WALL0 = 0;
/* --threads, else OMP_NUM_THREADS, else one per core capped at 8. The cap is
   measured, not guessed: encode is about half serial, so past 8 threads the
   extra ones only spin at barriers. On a 40M-row track, 16 cores:
     -t 4   5.09 s wall  13.7 s cpu
     -t 8   4.41 s wall  19.5 s cpu   <- the knee
     -t 16  4.69 s wall  35.3 s cpu
   An explicit --threads or OMP_NUM_THREADS overrides it. */
#define THRCAP 8
static int nthreads(void) {
    if (NTHREADS > 0) return NTHREADS;
#ifdef _OPENMP
    int m = omp_get_max_threads();
    if (m > THRCAP && !getenv("OMP_NUM_THREADS")) m = THRCAP;
    return m;
#else
    return 1;
#endif
}

static size_t rss_now(void) {          /* VmHWM: the true high-water mark */
    FILE *f = fopen("/proc/self/status", "r");
    if (!f) return 0;
    char l[256]; size_t kb = 0;
    while (fgets(l, sizeof l, f)) if (!strncmp(l, "VmHWM:", 6)) { kb = (size_t)strtoul(l + 6, NULL, 10); break; }
    fclose(f); return kb;
}

static void profdump(void) {
    /* WALL0 is still 0 when we die before the run starts, and the timings
       would then be the machine uptime. */
    if (!PROF || WALL0 == 0) return;
    double tot = 0;
    for (int i = 0; i < 10; i++) tot += PH[i];
    for (int i = 0; i < 10; i++)
        if (PH[i] > 0)
            fprintf(stderr, "  phase %-12s %8.3f s  %5.1f%%   hwm %5zu MB\n",
                    PHN[i], PH[i], 100 * PH[i] / (tot > 0 ? tot : 1), PHRSS[i] / 1024);

    struct rusage ru;
    if (!getrusage(RUSAGE_SELF, &ru)) {
        double cpu = ru.ru_utime.tv_sec + 1e-6 * ru.ru_utime.tv_usec
                   + ru.ru_stime.tv_sec + 1e-6 * ru.ru_stime.tv_usec;
        double wall = nowsec() - WALL0;
        fprintf(stderr, "  cpu %.3f s over %.3f s wall = %.2f cores busy"
                        " (%d allowed), peak rss %.2f GB\n",
                cpu, wall, wall > 0 ? cpu / wall : 0.0,
                nthreads(),
                ru.ru_maxrss / 1048576.0);
    }
}
static void die(const char *m) { fprintf(stderr, "karyozip: %s\n", m); exit(1); }
static void die2(const char *a, const char *b) { fprintf(stderr, "karyozip: %s%s\n", a, b); exit(1); }

static __thread size_t ALLOCMAX = 0;
static void ckalloc(size_t n) {
    if (ALLOCMAX && n > ALLOCMAX) die("corrupt file: an impossible allocation size");
}

/* the reporting form: a thread inside a parallel region must not exit */
static uint64_t bedcoord_e(const char *p, const char *end, const char *what,
                           const char **err) {
    if (p == end) { *err = "empty start or end column"; return 0; }
    uint64_t v = 0;
    for (; p < end; p++) {
        if (*p < '0' || *p > '9') { *err = what; return 0; }
        if (v > (UINT64_MAX - (uint64_t)(*p - '0')) / 10)
            { *err = "a coordinate does not fit in 64 bits"; return 0; }
        v = v * 10 + (uint64_t)(*p - '0');
    }
    return v;
}
static uint64_t bedcoord(const char *p, const char *end, const char *what) {
    const char *err = NULL;
    uint64_t v = bedcoord_e(p, end, what, &err);
    if (err) die(err);
    return v;
}

/* Every text input goes through zlib, which reads a plain file transparently
   and a .gz on the fly. fopencookie wraps it in a FILE*, so getline, fgets
   and fread below are unchanged. A path of NULL means stdin. */
static ssize_t zread_(void *c, char *b, size_t n) {
    int g = gzread((gzFile)c, b, (unsigned)n), e;
    gzerror((gzFile)c, &e);              /* a cut .gz hands over what it has with Z_BUF_ERROR: refuse it */
    return g < 0 || (e != Z_OK && e != Z_STREAM_END) ? -1 : g;
}
static int zclose_(void *c) { return gzclose((gzFile)c) == Z_OK ? 0 : -1; }
static FILE *zfopen(const char *path) {
    gzFile g = path ? gzopen(path, "rb") : gzdopen(dup(0), "rb");
    if (!g) return NULL;
    cookie_io_functions_t io = { zread_, NULL, NULL, zclose_ };
    FILE *f = fopencookie(g, "r", io);
    if (!f) gzclose(g);
    return f;
}

/* a line read by getline, its \n and a \r before it gone */
static void chomp(char *l, ssize_t *n) {
    if (*n && l[*n - 1] == '\n') l[--*n] = 0;
    if (*n && l[*n - 1] == '\r') l[--*n] = 0;
}

static uint64_t fnvbuf(uint64_t h, const void *p, size_t n) {
    for (size_t i = 0; i < n; i++) { h ^= ((const uint8_t *)p)[i]; h *= 0x100000001b3ULL; }
    return h;
}
#define FNV0 0xcbf29ce484222325ULL
/* one row into the row hash: name hash, start, end, label hash, as words */
static inline uint64_t rowmix(uint64_t h, uint64_t nh, uint64_t st, uint64_t en, uint64_t lh) {
    h = (h ^ nh) * 0x100000001b3ULL; h = (h ^ st) * 0x100000001b3ULL;
    h = (h ^ en) * 0x100000001b3ULL; return (h ^ lh) * 0x100000001b3ULL;
}
static uint64_t fnvstr(const char *s) { return fnvbuf(FNV0, s, strlen(s)); }
static uint64_t ld64(const void *p) { uint64_t v; memcpy(&v, p, 8); return v; }

static long int_arg(const char *p, long lo, long hi, const char *what) {   /* digits only, in range */
    if (!p || !*p || strspn(p, "0123456789") != strlen(p) || strlen(p) > 18) die(what);
    long v = strtol(p, NULL, 10);
    if (v < lo || v > hi) die(what);
    return v;
}

static uint64_t mem_arg(const char *p) {   /* bytes, with a K, M, G or T suffix */
    if (!p || !*p) die("--sort=MEM must be positive");
    /* Digits only, like every other number on the command line. This used
       strtod, which quietly took " 1M", "+1M", "1e9" and "1.5M", though the
       README says a sign and an exponent are refused. */
    uint64_t val = 0; const char *q = p;
    for (; *q >= '0' && *q <= '9'; q++) {
        if (val > (UINT64_MAX - (uint64_t)(*q - '0')) / 10)
            die("--sort=MEM is larger than the address space");
        val = val * 10 + (uint64_t)(*q - '0');
    }
    if (q == p || !val) die("--sort=MEM must be positive");
    uint64_t mul = 1;
    if (*q) {
        switch (*q) {
        case 'k': case 'K': mul = 1ULL << 10; break;
        case 'm': case 'M': mul = 1ULL << 20; break;
        case 'g': case 'G': mul = 1ULL << 30; break;
        case 't': case 'T': mul = 1ULL << 40; break;
        default: die("--sort=MEM: use a plain number of bytes or a K, M, G or T suffix");
        }
        if (q[1]) die("--sort=MEM: nothing may follow the unit");
    }
    if (val > UINT64_MAX / mul) die("--sort=MEM is larger than the address space");
    uint64_t m = val * mul;
    if (m < (1u << 20)) die("--sort=MEM must be at least 1M");
    return m;
}
static uint64_t num_arg(const char *p, uint64_t max, const char *what) {
    uint64_t v = bedcoord(p, p + strlen(p), what);
    if (v > max) die(what);
    return v;
}

static void *xmalloc(size_t n) { ckalloc(n); void *p = malloc(n); if (!p) die("oom"); return p; }

static size_t mul(size_t n, size_t sz) { if (sz && n > SIZE_MAX / sz) die("allocation size overflows"); return n * sz; }
static void *xcalloc(size_t n, size_t sz) {
    if (sz && n > SIZE_MAX / sz) die("allocation size overflows");
    ckalloc(n * sz); void *p = calloc(n, sz); if (!p) die("oom"); return p; }
static void *xrealloc(void *p, size_t n) { ckalloc(n); void *q = realloc(p, n); if (!q) die("oom"); return q; }
static char *xstrdup(const char *s) { char *p = strdup(s); if (!p) die("oom"); return p; }

/* ---------------------------------------------------------------- varint */
static inline void vput(uint8_t **p, uint64_t v) {
    while (v >= 0x80) { *(*p)++ = (uint8_t)(v | 0x80); v >>= 7; }
    *(*p)++ = (uint8_t)v;
}

#define VGUARD 16

static __thread const uint8_t *PAYEND;
static inline uint64_t vget(const uint8_t **p) {
    const uint8_t *q = *p;      /* VGUARD bytes spare: bound once, not per byte */
    uint64_t v = 0; int s = 0; uint8_t b;
    do { b = *q++;
         if (s < 64) v |= (uint64_t)(b & 0x7f) << s;
         s += 7;
         if (s > 70) die("corrupt stream: a varint does not terminate");
    } while (b & 0x80);
    if (PAYEND && q > PAYEND) die("corrupt block: a run reads past the payload");
    *p = q;
    return v;
}

static uint64_t vgetb(const uint8_t **p, const uint8_t *end) {
    uint64_t v = 0; int s = 0;
    for (;;) {
        if (*p >= end || s > 63) die("corrupt varint in metadata");
        uint8_t b = *(*p)++;
        v |= (uint64_t)(b & 0x7f) << s;
        if (!(b & 0x80)) return v;
        s += 7;
    }
}

/* ------------------------------------------------------------ growable buf */
typedef struct { uint8_t *d; size_t n, cap; } buf;
static void bput(buf *b, const void *s, size_t n) {
    if (!n) return;                 /* memcpy(dst, NULL, 0) is undefined */
    if (b->n + n > b->cap) {
        while (b->n + n > b->cap) b->cap = b->cap ? b->cap * 2 : (n > 4096 ? n : 4096);
        b->d = xrealloc(b->d, b->cap);
    }
    memcpy(b->d + b->n, s, n);
    b->n += n;
}
static void bvar(buf *b, uint64_t v) {
    uint8_t t[10], *p = t; vput(&p, v); bput(b, t, (size_t)(p - t));
}

/* --------------------------------------------------------------- encoding */
typedef struct { uint32_t blk, off, nruns; } rec;   /* 12 B, was 32 */

static char **labtab = NULL;
static int nlabtab = 0, labcap = 0;
static int32_t *labhash; static size_t LABHASH;   /* open addressed, grows at half full */
static int labhash_init = 0;

static size_t *labnlen;
static uint64_t *labh;             /* fnv of each label, for the row hash */
static uint64_t ROWHASH;           /* the rows as words, what binds a track used as --ref */
static int lastid = -1;

/* A track has a handful of distinct labels, 31 on a real chromosome track,
   but it cycles through them: measured on 4M real rows, the label equals the
   previous one 0.1% of the time, so the old single-slot cache never hit and
   every row paid a full FNV over the label plus a probe. This is a small
   direct-mapped cache on a four-byte key, checked with one memcmp. It changes
   no id: the table below still assigns them in first-seen order. */
#define LCACHE 512
static int lcache[LCACHE];
static int lcache_init = 0;
static int intern(const char *s, size_t n) {
    if (lastid >= 0 && labnlen[lastid] == n && !memcmp(labtab[lastid], s, n))
        return lastid;
    if (!lcache_init) { for (int i = 0; i < LCACHE; i++) lcache[i] = -1; lcache_init = 1; }
    uint32_t ck = (uint32_t)n * 2654435761u;
    if (n) ck ^= (uint32_t)(uint8_t)s[0] << 7 ^ (uint32_t)(uint8_t)s[n - 1] << 15
                ^ (uint32_t)(uint8_t)s[n >> 1] << 23;   /* the library allows an empty label */
    ck = (ck ^ (ck >> 13)) & (LCACHE - 1);
    int cid = lcache[ck];
    if (cid >= 0 && labnlen[cid] == n && !memcmp(labtab[cid], s, n)) return lastid = cid;
    if (!labhash_init) { LABHASH = 1 << 16; labhash = xmalloc(LABHASH * 4); memset(labhash, 0xff, LABHASH * 4); labhash_init = 1; }
    if ((size_t)nlabtab * 2 >= LABHASH) {           /* rehash at half full */
        LABHASH *= 4; free(labhash); labhash = xmalloc(LABHASH * 4); memset(labhash, 0xff, LABHASH * 4);
        for (int i = 0; i < nlabtab; i++) {
            uint64_t h = fnvbuf(FNV0, labtab[i], labnlen[i]) & (LABHASH - 1);
            while (labhash[h] >= 0) h = (h + 1) & (LABHASH - 1);
            labhash[h] = i;
        }
    }
    uint64_t k = fnvbuf(FNV0, s, n) & (LABHASH - 1);
    while (labhash[k] >= 0) {
        int id = labhash[k];
        if (labnlen[id] == n && !memcmp(labtab[id], s, n)) { lcache[ck] = id; return lastid = id; }
        k = (k + 1) & (LABHASH - 1);
    }
    if (nlabtab == labcap) {
        labcap = labcap ? labcap * 2 : 256;
        labtab = xrealloc(labtab, labcap * sizeof(char *));
        labnlen = xrealloc(labnlen, labcap * sizeof(size_t));
        labh = xrealloc(labh, labcap * sizeof(uint64_t));
    }
    labtab[nlabtab] = strndup(s, n);
    labnlen[nlabtab] = n;
    labh[nlabtab] = fnvbuf(FNV0, s, n);
    labhash[k] = nlabtab;
    lcache[ck] = nlabtab;
    return lastid = nlabtab++;
}

/* ZSTD_compress builds a level-17 match table per call, about 46 MB. Do NOT
   hold one context across calls: measured, that costs 42 MB more, because a
   freed table can be reused by the tokeniser and a live one cannot. */
static size_t zcl(buf *in, uint8_t **out, int level) {
    size_t bound = ZSTD_compressBound(in->n);
    *out = xmalloc(bound);
    size_t got = ZSTD_compress(*out, bound, in->d, in->n, level);
    if (ZSTD_isError(got)) die(ZSTD_getErrorName(got));
    return got;
}
static size_t zc(buf *in, uint8_t **out) { return zcl(in, out, LEVEL); }

/* one zstd frame, NUL-terminated; the refusals name what it held */
static char *unz(const void *src, size_t n, size_t *len, const char *what) {
    unsigned long long rb = ZSTD_getFrameContentSize(src, n);
    if (rb == ZSTD_CONTENTSIZE_UNKNOWN || rb == ZSTD_CONTENTSIZE_ERROR) die2(what, " has no usable frame size");
    char *o = xmalloc((size_t)rb + 1);
    if (ZSTD_isError(ZSTD_decompress(o, rb, src, n))) die2(what, " decompress failed");
    o[rb] = 0; *len = (size_t)rb; return o;
}

#define METALEVEL (LEVEL < 19 ? LEVEL : 19)

#define ICHUNK 65536

#define TOKFLAG (1ULL << 63)

#define CODEC_ARITH 0x10000
#define CODEC_BZIP3 0x20000           /* 0x20000|MB: bzip3, block size in MB */
#define CODEC_LZMA  0x40000           /* LZMA2, the settings lzmafilt pins */
#define DROPPED (1ULL << 60)          /* names deliberately not stored */
#define ALNFLAG (1ULL << 59)          /* rows deltaed against a reference */
#define DBFLAG  (1ULL << 58)          /* labels are featureIDs, not labels */
#define CARTFLAG (1ULL << 52)         /* some featureIDs are 2^26 + a vocabulary rank */
#define PERMFLAG (1ULL << 54)         /* names in input order + a permutation */
#define PERMTOK  (1ULL << 53)         /* that blob went through the tokeniser */
#define SIDEFLAG (1ULL << 57)
#define LENFLAG  (1ULL << 56)         /* lengths the reference implies are 0 */
#define TAILFLAG (1ULL << 55)
/* Every read's run count heads its own lengths, so the index blob holds
   nothing. One frame models the count beside the lengths it counts, which
   is smaller than a second frame of counts on its own. */
#define RUNFLAG  (1ULL << 61)
/* One frame, the label of a run beside its length: a repeated read pattern
   is one string the codec can match, not two halves in two frames. */
#define ILVFLAG  (1ULL << 51)
/* Lengths grouped by label id within a read. The labels stay in their own
   order and are the permutation, so nothing else is stored. */
#define GRPFLAG  (1ULL << 49)
static int RUNHEAD = 0;               /* the writer is using that layout */
static int ILV = 0;                   /* and whether the two frames are one */
static int GRP = 0;                   /* and whether the lengths are grouped */

static int NAMEMODE = 0;              /* 0 store, 1 drop */
static const char *NAMEFILE = NULL;
static const char *NAMESIN = NULL;
/* One flag, two kinds of file: the blob "--names drop:F" wrote, or the
 * sequences themselves. The track says which, so the file is never sniffed
 * (a pipe cannot be reopened). Resolved in ropen, see NAMESARG. */
static const char *NAMESARG = NULL;
/* the default searches hard; -e low trades size for encode time */
static int SLOW = 0;                  /* -e high: every codec over the whole file */
#ifndef DEFSAMPLE
#define DEFSAMPLE 16                  /* MB the default ranks every codec on */
#endif
static int FASTMODE = 0;              /* -e low: the sampled search */
static int SORT = 0;                  /* --sort: order records by label tuple */
static uint64_t SORTBP;               /* end of the run --sort last took */
static uint64_t SORTMEM = 0;
static const char *DBFEATS = NULL;    /* --db features.tsv */
static const char *HIERTAB = NULL;    /* --hier hierarchy.tsv, with --db */
static int DBIDS = 0;
static const char *ALNPAF = NULL;     /* --aln alignment.paf */
static uint64_t ALNHIT = 0;           /* rows the reference already answered */
static size_t RAWHELD = 0;            /* raw block payload held in memory at the peak */
static int RAW = 0;
static const char *OUTDIR = NULL;   /* --out-dir DIR, for decode --db */
static const char *SETNAME = NULL;  /* --set NAME: one set of a --db track */
static int INDEXING = 0;            /* the index command: names only, no companions */
static int ORDNAMES = 0;              /* --names ordinal: accept dropped names */
static const char *ALNREF = NULL;     /* --ref reference track */
static const char *ALNNAMES = NULL;   /* --names FILE: ordinal N is line N+1 */
static char **ORD; static size_t nord;

static uint64_t BYT_LENS = 0, BYT_LABS = 0;
static int IDCODEC = -1, LENCODEC = -1;   /* the search sets them */
#define INTFLAG (1ULL << 62)
#define FLDFLAG (1ULL << 50)   /* names as field deltas */

static size_t intnames(buf *names, uint64_t n, uint8_t **out) {
    const char *p = (const char *)names->d, *end = p + names->n;
    buf v = {0};
    int64_t prev = 0;
    for (uint64_t i = 0; i < n; i++) {
        if (p == end || *p == '\n' || (*p == '0' && p + 1 < end && p[1] != '\n'))
            { free(v.d); return 0; }             /* empty, or a leading zero */
        int64_t x = 0;
        while (p < end && *p != '\n') {
            if (*p < '0' || *p > '9' || x > (int64_t)1 << 50) { free(v.d); return 0; }
            x = x * 10 + (*p++ - '0');
        }
        if (p == end || x > (int64_t)1 << 50) { free(v.d); return 0; }   /* what the reader refuses */
        p++;
        int64_t d = x - prev; prev = x;

        bvar(&v, ((uint64_t)d << 1) ^ (uint64_t)(d >> 63));
    }
    if (p != end) { free(v.d); return 0; }
    size_t z = zc(&v, out); free(v.d);
    return z;
}

static char *intnamesd(const uint8_t *b, size_t n, uint64_t cnt, size_t *outn) {
    size_t rb; uint8_t *v = (uint8_t *)unz(b, n, &rb, "int names");
    char *o = xmalloc(mul(cnt, 21) + 1);

    const uint8_t *q = v, *qe = v + rb; size_t w = 0; uint64_t cur = 0;
    for (uint64_t i = 0; i < cnt; i++) {
        uint64_t u = vgetb(&q, qe);
        cur += (uint64_t)((int64_t)(u >> 1) ^ -(int64_t)(u & 1));
        if (cur > (1ull << 50)) die("corrupt names: an integer name past 2^50");
        w += (size_t)sprintf(o + w, "%llu\n", (unsigned long long)cur);
    }
    if (q != qe) die("corrupt names: the integer stream has more than the count");
    free(v); o[w] = 0; *outn = w; return o;
}

/* One list. The default tries all of it and spends whatever that costs;
   -e low is the subset marked 1, the codecs whose encode and decode stay
   reasonable, so it can never reach a candidate the default did not. */
typedef struct { int codec; int fast; } kcand;
static const kcand KCAND[] = {
    { 0x100, 1 }, { 0x80, 0 }, { 0x40, 0 }, { 0xc0, 0 },
    { 0x8, 0 }, { 0x88, 0 }, { 0x48, 0 }, { 0xc8, 0 },
    { 0x4, 0 }, { 0x84, 0 }, { 0x44, 0 }, { 0xc4, 0 },
    { 0xc, 0 }, { 0x8c, 0 }, { 0x4c, 0 }, { 0xcc, 0 },
    { 0x1, 1 }, { 0x81, 0 }, { 0x41, 1 }, { 0xc1, 0 },
    { 0x9, 0 }, { 0x89, 0 }, { 0x49, 0 }, { 0xc9, 0 },
    { 0x5, 0 }, { 0x85, 0 }, { 0x45, 0 }, { 0xc5, 0 },
    { 0xd, 0 }, { 0x8d, 0 }, { 0x4d, 0 }, { 0xcd, 0 },
    /* STRIPE carries its stripe count in bits 8 and up */
    { 0x208, 0 }, { 0x209, 1 }, { 0x408, 0 }, { 0x409, 0 },
    { CODEC_ARITH, 1 }, { CODEC_ARITH | 0x80, 0 },
    { CODEC_ARITH | 0x40, 0 }, { CODEC_ARITH | 0xc0, 0 },
    { CODEC_ARITH | 1, 1 }, { CODEC_ARITH | 0x81, 0 },
    { CODEC_ARITH | 0x41, 0 }, { CODEC_ARITH | 0xc1, 0 },
    { CODEC_BZIP3 | 8, 0 }, { CODEC_BZIP3 | 16, 1 }, { CODEC_BZIP3 | 32, 0 },
    { CODEC_BZIP3 | 64, 1 }, { CODEC_BZIP3 | 128, 0 }, { CODEC_BZIP3 | 256, 0 },
    { CODEC_LZMA, 0 },
};
#define KCAND_N (sizeof KCAND / sizeof *KCAND)

/* the mode's candidates, into `out`, and how many */
static size_t candlist(int *out) {
    size_t n = 0;
    for (size_t i = 0; i < KCAND_N; i++)
        if (!FASTMODE || KCAND[i].fast) out[n++] = KCAND[i].codec;
    return n;
}

static size_t paycomp(const uint8_t *src, size_t n, int codec, uint8_t **out);
static size_t paydecomp(const uint8_t *src, size_t n, int codec, uint8_t *dst, size_t raw);

/* Read names split on ':' into fields that are either the same in every name
   or a plain integer. The integers go down as zigzag deltas, interleaved so a
   read's fields sit together, and bzip3 takes the lot: on 3444234 Illumina
   names that is 2301419 B against the tokeniser's 2818102. */
static char FLDSEP = ':';   /* ':' or '-', whichever splits every name alike */
#define FLDMAX 16
static size_t fldnames1(buf *names, uint64_t n, uint8_t **out);

/* ':' splits an Illumina name, '-' a UUID: try both, keep the smaller */
static size_t fldnames(buf *names, uint64_t n, uint8_t **out) {
    uint8_t *a = NULL, *b = NULL;
    FLDSEP = ':'; size_t za = fldnames1(names, n, &a);
    FLDSEP = '-'; size_t zb = fldnames1(names, n, &b);
    if (zb && (!za || zb < za)) { free(a); *out = b; return zb; }
    free(b); *out = a; return za;
}

static size_t fldnames1(buf *names, uint64_t n, uint8_t **out) {
    if (!n) return 0;
    const char *p = (const char *)names->d, *end = p + names->n;
    const char *f0[FLDMAX]; size_t l0[FLDMAX]; int nf = 0, isnum[FLDMAX];
    {   const char *q = p, *nl = memchr(q, '\n', (size_t)(end - q));
        if (!nl) return 0;
        while (q <= nl) {
            const char *s = memchr(q, FLDSEP, (size_t)(nl - q));
            const char *e2 = s ? s : nl;
            if (nf == FLDMAX) return 0;
            f0[nf] = q; l0[nf] = (size_t)(e2 - q); nf++;
            if (!s) break;
            q = s + 1;
        }
    }
    int ishex[FLDMAX];
    for (int i = 0; i < nf; i++) { isnum[i] = 1; ishex[i] = !(l0[i] & 1) && l0[i] > 0 && l0[i] <= 64; }   /* the reader takes 64 */
    int *same = xmalloc((size_t)nf * sizeof *same);
    for (int i = 0; i < nf; i++) same[i] = 1;
    int64_t *prev = xcalloc((size_t)nf, sizeof *prev);
    buf v = {0};
    /* one pass: the fields of every name, checked and delta'd as we go */
    const char *q = p;
    for (uint64_t r = 0; r < n; r++) {
        const char *nl = q < end ? memchr(q, '\n', (size_t)(end - q)) : NULL;
        if (!nl) { free(v.d); free(same); free(prev); return 0; }
        const char *c = q; int k = 0;
        while (c <= nl && k < nf) {
            const char *s = memchr(c, FLDSEP, (size_t)(nl - c));
            const char *e2 = s ? s : nl;
            size_t len = (size_t)(e2 - c);
            if (same[k] && (len != l0[k] || memcmp(c, f0[k], len))) same[k] = 0;
            if (isnum[k]) {
                if (!len || len > 18 || (len > 1 && c[0] == '0')) isnum[k] = 0;
                else for (size_t j = 0; j < len; j++)
                         if (c[j] < '0' || c[j] > '9') { isnum[k] = 0; break; }
            }
            if (ishex[k]) {
                if (len != l0[k]) ishex[k] = 0;
                else for (size_t j = 0; j < len; j++)
                         if (!((c[j] >= '0' && c[j] <= '9') || (c[j] >= 'a' && c[j] <= 'f')))
                             { ishex[k] = 0; break; }
            }
            if (!same[k] && !isnum[k] && !ishex[k]) { free(v.d); free(same); free(prev); return 0; }
            c = e2 + 1; k++;
            if (!s) break;
        }
        if (k != nf || c != nl + 1) { free(v.d); free(same); free(prev); return 0; }
        q = nl + 1;
    }
    /* the template, then a read's varying fields together */
    bvar(&v, (uint64_t)nf);
    bvar(&v, (uint64_t)(unsigned char)FLDSEP);
    for (int i = 0; i < nf; i++) {
        int ty = same[i] ? 0 : isnum[i] ? 1 : 2;      /* constant, integer, hex */
        bvar(&v, (uint64_t)ty);
        if (ty == 0) { bvar(&v, l0[i]); bput(&v, f0[i], l0[i]); }
        else if (ty == 2) bvar(&v, l0[i]);
    }
    q = p;
    for (uint64_t r = 0; r < n; r++) {
        const char *nl = memchr(q, '\n', (size_t)(end - q));
        const char *c = q; int k = 0;
        while (k < nf) {
            const char *s = memchr(c, FLDSEP, (size_t)(nl - c));
            const char *e2 = s ? s : nl;
            if (!same[k] && isnum[k]) {
                int64_t x = 0;
                for (const char *j = c; j < e2; j++) x = x * 10 + (*j - '0');
                int64_t d = x - prev[k]; prev[k] = x;
                bvar(&v, ((uint64_t)d << 1) ^ (uint64_t)(d >> 63));
            } else if (!same[k]) {                     /* hex, two nibbles a byte */
                for (size_t j = 0; j + 1 < l0[k]; j += 2) {
                    int hi = c[j] <= '9' ? c[j] - '0' : c[j] - 'a' + 10;
                    int lo = c[j + 1] <= '9' ? c[j + 1] - '0' : c[j + 1] - 'a' + 10;
                    uint8_t b8 = (uint8_t)(hi << 4 | lo); bput(&v, &b8, 1);
                }
            }
            c = e2 + 1; k++;
        }
        q = nl + 1;
    }
    free(same); free(prev);
    uint8_t *z = NULL; size_t zn = paycomp(v.d, v.n, CODEC_BZIP3 | 64, &z);
    if (!zn) { free(v.d); free(z); return 0; }
    uint8_t hdr[10], *hw = hdr; vput(&hw, v.n);   /* bzip3 needs the size back */
    size_t hn = (size_t)(hw - hdr);
    *out = xmalloc(hn + zn);
    memcpy(*out, hdr, hn); memcpy(*out + hn, z, zn);
    free(z); free(v.d);
    return hn + zn;
}

static char *fldnamesd(const uint8_t *b, size_t n, uint64_t cnt, size_t *outn) {
    const uint8_t *hp = b, *he = b + n;
    size_t raw = (size_t)vgetb(&hp, he);
    if (!raw || raw > (1ull << 40)) die("corrupt names: a bad field blob size");
    uint8_t *v = xmalloc(raw + 16);
    size_t got = paydecomp(hp, (size_t)(he - hp), CODEC_BZIP3 | 64, v, raw);
    if (got != raw) die("field names decompress failed");
    const uint8_t *q = v, *qe = v + got;
    uint64_t nf = vgetb(&q, qe);
    if (!nf || nf > FLDMAX) die("corrupt names: a bad field count");
    char sep = (char)vgetb(&q, qe);
    const uint8_t *cf[FLDMAX]; size_t cl[FLDMAX]; int var[FLDMAX];
    size_t fixed = 0; int nvar = 0;
    for (uint64_t i = 0; i < nf; i++) {
        var[i] = (int)vgetb(&q, qe);
        if (var[i] > 2) die("corrupt names: a bad field type");
        if (var[i] == 1) { nvar++; cf[i] = NULL; cl[i] = 0; }
        else if (var[i] == 2) { cl[i] = (size_t)vgetb(&q, qe); cf[i] = NULL;
                                if (cl[i] & 1 || !cl[i] || cl[i] > 64)
                                    die("corrupt names: a bad hex field width");
                                fixed += cl[i]; }
        else { cl[i] = (size_t)vgetb(&q, qe);
               if ((size_t)(qe - q) < cl[i]) die("corrupt names: a field runs past the blob");
               cf[i] = q; q += cl[i]; fixed += cl[i]; }
    }
    size_t per = fixed + nf + (size_t)nvar * 20;
    char *o = xmalloc(mul(cnt, per) + 1); size_t w = 0;
    int64_t *cur = xcalloc(nf, sizeof *cur);
    for (uint64_t r = 0; r < cnt; r++) {
        for (uint64_t i = 0; i < nf; i++) {
            if (i) o[w++] = sep;
            if (var[i] == 1) {
                uint64_t u = vgetb(&q, qe);
                cur[i] += (int64_t)(u >> 1) ^ -(int64_t)(u & 1);
                if (cur[i] < 0) die("corrupt names: a negative field");
                w += (size_t)sprintf(o + w, "%lld", (long long)cur[i]);
            } else if (var[i] == 2) {
                static const char HEX[] = "0123456789abcdef";
                for (size_t j = 0; j < cl[i] / 2; j++) {
                    if (q >= qe) die("corrupt names: a hex field runs past the blob");
                    o[w++] = HEX[*q >> 4]; o[w++] = HEX[*q & 15]; q++;
                }
            } else { memcpy(o + w, cf[i], cl[i]); w += cl[i]; }
        }
        o[w++] = '\n';
    }
    if (q != qe) die("corrupt names: the field stream has more than the count");
    free(cur); free(v); o[w] = 0; *outn = w; return o;
}

#define TOKMAX 9000000u

static int V1 = 0;
static const char *TMPD = NULL;         /* --tmp-dir, else TMPDIR, else /tmp */

#include "htscodecs/tokenise_name3.h"
#include "htscodecs/version.h"
#include "htscodecs/rANS_static4x16.h"
#include "htscodecs/arith_dynamic.h"
#include "bzip3/libbz3.h"

static size_t tokc(buf *in, uint8_t **out, uint32_t chunk) {
    *out = NULL;
    if (!in->n) return 0;
    const char *base = (const char *)in->d;

    size_t cap = 64, nchunk = 0, *off = xmalloc((cap + 1) * sizeof *off);
    off[0] = 0;
    for (size_t p2 = 0; p2 < in->n; ) {
        size_t end = p2; uint32_t c = 0;
        while (end < in->n && c < chunk) { if (base[end] == '\n') c++; end++; }
        if (nchunk + 1 >= cap) { cap *= 2; off = xrealloc(off, (cap + 1) * sizeof *off); }
        off[++nchunk] = end; p2 = end;
    }

    uint8_t **blob = xcalloc(nchunk, sizeof *blob);
    uint32_t *clen = xmalloc(mul(nchunk, 4));
    _Atomic int failed = 0;
    int nthr = nthreads();
    if (nthr > (int)nchunk) nthr = nchunk ? (int)nchunk : 1;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthr)
#endif
    for (size_t i = 0; i < nchunk; i++) {
        if (failed) continue;
        int ol = 0;
        uint8_t *o = tok3_encode_names((char *)base + off[i],
                                       (int)(off[i + 1] - off[i]), 9, 0, &ol, NULL);
        if (!o) { failed = 1; continue; }
        blob[i] = o; clen[i] = (uint32_t)ol;
    }

    { uint8_t *d = in->d; for (size_t i = 0; i < in->n; i++) if (!d[i]) d[i] = 10; }
    if (failed) {
        for (size_t i = 0; i < nchunk; i++) free(blob[i]);
        free(blob); free(clen); free(off);
        fprintf(stderr, "name tokeniser refused a chunk, using zstd\n");
        return 0;
    }
    size_t hdr = 8 + 8 * nchunk, total = 0;
    for (size_t i = 0; i < nchunk; i++) total += clen[i];
    uint8_t *b = xmalloc(hdr + total);
    *out = b;
    uint32_t h0 = chunk, h1 = (uint32_t)nchunk;
    memcpy(b, &h0, 4); memcpy(b + 4, &h1, 4);
    uint32_t *rl = (uint32_t *)(b + 8);
    for (size_t i = 0; i < nchunk; i++) rl[i] = (uint32_t)(off[i + 1] - off[i]);
    memcpy(b + 8 + 4 * nchunk, clen, 4 * nchunk);
    uint8_t *q = b + hdr;
    for (size_t i = 0; i < nchunk; i++) { memcpy(q, blob[i], clen[i]); q += clen[i]; free(blob[i]); }
    free(blob); free(clen); free(off);
    return hdr + total;
}

static char *tokd(const uint8_t *b, size_t n, size_t *outn) {
    if (n < 8) die("truncated tokenised names");
    uint32_t nchunk; memcpy(&nchunk, b + 4, 4);

    size_t hdr = 8 + 8 * (size_t)nchunk, raw = 0;
    if (nchunk > (n - 8) / 8 || hdr > n) die("corrupt tokenised names header");
    size_t *woff = xmalloc(((size_t)nchunk + 1) * sizeof *woff);
    size_t *coff = xmalloc(((size_t)nchunk + 1) * sizeof *coff);
    woff[0] = coff[0] = 0;
    for (uint32_t i = 0; i < nchunk; i++) {   /* unaligned in a packed alignment */
        uint32_t rl, cl; memcpy(&rl, b + 8 + 4 * (size_t)i, 4); memcpy(&cl, b + 8 + 4 * ((size_t)nchunk + i), 4);
        raw += rl;
        woff[i + 1] = woff[i] + rl;
        coff[i + 1] = coff[i] + cl;
    }
    if (coff[nchunk] > n - hdr) die("corrupt tokenised names: chunks overrun the stream");
    char *o = xmalloc(raw + 1);
    const uint8_t *p = b + hdr;
    _Atomic int failed = 0;
    int nthr = nthreads();                 /* one chunk: no pool, a query pays for it */
    if (nthr > (int)nchunk) nthr = nchunk ? (int)nchunk : 1;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthr)
#endif
    for (uint32_t i = 0; i < nchunk; i++) {
        if (failed) continue;
        uint32_t dl = 0;
        uint8_t *d = tok3_decode_names((uint8_t *)p + coff[i], (uint32_t)(coff[i + 1] - coff[i]), &dl);
        if (!d || dl != woff[i + 1] - woff[i]) { failed = 1; free(d); continue; }

        for (uint32_t j = 0; j < dl; j++) if (!d[j]) d[j] = '\n';
        memcpy(o + woff[i], d, dl); free(d);
    }
    if (failed) die("name tokeniser decode failed");
    free(woff); free(coff);
    o[raw] = 0; *outn = raw; return o;
}

static const uint32_t TCAND[] = { 250000, 1000000, 3000000, TOKMAX };
static uint32_t TOKPICK = 250000;        /* the chunk the last tokbest chose */
static size_t tokbest(buf *in, uint8_t **out) {
    size_t n = sizeof TCAND / sizeof *TCAND, best = 0;
    for (size_t i = 0; i < n; i++) {
        if (FASTMODE && i) continue;    /* -e low: the smallest chunk only */
        uint8_t *o2 = NULL;
        size_t b = tokc(in, &o2, TCAND[i]);
        /* smallest wins; a tie keeps the smaller chunk */
        if (b && (!best || b < best)) { free(*out); *out = o2; best = b; TOKPICK = TCAND[i]; }
        else free(o2);
    }
    return best;
}

/* zstd and integer deltas, on their own thread beside the tokeniser */
typedef struct { buf *b; uint64_t n; uint8_t *a1, *a2; size_t s1, s2; } zjob;
static void *zjob_run(void *p) {
    zjob *j = p; j->s1 = zc(j->b, &j->a1); j->s2 = intnames(j->b, j->n, &j->a2); return NULL;
}

/* smallest of zstd, integer deltas, tokeniser, field deltas; losers freed */
static size_t bestnames(buf *b, uint64_t n, uint8_t **out, int *kind) {
    uint8_t *a3 = NULL; size_t s3 = 0;
    double _t0 = nowsec();
    /* A plain thread, not an OpenMP section: inside a section the
       tokeniser's own parallel loop would be nested and run on one thread.
       tok3 overwrites every separator with NUL and puts it back, so it
       cannot share the buffer with a concurrent reader: it gets a copy. */
    buf tb = { xmalloc(b->n ? b->n : 1), b->n, b->n };
    memcpy(tb.d, b->d, b->n);
    zjob zj = { b, n, NULL, NULL, 0, 0 };
    pthread_t th; int spawned = !pthread_create(&th, NULL, zjob_run, &zj);
    if (!spawned) zjob_run(&zj);
    s3 = tokbest(&tb, &a3);
    if (spawned) pthread_join(th, NULL);
    free(tb.d);
    uint8_t *a1 = zj.a1, *a2 = zj.a2; size_t s1 = zj.s1, s2 = zj.s2;
    uint8_t *a4 = NULL; size_t s4 = fldnames(b, n, &a4);
    if (PROF) fprintf(stderr, "  names %zu B: candidates %.3f s\n", b->n, nowsec() - _t0);
    size_t best = s1; *kind = 0; *out = a1;
    if (s2 && s2 < best) { best = s2; *kind = 1; *out = a2; }
    if (s3 && s3 < best) { best = s3; *kind = 2; *out = a3; }
    if (s4 && s4 < best) { best = s4; *kind = 3; *out = a4; }
    if (*out != a4) free(a4);
    if (*out != a1) free(a1);
    if (*out != a2) free(a2);
    if (*out != a3) free(a3);
    return best;
}


static int nosz(int c) { unsigned u = (unsigned)c; return u && !(u & CODEC_BZIP3) && (u & 0x10); }

/* every field pinned: a preset that moved between liblzma releases would
   change the bytes, and the C and the Rust must write the same ones */
static void lzmafilt(lzma_filter *f, lzma_options_lzma *o) {
    memset(o, 0, sizeof *o);
    o->dict_size = 64u << 20; o->lc = 3; o->lp = 0; o->pb = 0;
    o->mode = LZMA_MODE_NORMAL; o->nice_len = 273;
    o->mf = LZMA_MF_BT4; o->depth = 0;
    f[0].id = LZMA_FILTER_LZMA2; f[0].options = o;
    f[1].id = LZMA_VLI_UNKNOWN;  f[1].options = NULL;
}

static size_t paycomp(const uint8_t *src, size_t n, int codec, uint8_t **out) {
    unsigned int ol;
    if (codec & CODEC_LZMA) {
        lzma_filter f[2]; lzma_options_lzma o; lzmafilt(f, &o);
        size_t ob = n + n / 16 + 4096, op = 0;
        *out = malloc(ob);
        if (!*out || lzma_raw_buffer_encode(f, NULL, src, n, *out, &op, ob) != LZMA_OK) {
            free(*out); *out = NULL; return 0; }
        return op;
    }
    /* htscodecs is 32-bit throughout; bzip3 is not */
    if (!(codec & CODEC_BZIP3) && n >= 0x7fffffffu) { *out = NULL; return 0; }   /* rANS and arith decode below INT_MAX */
    if (codec & CODEC_BZIP3) {

        size_t bb = bz3_bound(n); *out = malloc(bb);
        size_t got = bb;
        if (!*out || bz3_compress((uint32_t)(codec & 0xffff) << 20,
                                  src, *out, n, &got)) { free(*out); *out = NULL; return 0; }
        return got;
    }
    if (codec & CODEC_ARITH) {

        uint64_t ob = (uint64_t)n * 2 + 1024;   /* 32-bit here wraps near 2 GB */
        if (ob > 0xffffffffu) { *out = NULL; return 0; }
        ol = (unsigned int)ob; *out = malloc(ol);
        if (!*out) return 0;
        if (!arith_compress_to((unsigned char *)src, (unsigned int)n,
                               *out, &ol, codec & 0xffff)) { if (!(codec & 0x04)) free(*out); *out = NULL; return 0; }
    } else {
        ol = rans_compress_bound_4x16((unsigned int)n, codec); *out = malloc(ol);
        if (!*out) return 0;

        if (!rans_compress_to_4x16((unsigned char *)src, (unsigned int)n,
                                   *out, &ol, codec)) { free(*out); *out = NULL; return 0; }
    }
    return ol;
}

static size_t paydecomp(const uint8_t *src, size_t n, int codec, uint8_t *dst, size_t raw) {
    if (codec & CODEC_LZMA) {
        lzma_filter f[2]; lzma_options_lzma o; lzmafilt(f, &o);
        size_t ip = 0, op = 0;
        if (lzma_raw_buffer_decode(f, NULL, src, &ip, n, dst, &op, raw) != LZMA_OK
            || op != raw) return 0;
        return op;
    }
    if (codec & CODEC_BZIP3) {
        size_t got = raw;
        if (bz3_decompress(src, dst, n, &got) || got != raw) return 0;
        return got;
    }
    if (raw > 0xffffffffu) return 0;
    unsigned int ul = (unsigned int)raw;
    unsigned char *q = (codec & CODEC_ARITH)
        ? arith_uncompress_to((unsigned char *)src, (unsigned int)n, dst, &ul)
        : rans_uncompress_to_4x16((unsigned char *)src, (unsigned int)n, dst, &ul);
    if (!q || ul != raw) return 0;
    return ul;
}

static const char *NAMES_OUT = NULL;

static uint64_t TAILPRED = 0;
static int TAILSEEN = 0;   /* how many totals sampled so far */

#define TAILSKIP 64
static int TAILON = -1;          /* -1 undecided, decided at sequence TAILSKIP */
static uint64_t TAILSAMP[TAILSKIP];
static __thread int TAILDEC = 0;          /* the file being read uses it */
static __thread uint64_t TAILPRED_R = 0;  /* its predictor, from the header */
static __thread uint64_t TAILSKIP_R = 0;  /* sequences before this store a plain length */

/* one pending block, while the raw payload waits in the staging file */
typedef struct {
    size_t off, len, lens, rows, bp;   /* into the staging file, then its two frames */
    size_t aoff[2], alen[2];           /* the same for the two --aln spellings */
} blkmeta;

typedef struct ktw {
    buf names, lens, ids, sizes;
    buf ibx; size_t nidx; int ordok;
    uint64_t nblk, pend, rows;
    FILE *tmp;
    char cur[512]; size_t curlen;
    uint64_t currun, curoff, curblk, curbp, blkbp0;
    FILE *out;                    /* where the finished track goes */
    int rawfd;                    /* raw block payloads, staged not held */
    blkmeta *bm; size_t nraw, rawcap;
    buf lalt[2], rawalt[2];       /* plain; plus-strand sentinels only */
    uint64_t plen, pval; int phas, prev;
    uint64_t curnh;               /* fnv of cur, for the row hash */
} ktw;

/* the same lengths with each read's run count in front of its own */
static size_t splice_heads(const uint8_t *lens, size_t nl, const uint8_t *ib,
                           size_t ibn, uint8_t **out) {
    *out = xmalloc(nl + ibn + 16);
    const uint8_t *p = lens, *e = lens + nl, *q = ib, *qe = ib + ibn;
    uint8_t *w = *out;
    while (q < qe) {
        uint64_t k = vgetb(&q, qe);
        vput(&w, k);
        for (uint64_t j = 0; j < k; j++) {
            const uint8_t *s0 = p;
            if (p >= e) { free(*out); *out = NULL; return 0; }
            vgetb(&p, e);
            size_t len = (size_t)(p - s0);
            memcpy(w, s0, len); w += len;
        }
    }
    if (p != e) { free(*out); *out = NULL; return 0; }
    return (size_t)(w - *out);
}

/* One stream a read: its run count, then a label beside its own length.
   A read that repeats another is then one string, not two halves in two
   frames, and the count heads it so the index holds nothing. */
/* Labels counted down from 127 so they do not share bytes with the small
   lengths beside them: the BWT cannot see which side of a pair a byte is on,
   and a byte that means one thing is worth more than one that means two.
   Its own inverse, and the varint keeps its width. */
static inline uint64_t ilv_flip(uint64_t v) { return v < 128 ? 127 - v : v; }

static size_t ilv_pack(const uint8_t *lens, size_t nl, const uint8_t *ids, size_t ni,
                       const uint8_t *ib, size_t ibn, uint8_t **out) {
    *out = xmalloc(nl + ni + ibn + 16);
    const uint8_t *p = lens, *e = lens + nl, *q = ids, *qe = ids + ni;
    const uint8_t *r = ib, *re = ib + ibn;
    uint8_t *w = *out;
    while (r < re) {
        uint64_t k = vgetb(&r, re);
        vput(&w, k);
        for (uint64_t j = 0; j < k; j++) {
            if (q >= qe || p >= e) { free(*out); *out = NULL; return 0; }
            vput(&w, ilv_flip(vgetb(&q, qe)));
            const uint8_t *s0 = p; vgetb(&p, e);
            memcpy(w, s0, (size_t)(p - s0)); w += p - s0;
        }
    }
    if (p != e || q != qe) { free(*out); *out = NULL; return 0; }
    return (size_t)(w - *out);
}

/* Each read's lengths in label order: all of label 0's, then label 1's. The
   label stream stays as it was, so it says where each one goes back. */
/* sort key: the label above, the first-seen order below */
static int gcmp(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return x < y ? -1 : x > y;
}

static size_t grp_pack(const uint8_t *lens, size_t nl, const uint8_t *ids, size_t ni,
                       const uint8_t *ib, size_t ibn, uint8_t **out) {
    *out = xmalloc(nl + ibn + 16);
    const uint8_t *p = lens, *e = lens + nl, *q = ids, *qe = ids + ni, *r = ib, *re = ib + ibn;
    uint8_t *w = *out;
    size_t cap = 64;
    uint64_t *lv = xmalloc(mul(cap, 8)), *kv = xmalloc(mul(cap, 8));
    int bad = 0;
    while (r < re && !bad) {
        uint64_t k = vgetb(&r, re);
        vput(&w, k);
        if (k >> 32) { bad = 1; break; }
        if (k > cap) { cap = (size_t)k; lv = xrealloc(lv, mul(cap, 8)); kv = xrealloc(kv, mul(cap, 8)); }
        for (uint64_t j = 0; j < k; j++) {
            if (p >= e || q >= qe) { bad = 1; break; }
            lv[j] = vgetb(&p, e);
            uint64_t id = vgetb(&q, qe);
            if (id >> 32) { bad = 1; break; }
            kv[j] = (id << 32) | j;
        }
        if (bad) break;
        qsort(kv, (size_t)k, 8, gcmp);
        for (uint64_t j = 0; j < k; j++) vput(&w, lv[kv[j] & 0xffffffffu]);
    }
    free(lv); free(kv);
    if (bad || p != e || q != qe) { free(*out); *out = NULL; return 0; }
    return (size_t)(w - *out);
}

static void rawput(int fd, off_t off, const void *p, size_t n) {
    while (n) {
        ssize_t k = pwrite(fd, p, n, off);
        if (k <= 0) die("cannot stage a block: is the temp filesystem full?");
        p = (const char *)p + k; off += k; n -= (size_t)k;
    }
}
static void rawget(int fd, off_t off, void *p, size_t n) {
    while (n) {
        ssize_t k = pread(fd, p, n, off);
        if (k <= 0) die("cannot read back a staged block");
        p = (char *)p + k; off += k; n -= (size_t)k;
    }
}

static void ktw_flush(ktw *w) {
    if (!w->pend) return;
    if (w->nraw == w->rawcap) {
        w->rawcap = w->rawcap ? w->rawcap * 2 : 1 << 16;
        w->bm = xrealloc(w->bm, mul(w->rawcap, sizeof *w->bm));
    }
    blkmeta *m = w->bm + w->nraw;
    for (int v = 0; v < 2; v++) {
        m->aoff[v] = w->rawalt[v].n; m->alen[v] = w->lalt[v].n;
        bput(&w->rawalt[v], w->lalt[v].d, w->lalt[v].n); w->lalt[v].n = 0;
    }
    m->off = w->nraw ? w->bm[w->nraw - 1].off + w->bm[w->nraw - 1].len : 0;
    m->len = w->lens.n + w->ids.n;
    m->lens = w->lens.n; m->rows = w->pend;
    m->bp = w->blkbp0; w->blkbp0 = 0; w->nraw++;
    rawput(w->rawfd, (off_t)m->off, w->lens.d, w->lens.n);
    rawput(w->rawfd, (off_t)(m->off + w->lens.n), w->ids.d, w->ids.n);
    w->nblk++;
    w->lens.n = 0; w->ids.n = 0; w->pend = 0;
}

static uint64_t zz(int64_t v) { return ((uint64_t)v << 1) ^ (uint64_t)(v >> 63); }
static int64_t unzz(uint64_t v) { return (int64_t)(v >> 1) ^ -(int64_t)(v & 1); }

static void lensput(ktw *w, uint64_t v) {
    bvar(&w->lens, v);
    if (ALNPAF) { bvar(&w->lalt[0], w->plen); bvar(&w->lalt[1], w->prev ? w->plen : v); }
}

static void ktw_tail(ktw *w) {
    if (!w->phas) return;
    if (TAILSEEN < TAILSKIP) {     /* still sampling: store the plain length */
        TAILSAMP[TAILSEEN++] = w->curbp;
        if (TAILSEEN == TAILSKIP) {

            uint64_t srt[TAILSKIP];
            memcpy(srt, TAILSAMP, sizeof srt);
            for (int i = 1; i < TAILSKIP; i++) {       /* 64 values, insertion */
                uint64_t v = srt[i]; int j = i - 1;
                while (j >= 0 && srt[j] > v) { srt[j + 1] = srt[j]; j--; }
                srt[j + 1] = v;
            }
            TAILPRED = srt[TAILSKIP / 2];
            uint64_t dev = 0;
            for (int i = 0; i < TAILSKIP; i++)
                dev += TAILSAMP[i] > TAILPRED ? TAILSAMP[i] - TAILPRED
                                              : TAILPRED - TAILSAMP[i];
            TAILON = dev * 12 <= TAILPRED * TAILSKIP;
        }
        lensput(w, w->pval); w->phas = 0; return;
    }
    if (!TAILON) { lensput(w, w->pval); w->phas = 0; return; }
    bvar(&w->lens, zz((int64_t)(w->curbp - TAILPRED)));
    if (ALNPAF) { bvar(&w->lalt[0], zz((int64_t)(w->curbp - TAILPRED))); bvar(&w->lalt[1], zz((int64_t)(w->curbp - TAILPRED))); }
    w->phas = 0;
}

static void ktw_pushread(ktw *w) {
    if (!w->cur[0]) return;
    if (w->curblk > 0xffffffffu || w->curoff > 0xffffffffu || w->currun > 0xffffffffu)
        die("sequence index does not fit 32 bits");
    bvar(&w->ibx, w->currun);
    if (w->ordok) {                     /* still 0,1,2,...: nothing to store */
        char want[24];
        int wl = snprintf(want, sizeof want, "%zu", w->nidx);
        if (w->curlen == (size_t)wl && !memcmp(w->cur, want, w->curlen)) {
            w->nidx++;
            return;
        }

        for (size_t j = 0; j < w->nidx; j++) {
            char t[24]; int tl = snprintf(t, sizeof t, "%zu", j);
            bput(&w->names, t, (size_t)tl); bput(&w->names, "\n", 1);
        }
        w->ordok = 0;
    }
    w->nidx++;
    bput(&w->names, w->cur, w->curlen); bput(&w->names, "\n", 1);
}

static void ktw_close(ktw *w);

static ktw *ktw_open(FILE *out) {
    ktw *w = xcalloc(1, sizeof(ktw));
    w->ordok = 1; ROWHASH = FNV0;
    const char *td = TMPD; if (!td || !*td) td = getenv("TMPDIR"); if (!td || !*td) td = "/tmp";
    char *tp; if (asprintf(&tp, "%s/karyozip.XXXXXX", td) < 0) die("oom");
    int tfd = mkstemp(tp);
    if (tfd < 0) { fprintf(stderr, "karyozip: cannot create a temp file under %s (set TMPDIR)\n", td); exit(1); }
    unlink(tp); free(tp);
    w->out = out;
    w->tmp = fdopen(tfd, "w+b");
    if (!w->tmp) die("tmpfile");
    if (asprintf(&tp, "%s/karyozip.raw.XXXXXX", td) < 0) die("oom");
    w->rawfd = mkstemp(tp);
    if (w->rawfd < 0) { fprintf(stderr, "karyozip: cannot create a temp file under %s (set TMPDIR)\n", td); exit(1); }
    unlink(tp); free(tp);
    return w;
}

static __thread int LENDEC = 0;      /* the file being read leaves lengths out */
static uint64_t refend_at(const char *nm, uint64_t pos);
static uint32_t refrow_at(const char *nm, uint64_t st, uint64_t len);
static __thread int LASTREV;   /* strand of the last refend answer */
static const char *ALNNAME = NULL;    /* --sort writes ranks; the PAF wants names */

/* same: the caller already knows this row's name equals the previous row's,
   which the parse threads work out in parallel. It saves one memcmp per row,
   and the name changes only once per sequence. */
static void ktw_row(ktw *w, const char *name, size_t nlen,
                    uint64_t run_len, const char *label, size_t llen, uint64_t start,
                    int same) {
    if (!same && (nlen != w->curlen || memcmp(name, w->cur, nlen) != 0)) {
        if (start != 0) die("a sequence does not start at 0");
        ktw_tail(w);            /* the held row was this sequence's last */
        ktw_pushread(w);
        if (nlen >= sizeof w->cur) die("name too long");
        memcpy(w->cur, name, nlen); w->cur[nlen] = 0; w->curlen = nlen;
        w->curnh = fnvbuf(FNV0, name, nlen);
        w->currun = 0; w->curbp = 0;
        if (w->pend && w->pend >= BLOCK) ktw_flush(w);
        w->curblk = w->nblk; w->curoff = w->pend;
    }

    if (start != w->curbp) die("gap or overlap: start is not the previous end");
    uint64_t id = (uint64_t)intern(label, llen);   /* once, not once per check */
    ROWHASH = rowmix(ROWHASH, w->curnh, start, start + run_len, labh[id]);
    /* A label the reference already gives is stored as "=". Done here, on
       the row as written, so a --sort track orders on the input's own labels
       and its order is the same with or without --aln. */
    if (ALNPAF && refrow_at(ALNNAME ? ALNNAME : w->cur, start, run_len) == (uint32_t)id) {
        id = (uint64_t)intern("=", 1); ALNHIT++;
    }

    if (w->phas) { lensput(w, w->pval); w->phas = 0; }
    if (w->pend >= BLOCK) ktw_flush(w);            /* mid-sequence is fine */
    if (!w->pend) w->blkbp0 = w->curbp;   /* where this block opens */
    w->plen = run_len; w->phas = 1;
    w->pval = ALNPAF
              && refend_at(ALNNAME ? ALNNAME : w->cur, start) == start + run_len
              ? 0 : run_len;
    w->prev = LASTREV;
    bvar(&w->ids, id);
    w->pend++; w->currun++; w->rows++; w->curbp += run_len;
}

typedef struct { uint64_t q; uint32_t t, n; uint64_t tp; int rev; } seg;
typedef struct { char *name; seg *s; size_t n, cap; } qseg;
typedef struct { char *name; uint64_t *st, end; uint32_t *lab, *st32; size_t n, cap; int64_t ri;
                 /* st32 holds the starts once a whole contig is in and fits 32 bits: 8 bytes a run, not 12 */
                 _Atomic int filled;              /* every run is in st and lab */
                 uint64_t bp0, row0, *bend;       /* .kst per block: bases and row before the contig, */
                 uint32_t firstblk, lastblk;      /* its block range, and the position after each */
                 uint8_t *blkin; uint32_t nblkin; } rcon;

static qseg *QS; static size_t nqs, cqs;
static rcon *RC; static size_t nrc, crc;

static uint64_t ALNHASH;
static int REFLAZY;                  /* readers: decode a reference contig when first asked */
static void loadref_kst(const char *path);
static void reffill(rcon *rc);
static void rcompact(rcon *rc);
static rcon *rcon_get(const char *nm, int make) {
    static size_t *idx; static size_t im;           /* open addressed by hash */
    if (!idx) { im = 1 << 16; idx = xcalloc(im, sizeof *idx); }
    size_t j = fnvstr(nm) & (im - 1);
    while (idx[j]) { if (!strcmp(RC[idx[j] - 1].name, nm)) return RC + idx[j] - 1;
                     j = (j + 1) & (im - 1); }
    if (!make) return NULL;
    if (nrc * 2 + 2 > im) {                          /* grow and rehash, then place */
        free(idx); im *= 4; idx = xcalloc(im, sizeof *idx);
        for (size_t k = 0; k < nrc; k++) {
            size_t h = fnvstr(RC[k].name) & (im - 1);
            while (idx[h]) h = (h + 1) & (im - 1);
            idx[h] = k + 1;
        }
        j = fnvstr(nm) & (im - 1); while (idx[j]) j = (j + 1) & (im - 1);
    }
    idx[j] = nrc + 1;
    if (nrc == crc) { crc = crc ? crc * 2 : 64; RC = xrealloc(RC, crc * sizeof *RC); }
    RC[nrc] = (rcon){ .name = xstrdup(nm), .ri = -1, .filled = 1 };
    return RC + nrc++;
}
static qseg *qseg_get(const char *nm, int make) {
    static size_t *idx; static size_t im;           /* open addressed by hash */
    if (!idx) { im = 1 << 20; idx = xcalloc(im, sizeof *idx); }
    size_t j = fnvstr(nm) & (im - 1);
    while (idx[j]) { if (!strcmp(QS[idx[j] - 1].name, nm)) return QS + idx[j] - 1;
                     j = (j + 1) & (im - 1); }
    if (!make) return NULL;
    if (nqs == cqs) { cqs = cqs ? cqs * 2 : 1 << 16; QS = xrealloc(QS, cqs * sizeof *QS); }
    QS[nqs] = (qseg){ xstrdup(nm), NULL, 0, 0 };
    idx[j] = ++nqs;
    if (nqs * 2 > im) {                              /* grow and rehash */
        free(idx); im *= 4; idx = xcalloc(im, sizeof *idx);
        for (size_t k = 0; k < nqs; k++) {
            size_t h = fnvstr(QS[k].name) & (im - 1);
            while (idx[h]) h = (h + 1) & (im - 1);
            idx[h] = k + 1;
        }
    }
    return QS + nqs - 1;
}

static void loadref(const char *path) {
    {   FILE *k = fopen(path, "rb"); char m[8];
        int kst = k && fread(m, 1, 8, k) == 8 && !memcmp(m, MAGIC, 8);
        if (k) fclose(k);
        if (kst) { loadref_kst(path); return; }
    }
    FILE *f = zfopen(path);
    if (!f) die("cannot open the reference track");
    char *ln = NULL; size_t cap = 0; ssize_t n;
    rcon *cur = NULL; uint64_t badln = 0;
    uint64_t rh = FNV0, nh = 0;
    while ((n = getline(&ln, &cap, f)) > 0) {
        chomp(ln, &n);
        char *t1 = memchr(ln, '\t', (size_t)n); if (!t1) { badln++; continue; }
        char *t2 = memchr(t1 + 1, '\t', (size_t)(n - (t1 + 1 - ln))); if (!t2) { badln++; continue; }
        char *t3 = memchr(t2 + 1, '\t', (size_t)(n - (t2 + 1 - ln))); if (!t3) { badln++; continue; }
        *t1 = 0;
        if (t1 == ln) die("reference track: empty contig name");
        if (!cur || strcmp(ln, cur->name)) {
            if (rcon_get(ln, 0)) die("reference track: a contig's rows are not together");
            cur = rcon_get(ln, 1); nh = fnvbuf(FNV0, ln, (size_t)(t1 - ln));
        }
        {   uint64_t st0 = bedcoord(t1 + 1, t2, "reference start is not a number");
            if (st0 != (cur->n ? cur->end : 0)) die("reference track: rows are not sorted and gapless from 0");
        }
        if (cur->n == cur->cap) {
            cur->cap = cur->cap ? cur->cap * 2 : 1024;
            cur->st = xrealloc(cur->st, cur->cap * sizeof *cur->st);
            cur->lab = xrealloc(cur->lab, cur->cap * sizeof *cur->lab);
        }
        cur->st[cur->n] = bedcoord(t1 + 1, t2, "reference start is not a number");
        cur->lab[cur->n] = (uint32_t)intern(t3 + 1, strlen(t3 + 1));
        cur->end = bedcoord(t2 + 1, t3, "reference end is not a number");   /* the last run's end */
        if (cur->end <= cur->st[cur->n]) die("reference track: a run of no length");
        rh = rowmix(rh, nh, cur->st[cur->n], cur->end, labh[cur->lab[cur->n]]);
        cur->n++;
    }
    if (ferror(f)) die("read error on the reference track");
    free(ln); fclose(f); ALNHASH ^= rh;
    for (size_t i = 0; i < nrc; i++) rcompact(RC + i);
    if (badln) fprintf(stderr, "reference: %llu lines without four columns skipped\n", (unsigned long long)badln);
    size_t tot = 0; for (size_t i = 0; i < nrc; i++) tot += RC[i].n;
    fprintf(stderr, "reference: %zu contigs, %zu runs\n", nrc, tot);
}

static int segcmp(const void *a, const void *b) {
    const seg *x = a, *y = b;
    if (x->q != y->q) return x->q < y->q ? -1 : 1;
    if (x->t != y->t) return x->t < y->t ? -1 : 1;
    if (x->tp != y->tp) return x->tp < y->tp ? -1 : 1;
    if (x->rev != y->rev) return x->rev - y->rev;
    return x->n < y->n ? -1 : x->n > y->n;
}

static void seghash(void) {
    uint64_t sum = 0;
    for (size_t i = 0; i < nqs; i++) {
        uint64_t h = FNV0;
        h = fnvbuf(h, QS[i].name, strlen(QS[i].name) + 1);
        for (size_t k = 0; k < QS[i].n; k++) { const seg *g = QS[i].s + k;
            h = fnvbuf(h, &g->q, 8); h = fnvbuf(h, &g->t, 4); h = fnvbuf(h, &g->n, 4); h = fnvbuf(h, &g->tp, 8); h = fnvbuf(h, &g->rev, 4); }
        if (QS[i].n) sum += h;
    }
    ALNHASH ^= sum;
}

static void loadpaf(const char *path) {
    FILE *f = zfopen(path);
    if (!f) die("cannot open the PAF");
    char *ln = NULL; size_t cap = 0; ssize_t n;
    size_t nblk = 0; uint64_t badpaf = 0, nocg = 0, noref = 0;
    while ((n = getline(&ln, &cap, f)) > 0) {
        chomp(ln, &n);

        char *fl[12]; int nf = 0;
        char *p = ln;
        for (; nf < 12 && p; nf++) { fl[nf] = p; char *q = strchr(p, '\t');
                                     if (q) { *q = 0; p = q + 1; } else p = NULL; }
        if (nf < 12 || !p) { badpaf++; continue; }
        const char *cg = NULL;
        for (char *q = p; q; ) {                     /* cg:Z: is in the tags */
            if (!strncmp(q, "cg:Z:", 5)) { cg = q + 5; break; }
            char *t = strchr(q, '\t'); q = t ? t + 1 : NULL;
        }
        if (!cg) { nocg++; continue; }
        if ((fl[4][0] != '+' && fl[4][0] != '-') || fl[4][1]) die("PAF strand is not + or -");
        int rev = fl[4][0] == '-';

        uint64_t qs = bedcoord(fl[2], fl[2] + strlen(fl[2]), "PAF query start is not a number");
        uint64_t qe = bedcoord(fl[3], fl[3] + strlen(fl[3]), "PAF query end is not a number");
        uint64_t tp = bedcoord(fl[7], fl[7] + strlen(fl[7]), "PAF target start is not a number");
        if (qe < qs) die("PAF query end is before its start");
        uint64_t qp = rev ? qe : qs;
        rcon *rc = rcon_get(fl[5], 0);
        if (!rc) { noref++; continue; }               /* target absent from ref */
        uint32_t tix = (uint32_t)(rc - RC);
        qseg *Q = NULL;
        /* A CIGAR is checked, not trusted: */
        uint64_t qcons = 0, tcons = 0;
        for (const char *c = cg; *c && *c != '\t'; ) {
            if (*c < '0' || *c > '9') die("PAF CIGAR: an operation has no length");
            uint64_t v = 0;
            while (*c >= '0' && *c <= '9') {
                if (v > (UINT64_MAX - 9) / 10) die("PAF CIGAR: a length overflows");
                v = v * 10 + (uint64_t)(*c++ - '0'); }
            if (!v) die("PAF CIGAR: a zero-length operation");
            if (v > 0xffffffffu) die("PAF CIGAR: an operation longer than 4 GB");
            char op = *c++;
            if (op == 'M' || op == '=' || op == 'X') {
                if (rev && qp < v) die("PAF CIGAR runs before the query start");
                if (!Q) Q = qseg_get(fl[0], 1);
                if (Q->n == Q->cap) { Q->cap = Q->cap ? Q->cap * 2 : 8;
                                      Q->s = xrealloc(Q->s, Q->cap * sizeof *Q->s); }
                Q->s[Q->n++] = (seg){ rev ? qp - v : qp, tix, (uint32_t)v, tp, rev };
                qp = rev ? qp - v : qp + v; tp += v; qcons += v; tcons += v; nblk++;
            } else if (op == 'I' || op == 'S') {
                if (rev && qp < v) die("PAF CIGAR runs before the query start");
                qp = rev ? qp - v : qp + v; qcons += v;
            } else if (op == 'D' || op == 'N') { tp += v; tcons += v; }
            else if (op == 'H' || op == 'P') { }
            else die("PAF CIGAR: unknown operation");
        }
        if (qcons != qe - qs)
            die("PAF CIGAR consumes a different query length than the record says");
        {   uint64_t te = bedcoord(fl[8], fl[8] + strlen(fl[8]), "PAF target end is not a number");
            uint64_t ts = tp - tcons;
            if (te < ts || tcons != te - ts)
                die("PAF CIGAR consumes a different target length than the record says"); }
    }
    if (ferror(f)) die("read error on a companion file");
    free(ln); fclose(f);
    for (size_t i = 0; i < nqs; i++) qsort(QS[i].s, QS[i].n, sizeof *QS[i].s, segcmp);
    seghash();
    if (badpaf || nocg || noref)
        fprintf(stderr, "alignment: skipped %llu malformed PAF lines, %llu without a cg:Z: tag, "
                        "%llu whose target is not in the reference\n",
                (unsigned long long)badpaf, (unsigned long long)nocg, (unsigned long long)noref);
    if (!nblk)
        die(noref ? "alignment: no PAF record maps to a contig of the reference track"
                  : "alignment: no usable PAF record, so --aln would silently do "
                    "nothing. Check the PAF carries cg:Z: CIGAR tags");
    fprintf(stderr, "alignment: %zu query contigs, %zu collinear blocks\n", nqs, nblk);
}

/* --------------------------------------------------- packed alignment Only what --aln reads: */
#define ALNMAGIC "KSTAL3\n"   /* carries the reference contig-order hash */
#define ALNMAGLEN 6
#define ALNTOK 0x7fff0000   /* the name stream went through the tokeniser */

static uint64_t refhash(void) {   /* the reference contig names, in order */
    uint64_t h = FNV0;
    for (size_t i = 0; i < nrc; i++) { h = fnvbuf(h, RC[i].name, strlen(RC[i].name)); h = fnvbuf(h, "\n", 1); }
    return h;
}

/* A gap pair as one varint: 0 none, 3k+1 an insertion of k+1, 3k+2 a
   deletion of k+1, 3 both or a negative one, two zigzags after. */
static void gapput(buf *b, int64_t dq, int64_t dt) {
    if (!dq && !dt) bvar(b, 0);
    else if (dq > 0 && !dt && dq < INT64_MAX / 3) bvar(b, 3 * (uint64_t)(dq - 1) + 1);
    else if (dt > 0 && !dq && dt < INT64_MAX / 3) bvar(b, 3 * (uint64_t)(dt - 1) + 2);
    else { bvar(b, 3); bvar(b, zz(dq)); bvar(b, zz(dt)); }
}
static void gapget(const uint8_t **p, const uint8_t *e, int64_t *dq, int64_t *dt) {
    uint64_t v = vgetb(p, e);
    *dq = *dt = 0;
    if (v == 3) { *dq = unzz(vgetb(p, e)); *dt = unzz(vgetb(p, e)); }
    else if (v % 3 == 1) *dq = (int64_t)(v / 3) + 1;
    else if (v % 3 == 2) *dt = (int64_t)(v / 3) + 1;
    else if (v) die("corrupt packed alignment: a bad gap code");
}

static void alnpack(const char *inpaf, const char *out) {
    loadpaf(inpaf);
    /* A block that starts a target run carries an absolute position: its
       query offset as a varint and its target position as four byte planes,
       since a random genome coordinate does not shrink as a varint. A block
       within a run carries the gap since the last: the query bases skipped
       and the target bases skipped, measured from the block's own end on the
       minus strand, so both are the indel sizes. One varint codes them,
       since one is nearly always 0: gapcode. */
    buf names = {0}, cnt = {0}, pos = {0}, len = {0}, seg_ = {0}, fq = {0}, pl[4] = {{0}, {0}, {0}, {0}};
    for (size_t i = 0; i < nqs; i++) {
        const qseg *Q = QS + i;
        if (!Q->n) continue;   /* no aligned base: nothing to pack */
        bput(&names, Q->name, strlen(Q->name)); bput(&names, "\n", 1);
        bvar(&cnt, Q->n);

        uint64_t pq = 0, pt = 0, ps = 0; uint32_t px = 0; int prev_rev = -1; uint64_t run = 0;
        for (size_t j = 0; j < Q->n; j++) {
            const seg *g = Q->s + j;
            if (!j || g->t != px || g->rev != prev_rev) {   /* absolute */
                if (g->tp > UINT32_MAX) die("a target position past 32 bits cannot be packed");
                bvar(&fq, g->q);
                for (int k = 0; k < 4; k++) { uint8_t b = (uint8_t)(g->tp >> (8 * k)); bput(&pl[k], &b, 1); }
            } else gapput(&pos, (int64_t)(g->q - pq), (int64_t)(g->rev ? ps - (g->tp + g->n) : g->tp - pt));
            if (g->t != px || g->rev != prev_rev) {
                if (j) bvar(&seg_, run);
                bvar(&seg_, (uint64_t)g->t * 2 + (g->rev ? 1 : 0));
                px = g->t; prev_rev = g->rev; run = 0;
            }
            run++;
            bvar(&len, g->n);
            pq = g->q + g->n; pt = g->tp + g->n; ps = g->tp;
        }
        bvar(&seg_, run);                       /* close the last run */
    }
    /* each stream gets the codec search */
    buf all = {0};
    buf *st[10] = { &names, &cnt, &pos, &len, &seg_, &fq, &pl[0], &pl[1], &pl[2], &pl[3] };
    static const char *SNAME[10] = { "names", "counts", "deltas", "lengths", "targets", "qstarts", "tpos0", "tpos1", "tpos2", "tpos3" };
    for (int k = 0; k < 10; k++) {
        uint8_t *best = NULL; size_t bn = 0; int bc = -1;
        if (st[k]->n) {
            if (k == 0) {
                uint8_t *t = NULL; size_t g = tokbest(st[k], &t);
                if (g) { best = t; bn = g; bc = ALNTOK; }
            }
            if (!best) {
            ZSTD_CCtx *zc = ZSTD_createCCtx();
            size_t zb = ZSTD_compressBound(st[k]->n);
            uint8_t *zo = xmalloc(zb);
            bn = ZSTD_compressCCtx(zc, zo, zb, st[k]->d, st[k]->n, 19);
            if (ZSTD_isError(bn)) die("packing the alignment failed");
            best = zo; bc = 0;
            ZSTD_freeCCtx(zc);
            /* read whole at every query: bzip3 decodes at 8 MB/s, so no */
            static const int C[] = { 0x100, 1, 65, 521, CODEC_ARITH, CODEC_ARITH | 1 };
            for (unsigned c = 0; c < sizeof C / sizeof *C; c++) {
                uint8_t *t = NULL; size_t g = paycomp(st[k]->d, st[k]->n, C[c], &t);
                if (g && g < bn) { free(best); best = t; bn = g; bc = C[c]; }
                else free(t);
            }
            }
        }
        fprintf(stderr, "  aln %-9s raw %10zu  coded %10zu  codec 0x%x\n",
                SNAME[k], st[k]->n, bn, (unsigned)bc);

        bvar(&all, bc < 0 ? 0 : bc == ALNTOK ? 2 : bc == 0 ? 1 : (uint64_t)bc + 3);
        bvar(&all, st[k]->n); bvar(&all, bn);
        if (bn) bput(&all, best, bn);
        free(best);
    }
    FILE *f = fopen(out, "wb");
    if (!f) die("cannot write the packed alignment");
    uint64_t rh = refhash();
    if (fwrite(ALNMAGIC, 1, ALNMAGLEN, f) != ALNMAGLEN || fwrite(&rh, 8, 1, f) != 1 ||
        (all.n && fwrite(all.d, 1, all.n, f) != all.n) || fclose(f)) die("packed alignment write failed");
    size_t nb = 0; for (size_t i = 0; i < nqs; i++) nb += QS[i].n;
    fprintf(stderr, "packed alignment: %zu queries, %zu blocks, %zu B (%.2f B/query)\n",
            nqs, nb, all.n + ALNMAGLEN + 8, nqs ? (double)(all.n + ALNMAGLEN + 8) / (double)nqs : 0.0);
    free(names.d); free(cnt.d); free(pos.d); free(len.d); free(seg_.d); free(fq.d); free(all.d);
    for (int k = 0; k < 4; k++) free(pl[k].d);
}

static int alnload(const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f) return 0;
    char mg[ALNMAGLEN];
    if (fread(mg, 1, ALNMAGLEN, f) != ALNMAGLEN || memcmp(mg, ALNMAGIC, ALNMAGLEN)) {
        if (!memcmp(mg, "KSTAL", 5)) die("this packed alignment is an older format: pack it again with this version");
        fclose(f); return 0;                       /* not packed: it is a PAF */
    }
    uint64_t rh;
    if (fread(&rh, 8, 1, f) != 1) die("truncated packed alignment");
    if (rh != refhash()) die("this packed alignment was made against another reference, or another contig order");
    if (fseek(f, 0, SEEK_END)) die("cannot seek the packed alignment");
    long fsz = ftell(f);
    if (fsz < ALNMAGLEN + 8) die("truncated packed alignment");
    size_t n = (size_t)fsz - ALNMAGLEN - 8;
    uint8_t *d = xmalloc(n ? n : 1);
    if (fseek(f, ALNMAGLEN + 8, SEEK_SET)) die("cannot seek the packed alignment");
    if (n && fread(d, 1, n, f) != n) die("short read on the packed alignment");
    if (ferror(f)) die("read error on the PAF");
    fclose(f);
    const uint8_t *p = d, *e = d + n;
    const uint8_t *sec[10]; size_t slen[10]; uint8_t *raw[10] = {0};
    for (int k = 0; k < 10; k++) {
        uint64_t cd = vgetb(&p, e), rl = vgetb(&p, e), cl2 = vgetb(&p, e);
        if (cl2 > (size_t)(e - p)) die("corrupt packed alignment: a stream overruns the file");
        if (!cd) { sec[k] = p; slen[k] = 0; continue; }
        if (cd == 2) {                                   /* tokenised names */
            size_t on = 0;
            raw[k] = (uint8_t *)tokd(p, (size_t)cl2, &on);
            sec[k] = raw[k]; slen[k] = on; p += cl2; continue;
        }
        raw[k] = xmalloc(rl ? rl : 1);
        if (cd == 1) {                                            /* zstd */
            size_t g = ZSTD_decompress(raw[k], (size_t)rl, p, (size_t)cl2);
            if (ZSTD_isError(g) || g != rl)
                die("corrupt packed alignment: a stream does not decompress");
        } else if (!paydecomp(p, (size_t)cl2, (int)(cd - 3), raw[k], (size_t)rl))
            die("corrupt packed alignment: a stream does not decompress");
        sec[k] = raw[k]; slen[k] = (size_t)rl; p += cl2;
    }
    const uint8_t *cp = sec[1], *cpe = sec[1] + slen[1];
    const uint8_t *pp = sec[2], *ppe = sec[2] + slen[2];
    const uint8_t *lp = sec[3], *lpe = sec[3] + slen[3];
    const uint8_t *sp = sec[4], *spe = sec[4] + slen[4];
    const uint8_t *fqp = sec[5], *fqe = sec[5] + slen[5];
    size_t nabs = slen[6], m = 0;
    for (int k = 7; k < 10; k++) if (slen[k] != nabs) die("corrupt packed alignment: the position planes differ in length");
    const char *nm = (const char *)sec[0], *nme = nm + slen[0];
    while (nm < nme) {
        const char *nl = memchr(nm, '\n', (size_t)(nme - nm));
        if (!nl) die("corrupt packed alignment: a name has no newline");
        char *name = xmalloc((size_t)(nl - nm) + 1);
        memcpy(name, nm, (size_t)(nl - nm)); name[nl - nm] = 0;
        qseg *Q = qseg_get(name, 1);
        free(name);
        if (Q->s) die("corrupt packed alignment: a query name repeats");
        uint64_t cn = vgetb(&cp, cpe);
        if (!cn || cn > (uint64_t)(lpe - lp)) die("corrupt packed alignment: more blocks than lengths");
        Q->s = xrealloc(Q->s, mul((size_t)cn, sizeof *Q->s)); Q->cap = cn; Q->n = 0;
        uint64_t pq = 0, pt = 0, ps = 0; uint32_t tx = 0; int rev = 0, fresh = 1;
        uint64_t left = 0;
        for (uint64_t j = 0; j < cn; j++) {
            uint64_t q, tp, ln = vgetb(&lp, lpe);
            if (!ln || ln > UINT32_MAX) die("corrupt packed alignment: a block of no length, or past 32 bits");
            if (!left) {                          /* next target/strand run */
                uint64_t xr = vgetb(&sp, spe);
                rev = (int)(xr & 1); tx = (uint32_t)(xr >> 1);
                if (tx >= nrc) die("corrupt packed alignment: a target is not in the reference");
                left = vgetb(&sp, spe);
                if (!left) die("corrupt packed alignment: an empty target run");
                fresh = 1;
            }
            if (fresh) {
                if (m >= nabs) die("corrupt packed alignment: more target runs than positions");
                q = vgetb(&fqp, fqe);
                tp = (uint64_t)sec[6][m] | (uint64_t)sec[7][m] << 8 | (uint64_t)sec[8][m] << 16 | (uint64_t)sec[9][m] << 24;
                m++; fresh = 0;
            }
            else {
                int64_t dq, dt; gapget(&pp, ppe, &dq, &dt);
                q = pq + (uint64_t)dq;
                tp = rev ? ps - ln - (uint64_t)dt : pt + (uint64_t)dt;
            }
            Q->s[Q->n++] = (seg){ q, tx, (uint32_t)ln, tp, rev };
            pq = q + ln; pt = tp + ln; ps = tp; left--;
        }
        if (left) die("corrupt packed alignment: a target run outlives its query");
        nm = nl + 1;
    }
    if (cp != cpe || pp != ppe || lp != lpe || sp != spe || fqp != fqe || m != nabs || p != e)
        die("corrupt packed alignment: bytes left over");
    free(d);
    for (int k = 0; k < 10; k++) free(raw[k]);
    size_t nb2 = 0; for (size_t i = 0; i < nqs; i++) nb2 += QS[i].n;
    fprintf(stderr, "packed alignment: %zu query contigs, %zu collinear blocks\n", nqs, nb2);
    seghash();
    return 1;
}

static void loadaln(const char *path) {
    if (!alnload(path)) loadpaf(path);
}

static void loadordinals(const char *path) {
    FILE *f = zfopen(path);
    if (!f) die("cannot open the file given to --aln-names");
    char *ln = NULL; size_t cap = 0; ssize_t n; size_t c = 0;
    while ((n = getline(&ln, &cap, f)) > 0) {
        chomp(ln, &n);
        if (!c && (ln[0] == '>' || ln[0] == '@')) die("--aln-names takes a plain list, one name a line");
        if (!n) die("--aln-names: an empty line: every line is one name");
        if (c == nord) { nord = nord ? nord * 2 : 1 << 16;
                         ORD = xrealloc(ORD, nord * sizeof *ORD); }
        ORD[c++] = xstrdup(ln);
    }
    if (ferror(f)) die("read error on a companion file");
    free(ln); fclose(f); nord = c;
    {   uint64_t h = FNV0;
        for (size_t i = 0; i < nord; i++) h = fnvbuf(h, ORD[i], strlen(ORD[i]) + 1);
        ALNHASH ^= h; }
    fprintf(stderr, "ordinal map: %zu names\n", nord);
}

static const qseg *qlookup(const char *nm) {
    if (ORD) {
        const char *p = nm;
        while (*p >= '0' && *p <= '9') p++;
        if (!*p && p != nm) {                     /* it is a bare ordinal */
            unsigned long i = strtoul(nm, NULL, 10);
            if (i >= nord) return NULL;
            nm = ORD[i];
        }
    }
    return qseg_get(nm, 0);
}

static __thread const qseg *CQ; static __thread size_t CQI;
static __thread const rcon *CR; static __thread size_t CRI;
static const seg *segat(const qseg *Q, uint64_t pos) {   /* last seg with q <= pos */
    if (!Q || !Q->n) return NULL;
    size_t lo;
    if (Q == CQ) {
        lo = CQI;
        int k = 0;
        while (lo < Q->n && Q->s[lo].q <= pos && k++ < 8) lo++;
        while (lo > 0 && Q->s[lo - 1].q > pos && k++ < 16) lo--;
        if (k >= 16 || (lo < Q->n && Q->s[lo].q <= pos)) lo = Q->n + 1;   /* give up */
    } else lo = Q->n + 1;
    if (lo == Q->n + 1) {
        size_t hi = Q->n; lo = 0;
        while (lo < hi) { size_t m = (lo + hi) / 2;
                          if (Q->s[m].q <= pos) lo = m + 1; else hi = m; }
    }
    CQ = Q; CQI = lo;
    if (!lo) return NULL;
    const seg *g = Q->s + lo - 1;
    return pos < g->q + g->n ? g : NULL;
}
static inline uint64_t rst(const rcon *rc, size_t i) { return rc->st32 ? rc->st32[i] : rc->st[i]; }
static void rcompact(rcon *rc) {   /* a whole contig under 2^32: starts to 4 bytes */
    if (rc->st32 || rc->end >= (1ull << 32)) return;
    uint32_t *w = xmalloc(mul(rc->n ? rc->n : 1, 4));
    for (size_t i = 0; i < rc->n; i++) w[i] = (uint32_t)rc->st[i];
    free(rc->st); rc->st = NULL; rc->st32 = w;
}
/* count of starts <= tpos, inside the rows [lo, hi) that are decoded and
   whose last run ends at wend. A whole contig is lo = 0, hi = n. */
static size_t runat(const rcon *rc, uint64_t tpos, size_t lo, size_t hi, uint64_t wend) {
    if (tpos >= wend) return 0;         /* past the last run: no answer, not the last run's */
    size_t a;
    if (rc == CR && CRI >= lo && CRI <= hi) {
        a = CRI;
        int k = 0;
        while (a < hi && rst(rc, a) <= tpos && k++ < 8) a++;
        while (a > lo && rst(rc, a - 1) > tpos && k++ < 16) a--;
        if (k >= 16 || (a < hi && rst(rc, a) <= tpos)) a = hi + 1;
    } else a = hi + 1;
    if (a == hi + 1) {
        size_t b = hi; a = lo;
        while (a < b) { size_t m = (a + b) / 2;
                        if (rst(rc, m) <= tpos) a = m + 1; else b = m; }
    }
    CR = rc; CRI = a;
    return a;
}
static void refwin(rcon *rc, uint64_t tpos, size_t *lo, size_t *hi, uint64_t *wend);
static int REFPERBLK;   /* one sequence asked for, no tail-encoded row: fill by block */
static uint32_t reflabel(const qseg *Q, uint64_t pos) {
    const seg *g = segat(Q, pos);
    if (!g) return UINT32_MAX;
    uint64_t tpos = g->rev ? g->tp + (g->q + g->n - 1 - pos) : g->tp + (pos - g->q);
    rcon *rc = RC + g->t;
    size_t lo = 0, hi; uint64_t wend;
    if (REFPERBLK) refwin(rc, tpos, &lo, &hi, &wend);
    else { reffill(rc); hi = rc->n; wend = rc->end; }
    size_t a = runat(rc, tpos, lo, hi, wend);
    return a ? rc->lab[a - 1] : UINT32_MAX;
}

static uint64_t refend(const qseg *Q, uint64_t pos) {
    const seg *g = segat(Q, pos);
    LASTREV = g ? g->rev : 0;
    if (!g) return 0;
    uint64_t tpos = g->rev ? g->tp + (g->q + g->n - 1 - pos) : g->tp + (pos - g->q);
    rcon *rc = RC + g->t;
    size_t lo = 0, hi; uint64_t wend;
    if (REFPERBLK) refwin(rc, tpos, &lo, &hi, &wend);
    else { reffill(rc); hi = rc->n; wend = rc->end; }
    size_t a = runat(rc, tpos, lo, hi, wend);
    if (!a) return 0;
    uint64_t tend = a < hi ? rst(rc, a) : wend;        /* target end of that run */
    if (g->rev) { uint64_t ts = rst(rc, a - 1);
                  if (ts < g->tp) return 0;
                  return g->q + g->n - (ts - g->tp); }
    if (tend < g->tp || tend - g->tp > g->n) return 0;
    return g->q + (tend - g->tp);
}

static uint64_t refend_at(const char *nm, uint64_t pos) {
    /* __thread like every other reader cursor: smooth refuses --aln tracks
       today, so no threaded path reaches this, but a shared cursor here
       would be a silent wrong-length race the day one does. */
    static __thread const qseg *cq = NULL; static __thread char cn[512];
    if (strcmp(nm, cn)) { cq = qlookup(nm); snprintf(cn, sizeof cn, "%s", nm); }
    return refend(cq, pos);
}

static uint32_t refrow(const qseg *Q, uint64_t st, uint64_t len) {
    if (len < 2) return reflabel(Q, st);
    uint32_t a = reflabel(Q, st + len / 4), b = reflabel(Q, st + len / 2);
    if (a == b) return a;
    return a == reflabel(Q, st + 3 * len / 4) ? a : b;
}
static uint32_t refrow_at(const char *nm, uint64_t st, uint64_t len) {
    static __thread const qseg *cq = NULL; static __thread char cn[512];
    if (strcmp(nm, cn)) { cq = qlookup(nm); snprintf(cn, sizeof cn, "%s", nm); }
    return refrow(cq, st, len);
}

typedef struct { char *key; uint32_t id; } fent;
static fent *ftab; static size_t fcap, fn;
static char **fset;                     /* feature set names, in table order */
static size_t nfset;
static char ***fid2lab;                 /* featureID -> label per set */
static size_t nfid;

static size_t fhash(const char *s) {
    size_t h = 1469598103934665603ULL;
    for (; *s; s++) { h ^= (unsigned char)*s; h *= 1099511628211ULL; }
    return h;
}
static void fput(const char *key, uint32_t id) {
    if (fn * 2 >= fcap) {
        size_t nc = fcap ? fcap * 2 : 1 << 16;
        fent *nt = xcalloc(nc, sizeof *nt);
        for (size_t i = 0; i < fcap; i++)
            if (ftab[i].key) {
                size_t j = fhash(ftab[i].key) & (nc - 1);
                while (nt[j].key) j = (j + 1) & (nc - 1);
                nt[j] = ftab[i];
            }
        free(ftab); ftab = nt; fcap = nc;
    }
    size_t j = fhash(key) & (fcap - 1);
    while (ftab[j].key) { if (!strcmp(ftab[j].key, key)) return; j = (j + 1) & (fcap - 1); }
    ftab[j].key = xstrdup(key); ftab[j].id = id; fn++;
}
static uint32_t fget(const char *key) {
    if (!fcap) return UINT32_MAX;
    size_t j = fhash(key) & (fcap - 1);
    while (ftab[j].key) { if (!strcmp(ftab[j].key, key)) return ftab[j].id;
                          j = (j + 1) & (fcap - 1); }
    return UINT32_MAX;
}

/* --db hot loop: */
#define DBMAXSET 32
typedef struct { char **k; uint32_t *v; size_t cap, n; } sintern;
static sintern SI[DBMAXSET];
static uint32_t *TUPK = NULL;      /* nfset ids per slot */
static uint32_t *TUPV = NULL;
static uint8_t  *TUPU = NULL;      /* slot used */
static size_t TUPCAP = 0;
static int DBFAST = 0;

static uint32_t si_get(sintern *t, const char *k, int add) {
    if (!t->cap) {
        if (!add) return UINT32_MAX;
        t->cap = 256; t->k = xcalloc(t->cap, sizeof *t->k);
        t->v = xcalloc(t->cap, sizeof *t->v);
    }
    size_t j = fhash(k) & (t->cap - 1);
    while (t->k[j]) { if (!strcmp(t->k[j], k)) return t->v[j];
                      j = (j + 1) & (t->cap - 1); }
    if (!add) return UINT32_MAX;
    if ((t->n + 1) * 4 > t->cap * 3) {           /* grow, rehash */
        size_t oc = t->cap; char **ok = t->k; uint32_t *ov = t->v;
        t->cap *= 2; t->k = xcalloc(t->cap, sizeof *t->k);
        t->v = xcalloc(t->cap, sizeof *t->v);
        for (size_t i = 0; i < oc; i++) if (ok[i]) {
            size_t m = fhash(ok[i]) & (t->cap - 1);
            while (t->k[m]) m = (m + 1) & (t->cap - 1);
            t->k[m] = ok[i]; t->v[m] = ov[i];
        }
        free(ok); free(ov);
        j = fhash(k) & (t->cap - 1);
        while (t->k[j]) j = (j + 1) & (t->cap - 1);
    }
    t->k[j] = xstrdup(k); t->v[j] = (uint32_t)t->n++;
    return t->v[j];
}
static uint64_t tuphash(const uint32_t *v, size_t n) {
    uint64_t h = 1469598103934665603ULL;
    for (size_t i = 0; i < n; i++) {
        h ^= v[i]; h *= 1099511628211ULL;
    }
    return h;
}
static void tup_put(const uint32_t *v, size_t n, uint32_t id) {
    size_t j = tuphash(v, n) & (TUPCAP - 1);
    while (TUPU[j]) {
        if (!memcmp(TUPK + j * n, v, n * sizeof *v)) return;
        j = (j + 1) & (TUPCAP - 1);
    }
    memcpy(TUPK + j * n, v, n * sizeof *v); TUPV[j] = id; TUPU[j] = 1;
}
static uint32_t tup_get(const uint32_t *v, size_t n) {
    if (!TUPCAP) return UINT32_MAX;
    for (size_t i = 0; i < n; i++) if (v[i] == UINT32_MAX) return UINT32_MAX;
    size_t j = tuphash(v, n) & (TUPCAP - 1);
    while (TUPU[j]) {
        if (!memcmp(TUPK + j * n, v, n * sizeof *v)) return TUPV[j];
        j = (j + 1) & (TUPCAP - 1);
    }
    return UINT32_MAX;
}

static uint64_t FEATHASH;
static void loadfeatures(const char *path) {
    FILE *f = zfopen(path);
    if (!f) die("cannot open the features table");
    char *ln = NULL; size_t cap = 0; ssize_t n;
    uint32_t maxid = 0; FEATHASH = FNV0;
    while ((n = getline(&ln, &cap, f)) > 0) {
        FEATHASH = fnvbuf(FEATHASH, ln, (size_t)n);
        chomp(ln, &n);
        if (!nfset) {                                  /* header */
            for (char *p = ln; (p = strchr(p, '\t')); p++) nfset++;
            if (!nfset) die("features.tsv: the first line must be the tab-separated header");
            fset = xcalloc(nfset, sizeof *fset);
            char *p = strchr(ln, '\t'); size_t k = 0;
            while (p) { char *q = strchr(++p, '\t'); if (q) *q = 0;
                        fset[k++] = xstrdup(p); p = q; }
            for (size_t i = 0; i < nfset; i++) {
                if (!fset[i][0] || strchr(fset[i], '/')) die("features.tsv: a set name is empty or has a slash");
                for (size_t j = 0; j < i; j++) if (!strcmp(fset[i], fset[j])) die("features.tsv names a set twice");
            }
            continue;
        }
        char *tab = strchr(ln, '\t');
        if (!tab) { if (n) die("features.tsv: a row without a tab"); continue; }
        *tab = 0;
        if (ln[0] == '0' && ln[1]) die("features.tsv: a featureID with a leading zero");
        uint32_t id = (uint32_t)num_arg(ln, (1u << 26) - 1, "features.tsv: a featureID is not a number below 2^26");
        *tab = '\t';
        if (!id) die("features.tsv defines featureID 0, which is reserved for novel");
        if (id > maxid) maxid = id;
        fput(tab + 1, id);                             /* the tuple, verbatim */
        if (id >= nfid) {
            size_t nn = (size_t)id + 1024;
            fid2lab = xrealloc(fid2lab, nn * sizeof *fid2lab);
            memset(fid2lab + nfid, 0, (nn - nfid) * sizeof *fid2lab);
            nfid = nn;
        }
        if (fid2lab[id]) die("features.tsv defines a featureID twice");
        fid2lab[id] = xcalloc(nfset, sizeof **fid2lab);
        char *p = tab; size_t k = 0;
        while (p && k < nfset) { char *q = strchr(++p, '\t'); if (q) *q = 0;
                                 fid2lab[id][k++] = xstrdup(p); p = q; }
        if (k < nfset) die("a features.tsv row has fewer labels than sets");
        if (p) die("a features.tsv row has more labels than sets");
    }
    if (ferror(f)) die("read error on a companion file");
    free(ln); fclose(f);
    if (!nfset) die("features table has no header");
    fprintf(stderr, "features: %zu ids, %zu sets\n", fn, nfset);
    if (nfset <= DBMAXSET) {                  /* build the integer-tuple index */
        TUPCAP = 1; while (TUPCAP < fn * 2) TUPCAP <<= 1;
        TUPK = xcalloc(TUPCAP * nfset, sizeof *TUPK);
        TUPV = xcalloc(TUPCAP, sizeof *TUPV);
        TUPU = xcalloc(TUPCAP, sizeof *TUPU);
        uint32_t v[DBMAXSET];
        for (uint32_t id = 0; id < nfid; id++) {
            if (!fid2lab[id]) continue;
            int ok = 1;
            for (size_t k = 0; k < nfset && ok; k++) {
                if (!fid2lab[id][k]) ok = 0;
                else v[k] = si_get(&SI[k], fid2lab[id][k], 1);
            }
            if (ok) tup_put(v, nfset, id);
        }
        DBFAST = 1;
        fprintf(stderr, "features: integer-tuple index on\n");
    }
}

/* --hier: the label vocabulary of each set, features.tsv column plus the
   hierarchy nodes of that set plus "novel", sorted. A tuple outside
   features.tsv gets 2^26 + its mixed-radix rank, the same in every file. */
static char ***VOC; static size_t *NVOC; static uint64_t *VSTRIDE, VPROD, HIERHASH;
static int CARTON, CARTUSED;
#define CARTBASE (1u << 26)
static int cmpstr(const void *a, const void *b) { return strcmp(*(char *const *)a, *(char *const *)b); }
static void vocadd(size_t k, const char *w, size_t n) {
    if (NVOC[k] % 1024 == 0) VOC[k] = xrealloc(VOC[k], (NVOC[k] + 1024) * sizeof **VOC);
    VOC[k][NVOC[k]++] = strndup(w, n);
}
static void loadvocab(const char *path) {
    FILE *f = zfopen(path);
    if (!f) die("cannot open the hierarchy file");
    VOC = xcalloc(nfset, sizeof *VOC); NVOC = xcalloc(nfset, sizeof *NVOC); VSTRIDE = xcalloc(nfset, sizeof *VSTRIDE);
    for (size_t k = 0; k < nfset; k++) {
        vocadd(k, "novel", 5);
        for (uint32_t id = 0; id < nfid; id++) if (fid2lab[id]) vocadd(k, fid2lab[id][k], strlen(fid2lab[id][k]));
    }
    char *line = NULL; size_t cap = 0; ssize_t len; HIERHASH = FEATHASH;
    while ((len = getline(&line, &cap, f)) > 0) {
        HIERHASH = fnvbuf(HIERHASH, line, (size_t)len);
        chomp(line, &len);
        if (!len) continue;
        char *t1 = memchr(line, '\t', (size_t)len); if (!t1) continue;
        char *t2 = memchr(t1 + 1, '\t', (size_t)(len - (t1 + 1 - line))); if (!t2) continue;
        if (memchr(t2 + 1, '\t', (size_t)(len - (t2 + 1 - line)))) die("hierarchy: a row has more than three columns");
        if (!strncmp(line, "feature_set\t", 12)) continue;
        for (size_t k = 0; k < nfset; k++)
            if ((size_t)(t1 - line) == strlen(fset[k]) && !memcmp(line, fset[k], strlen(fset[k]))) {
                vocadd(k, t1 + 1, (size_t)(t2 - t1 - 1)); vocadd(k, t2 + 1, (size_t)(len - (t2 + 1 - line))); }
    }
    if (ferror(f)) die("read error on a companion file");
    free(line); fclose(f);
    VPROD = 1;
    for (size_t k = nfset; k-- > 0;) {
        qsort(VOC[k], NVOC[k], sizeof **VOC, cmpstr);
        size_t u = 0;
        for (size_t i = 0; i < NVOC[k]; i++) if (!u || strcmp(VOC[k][u - 1], VOC[k][i])) VOC[k][u++] = VOC[k][i]; else free(VOC[k][i]);
        NVOC[k] = u; VSTRIDE[k] = VPROD;
        if (VPROD > (UINT32_MAX - CARTBASE) / u) die("the label vocabulary is too large for featureIDs");
        VPROD *= u;
    }
    CARTON = 1;
    fprintf(stderr, "hierarchy: %llu label tuples possible\n", (unsigned long long)VPROD);
}
static size_t vocidx(size_t k, const char *w) {   /* the rank of a label in its set, or SIZE_MAX */
    size_t lo = 0, hi = NVOC[k];
    while (lo < hi) { size_t m = (lo + hi) / 2; int c = strcmp(VOC[k][m], w); if (!c) return m; if (c < 0) lo = m + 1; else hi = m; }
    return SIZE_MAX;
}
static const char *flab(uint32_t id, size_t k) {   /* the label of set k for a featureID */
    if (!id) return "novel";
    if (id < CARTBASE) return fid2lab[id][k];
    return VOC[k][(id - CARTBASE) / VSTRIDE[k] % NVOC[k]];
}

/* --sort  Records are a set: */
typedef struct { uint64_t nam; uint32_t off, nrun, idx; } srec;  /* name in rbuf, runs in rrun */
typedef struct { uint32_t len; uint32_t lab; } srun;

static char *rbuf; static size_t rbufn, rbufc;      /* names, NUL separated */
static srun *rrun; static size_t rrunn, rrunc;      /* runs of every record */
static srec *rrec; static size_t rrecn, rrecc;
static buf natbuf; static uint64_t NATN;            /* names, input order */
static int PERMON;                    /* --sort keeps names: store a permutation */
static int NATORD = 1;
static uint8_t *PERMBLOB; static size_t PERMBLOBLEN; static int PERMBLOBTOK;

static void psort(srec *v, size_t n);

/* labels lead, then the run count, then lengths; 0 on a tie */
static int runcmp(const srun *p, uint32_t pn, const srun *q, uint32_t qn) {
    uint32_t n = pn < qn ? pn : qn;
    for (uint32_t i = 0; i < n; i++)
        if (p[i].lab != q[i].lab) return p[i].lab < q[i].lab ? -1 : 1;
    if (pn != qn) return pn < qn ? -1 : 1;
    for (uint32_t i = 0; i < n; i++)
        if (p[i].len != q[i].len) return p[i].len < q[i].len ? -1 : 1;
    return 0;
}
static int srec_cmp(const void *a, const void *b) {   /* a tie keeps input order */
    const srec *x = a, *y = b;
    int c = runcmp(rrun + x->off, x->nrun, rrun + y->off, y->nrun);
    return c ? c : (x->nam > y->nam) - (x->nam < y->nam);
}
/* one sorted record to the writer, under its rank when names are permuted */
static void sortemit(ktw *w, const char *nm, uint32_t idx, const srun *rs, uint32_t n) {
    char ib[24]; size_t nl;
    ALNNAME = nm;
    if (PERMON) { nl = (size_t)sprintf(ib, "%u", idx); nm = ib; } else nl = strlen(nm);
    uint64_t pos = 0;
    for (uint32_t k = 0; k < n; k++) {
        const char *lb = labtab[rs[k].lab];
        ktw_row(w, nm, nl, rs[k].len, lb, strlen(lb), pos, 0);
        pos += rs[k].len;
    }
}

static void psort(srec *v, size_t n) {
    int P = nthreads();
    if (P > 64) P = 64;
    while (P & (P - 1)) P &= P - 1;        /* the merge needs a power of two */
    if (n < 1u << 16 || P < 2) { qsort(v, n, sizeof *v, srec_cmp); return; }
    size_t bnd[65];
    for (int i = 0; i <= P; i++) bnd[i] = n * (size_t)i / (size_t)P;
#pragma omp parallel for schedule(dynamic, 1) num_threads(P)
    for (int i = 0; i < P; i++)
        qsort(v + bnd[i], bnd[i + 1] - bnd[i], sizeof *v, srec_cmp);
    srec *tmp = malloc(n * sizeof *tmp);
    if (!tmp) { qsort(v, n, sizeof *v, srec_cmp); return; }   /* fall back */
    srec *src = v, *dst = tmp;
    for (int w = 1; w < P; w *= 2) {
#pragma omp parallel for schedule(dynamic, 1) num_threads(P)
        for (int i = 0; i < P; i += 2 * w) {
            size_t a = bnd[i], m = bnd[i + w < P ? i + w : P], e = bnd[i + 2 * w < P ? i + 2 * w : P];
            size_t p = a, q = m, o = a;
            while (p < m && q < e)
                dst[o++] = srec_cmp(src + p, src + q) <= 0 ? src[p++] : src[q++];
            while (p < m) dst[o++] = src[p++];
            while (q < e) dst[o++] = src[q++];
        }
        srec *t = src; src = dst; dst = t;
    }
    if (src != v) memcpy(v, src, n * sizeof *v);
    free(tmp);
}

static void sortflush(ktw *w) {
    psort(rrec, rrecn);
    for (size_t i = 0; i < rrecn; i++)
        sortemit(w, rbuf + rrec[i].nam, rrec[i].idx, rrun + rrec[i].off, rrec[i].nrun);
    fprintf(stderr, "sorted %zu records by label tuple then lengths\n", rrecn);
    rrecn = 0; rrunn = 0; rbufn = 0;
}

typedef struct { FILE *f; uint32_t nrun, idx; char nm[512];
                 srun *rs; size_t cap; int eof; } srun_rd;

static FILE **SPILL;
static size_t nspill, cspill;

static void spillwrite(void) {
    if (!rrecn) return;
    psort(rrec, rrecn);
    if (nspill == cspill) {
        cspill = cspill ? cspill * 2 : 16;
        SPILL = xrealloc(SPILL, cspill * sizeof *SPILL);
    }
    const char *sd = TMPD; if (!sd || !*sd) sd = getenv("TMPDIR"); if (!sd || !*sd) sd = "/tmp";
    char *tp; if (asprintf(&tp, "%s/karyozip.sort.XXXXXX", sd) < 0) die("oom");
    int fd = mkstemp(tp);
    if (fd < 0) die("cannot create a sort temp file: check --tmp-dir");
    unlink(tp); free(tp);
    FILE *f = fdopen(fd, "w+b");
    if (!f) die("cannot create a sort temp file: check --tmp-dir");
    SPILL[nspill++] = f;
    for (size_t i = 0; i < rrecn; i++) {
        const char *nm = rbuf + rrec[i].nam;
        uint32_t nl = (uint32_t)strlen(nm), nr = rrec[i].nrun;
        if (fwrite(&nr, 4, 1, f) != 1 || fwrite(&nl, 4, 1, f) != 1 ||
            fwrite(&rrec[i].idx, 4, 1, f) != 1 ||
            fwrite(nm, 1, nl, f) != nl ||
            fwrite(rrun + rrec[i].off, sizeof(srun), nr, f) != nr)
            die("sort spill write failed: is --tmp-dir full?");
    }
    if (fflush(f)) die("sort spill write failed: is --tmp-dir full?");
    fprintf(stderr, "spill %zu: %zu records\n", nspill - 1, rrecn);
    rrecn = 0; rrunn = 0; rbufn = 0;
}

static int srd_next(srun_rd *r) {
    uint32_t nl;
    if (fread(&r->nrun, 4, 1, r->f) != 1) { r->eof = 1; return 0; }
    if (fread(&nl, 4, 1, r->f) != 1 || nl >= sizeof r->nm) die("sort temp is corrupt");
    if (fread(&r->idx, 4, 1, r->f) != 1) die("sort temp is corrupt");
    if (fread(r->nm, 1, nl, r->f) != nl) die("sort temp is corrupt");
    r->nm[nl] = 0;
    if (r->nrun > r->cap) { r->cap = r->nrun; r->rs = xrealloc(r->rs, r->cap * sizeof *r->rs); }
    if (fread(r->rs, sizeof(srun), r->nrun, r->f) != r->nrun) die("sort temp is corrupt");
    return 1;
}

static int srd_cmp(const srun_rd *a, const srun_rd *b) {
    int c = runcmp(a->rs, a->nrun, b->rs, b->nrun);
    return c ? c : (a->idx > b->idx) - (a->idx < b->idx);
}

static void spillmerge(ktw *w) {
    srun_rd *rd = xcalloc(nspill, sizeof *rd);
    for (size_t i = 0; i < nspill; i++) {
        rd[i].f = SPILL[i]; rewind(rd[i].f);
        srd_next(&rd[i]);
    }
    uint64_t out = 0;
    for (;;) {
        size_t best = nspill;
        for (size_t i = 0; i < nspill; i++) {
            if (rd[i].eof) continue;
            if (best == nspill || srd_cmp(&rd[i], &rd[best]) < 0) best = i;
        }
        if (best == nspill) break;
        srun_rd *r = &rd[best];
        sortemit(w, r->nm, r->idx, r->rs, r->nrun);
        out++;
        srd_next(r);
    }
    for (size_t i = 0; i < nspill; i++) {
        fclose(rd[i].f); free(rd[i].rs);
    }
    fprintf(stderr, "external sort: %llu records merged from %zu runs\n",
            (unsigned long long)out, nspill);
    free(rd); free(SPILL); SPILL = NULL; nspill = 0;
}

/* Both encode paths end the same way: spill what is held, merge, free. */
static void sortfinish(ktw *w) {
    if (SORTMEM) {
        if (rrecn) spillwrite();
        if (nspill) spillmerge(w);
    } else if (SORT && rrecn) sortflush(w);
    if (SORT) { free(rbuf); free(rrun); free(rrec); }
}

static void sortpush(const char *nm, size_t nlen, uint64_t st, uint64_t len,
                     const char *lab, size_t llen) {

    if (nlen >= 512) die("name too long");
    int newrec = st == 0;
    if (!newrec && rrecn) {
        const char *cur = rbuf + rrec[rrecn - 1].nam;
        if (strlen(cur) != nlen || memcmp(cur, nm, nlen)) newrec = 1;
    }
    if (newrec && st != 0) die("a sequence does not start at 0");

    if (!newrec && st != SORTBP) die("gap or overlap: start is not the previous end");
    SORTBP = newrec ? len : SORTBP + len;
    if (newrec) {
        if (rrecn == rrecc) { rrecc = rrecc ? rrecc * 2 : 1 << 16;
            rrec = xrealloc(rrec, rrecc * sizeof *rrec); }
        if (rbufn + nlen + 1 > rbufc) { rbufc = (rbufn + nlen + 1) * 2;
            rbuf = xrealloc(rbuf, rbufc); }
        rrec[rrecn].off = (uint32_t)rrunn;
        rrec[rrecn].nrun = 0;
        rrec[rrecn].nam = rbufn;
        memcpy(rbuf + rbufn, nm, nlen); rbuf[rbufn + nlen] = 0;
        rbufn += nlen + 1;
        if (NATN > UINT32_MAX) die("--sort takes at most 2^32 records");
        rrec[rrecn].idx = (uint32_t)NATN;
        if (PERMON) {
            char ob[24];
            int on = sprintf(ob, "%llu", (unsigned long long)NATN);
            if (NATORD && ((size_t)on != nlen || memcmp(ob, nm, nlen))) NATORD = 0;
            bput(&natbuf, nm, nlen); bput(&natbuf, "\n", 1);
        }
        NATN++; rrecn++;
    }
    uint64_t held = rbufn + rrunn * sizeof(srun) + rrecn * sizeof(srec);
    if (SORTMEM && held > SORTMEM && st == 0) {

        srec last = rrec[rrecn - 1];
        char nmb[512];
        size_t nl3 = strlen(rbuf + last.nam);
        if (nl3 >= sizeof nmb) die("name too long");
        memcpy(nmb, rbuf + last.nam, nl3 + 1);
        rrecn--;
        spillwrite();
        memcpy(rbuf, nmb, nl3 + 1); rbufn = nl3 + 1;
        rrec[0].nam = 0; rrec[0].off = 0; rrec[0].nrun = 0;
        rrec[0].idx = last.idx;          /* the rank travels with the record */
        rrecn = 1; rrunn = 0;
    }
    if (!rrecn) die("a sequence does not start at 0");
    if (rrunn == UINT32_MAX) die("--sort=MEM: more than 2^32 runs held; give a smaller MEM");
    if (rrunn == rrunc) { rrunc = rrunc ? rrunc * 2 : 1 << 20;
        rrun = xrealloc(rrun, rrunc * sizeof *rrun); }
    if (len > 0xffffffffu) die("run longer than 4 GB");
    rrun[rrunn].len = (uint32_t)len;
    rrun[rrunn].lab = (uint32_t)intern(lab, llen);
    rrunn++; rrec[rrecn - 1].nrun++;
}

static void encode(void) {
    ktw *w = ktw_open(stdout);

    typedef struct { const char *nm; uint32_t nlen, llen; const char *lab;
                     uint64_t st, len; int same; } prow;   /* same: name equals the row before */
    const size_t SLAB = 1 << 24;
    char *slab = xmalloc(SLAB + 1);
    int T = nthreads();
    if (T > 32) T = 32;
    prow **pr = xcalloc((size_t)T, sizeof *pr);
    size_t *prn = xcalloc((size_t)T, sizeof *prn), *prc = xcalloc((size_t)T, sizeof *prc);
    size_t have = 0; int eof = 0;
    FILE *IN = zfopen(NULL);            /* plain or .gz, decided by the bytes */
    if (!IN) die("read error on stdin");
    double _tp = nowsec();
    while (!eof || have) {
        if (!eof) {
            size_t got = fread(slab + have, 1, SLAB - have, IN);
            have += got;
            if (got == 0) { if (ferror(IN)) die("read error on stdin"); eof = 1; }
        }
        char *end = slab + have;
        char *last = have ? memrchr(slab, '\n', have) : NULL;
        if (!last) {
            if (eof && have) die("input does not end with a newline");
            if (!eof && have == SLAB) die("a BED line is longer than the read buffer");

            if (!eof) continue;
            break;
        }
        end = last + 1;                       /* parse whole lines only */
        size_t span2 = (size_t)(end - slab);

        size_t cb[33];
        cb[0] = 0;
        for (int t = 1; t < T; t++) {
            size_t g = span2 * (size_t)t / (size_t)T;
            while (g < span2 && slab[g] != '\n') g++;
            cb[t] = g < span2 ? g + 1 : span2;
        }
        cb[T] = span2;
        for (int t = 1; t <= T; t++) if (cb[t] < cb[t - 1]) cb[t] = cb[t - 1];
        const char *perr[33] = {0};   /* a worker records, main dies after the region */
        double _tpar = nowsec();
#pragma omp parallel for schedule(static) num_threads(T)
        for (int t = 0; t < T; t++) {
            size_t n = 0;
            const char *e = NULL;
            char *p = slab + cb[t], *e2 = slab + cb[t + 1], *nlp;
            while (p < e2 && (nlp = memchr(p, '\n', (size_t)(e2 - p)))) {
                *nlp = 0;
                ssize_t len = (ssize_t)(nlp - p);
                char *line = p;
                p = nlp + 1;
                if (memchr(line, 0, (size_t)len)) { e = "a NUL byte inside a row"; break; }
                char *t1 = memchr(line, '\t', (size_t)len); if (!t1) { e = "bad line"; break; }
                char *t2 = memchr(t1 + 1, '\t', (size_t)(len - (t1 + 1 - line))); if (!t2) { e = "bad line"; break; }
                char *t3 = memchr(t2 + 1, '\t', (size_t)(len - (t2 + 1 - line))); if (!t3) { e = "bad line"; break; }
                if (t1 + 1 == t2 || t2 + 1 == t3) { e = "empty start or end column"; break; }
                uint64_t st = 0, en = 0;
                st = bedcoord_e(t1 + 1, t2, "start is not a number", &e); if (e) break;
                en = bedcoord_e(t2 + 1, t3, "end is not a number", &e); if (e) break;
                if (en < st) { e = "end is before start"; break; }

                if (en == st) { e = "end equals start: a run must cover at least one base"; break; }
                if (t1 == line) { e = "empty sequence name"; break; }
                if (n == prc[t]) {
                    prc[t] = prc[t] ? prc[t] * 2 : 1 << 16;
                    pr[t] = realloc(pr[t], prc[t] * sizeof **pr);
                    if (!pr[t]) { e = "oom"; break; }
                }
                pr[t][n].same = n && pr[t][n - 1].nlen == (uint32_t)(t1 - line)
                                && !memcmp(pr[t][n - 1].nm, line, (size_t)(t1 - line));
                pr[t][n].nm = line; pr[t][n].nlen = (uint32_t)(t1 - line);
                pr[t][n].lab = t3 + 1; pr[t][n].llen = (uint32_t)(len - (t3 + 1 - line));
                pr[t][n].st = st; pr[t][n].len = en - st;
                n++;
            }
            prn[t] = n; perr[t] = e;
        }
        PH_ADD(5, _tpar);                       /* threaded: split the lines */
        /* Pipelining this loop, so the workers parse slab N+1 while one
           thread drains slab N, has a ceiling of parse:lines, 0.25 s against
           parse:append's 1.004 s. Append is the larger of the pair, so it
           stays on the critical path whatever overlaps it. 3% of encode, on
           the loop that carries every input-side refusal. Declined. */
        for (int t = 0; t < T; t++) if (perr[t]) die(perr[t]);
        double _tser = nowsec();
        for (int t = 0; t < T; t++)
            for (size_t k = 0; k < prn[t]; k++) {
                prow *q = &pr[t][k];
                const char *lab = q->lab; size_t llen = q->llen;
                if (ALNPAF && llen == 1 && *lab == '=') die("a label \"=\" is reserved under --aln");
                if (SORT) sortpush(q->nm, q->nlen, q->st, q->len, lab, llen);
                else ktw_row(w, q->nm, q->nlen, q->len, lab, llen, q->st, q->same);
            }
        PH_ADD(6, _tser);                       /* serial: append the rows */
        have -= span2;
        if (have) memmove(slab, slab + span2, have);
        if (!eof && have == SLAB) die("a BED line is longer than the read buffer");
    }
    for (int t = 0; t < T; t++) free(pr[t]);
    free(pr); free(prn); free(prc); free(slab); fclose(IN);
    sortfinish(w);
    PH_ADD(0, _tp);
    ktw_tail(w);
    ktw_pushread(w);
    ktw_flush(w);

    ktw_close(w);
}


/* bzip3 splits the input into blocks of its MB argument, so every size at
   or above the stream gives one block and the same bytes: the smallest
   such size stands for all of them, and the rest are searched for nothing */
/* hts_pack and its arith twin refuse above 16 distinct bytes and the encoder
   then clears the PACK bit, so such a candidate is its plain sibling */
static int packdup(int codec, const uint8_t *p, size_t n) {
    if (!(codec & 0x80) || (codec & CODEC_BZIP3)) return 0;
    uint8_t seen[256] = {0}; int d = 0;
    for (size_t i = 0; i < n; i++) if (!seen[p[i]]) { seen[p[i]] = 1; if (++d > 16) return 1; }
    return 0;
}

/* bzip3 blocks at least as big as the m bytes ranked give the same bytes
   there, so only one of them runs: the one right for the whole n, the
   smallest that covers it, else the largest listed. A 16 MB sample cannot
   tell 16 MB blocks from 64 MB ones, and the tie went to 16. */
static int bz3dup(int codec, size_t m, size_t n, const int *cand, size_t ncand) {
    if (!(codec & CODEC_BZIP3)) return 0;
    uint64_t mb = (uint64_t)(codec & 0xffff) << 20, keep = 0, top = 0;
    for (size_t j = 0; j < ncand; j++)
        if (cand[j] & CODEC_BZIP3) {
            uint64_t o = (uint64_t)(cand[j] & 0xffff) << 20;
            if (o >= n && (!keep || o < keep)) keep = o;
            if (o > top) top = o;
        }
    return mb >= m && mb != (keep ? keep : top);
}

/* a query decodes a whole block: bzip3 at 8 MB/s, lzma at 44, arith at
   2.5x rANS */
/* a codec id as a word, for the diagnostics */
static const char *cname(int c) {
    static __thread char b[32];
    if (!c) return "zstd";
    if (c & CODEC_LZMA) return "lzma";
    if (c & CODEC_BZIP3) snprintf(b, sizeof b, "bzip3-%dMB", c & 0xffff);
    else snprintf(b, sizeof b, "%s-o%d%s%s", c & CODEC_ARITH ? "arith" : "rans", c & 1, c & 0x80 ? "+pack" : "", c & 0x40 ? "+rle" : "");
    return b;
}

static int fastenough(int codec) {
    return (codec & (CODEC_BZIP3 | CODEC_ARITH | CODEC_LZMA)) == 0;
}

static void pickcodecs(int rawfd, const blkmeta *bm, size_t nraw) {

        int CANDBUF[KCAND_N]; size_t NC = candlist(CANDBUF);
        const int *CAND = CANDBUF;
        const size_t NCAND = NC;

        size_t rawtot = 0, maxstream = 0;
        for (size_t i2 = 0; i2 < nraw; i2++) { rawtot += bm[i2].len;
            if (bm[i2].lens > maxstream) maxstream = bm[i2].lens;
            if (bm[i2].len - bm[i2].lens > maxstream) maxstream = bm[i2].len - bm[i2].lens; }
        int big = maxstream >= 0x7fffffffu;   /* rANS and arith decode below INT_MAX */
        int seekable = nraw > 1 || (uint64_t)rawtot > BLOCK;
        /* one rule, three samples: -e low its codecs on 4 MB, the default all
           of them on DEFSAMPLE MB, -e high all of them on the whole file */
        size_t SCAP = FASTMODE ? (4u << 20) : ((size_t)DEFSAMPLE << 20);

        if (SLOW) {
            size_t want = 0;
            for (size_t i2 = 0; i2 < nraw; i2++) {
                size_t nl2 = bm[i2].lens, ni2 = bm[i2].len - nl2;
                if (nl2 > want) want = nl2;
                if (ni2 > want) want = ni2;
            }
            SCAP = want > (2048u << 20) ? (size_t)(2048u << 20) : want;
            if (SCAP < (4u << 20)) SCAP = 4u << 20;   /* never below the default */
        }

        size_t pick[3], np = 0;
        /* Sampling the first, middle and last block is why compression cannot
           start until parsing has finished: the choice needs the last block.
           Overlapping the two phases would be worth about 26% of encode, they
           are 2.24 s each at -t 8 and strictly sequential, but sampling the
           first three blocks instead changes which codec wins on every encode
           and changes what -e high means. Measured and declined. */
        if (nraw) { pick[np++] = 0;
                    if (nraw > 2) pick[np++] = nraw / 2;
                    if (nraw > 1) pick[np++] = nraw - 1; }
        size_t zi = 0, zl = 0, ci[300] = {0}, cl[300] = {0};

        size_t z22i = 0, z22l = 0;

        uint8_t *sbuf = NULL; size_t sbcap = 0;
        for (size_t j = 0; j < np; j++) {
            size_t need = bm[pick[j]].len;
            if (need > sbcap) { sbcap = need; sbuf = xrealloc(sbuf, sbcap ? sbcap : 1); }
            rawget(rawfd, (off_t)bm[pick[j]].off, sbuf, need);
            const uint8_t *src = sbuf;
            size_t nl = bm[pick[j]].lens, ni = bm[pick[j]].len - nl;
            const uint8_t *ids = src + nl;
            if (nl > SCAP) nl = SCAP;
            if (ni > SCAP) ni = SCAP;
            /* The two zstd baselines are two more tasks in the same sweep,
               so the search is parallel from its first block. */
            size_t zb[2][2] = {{0, 0}, {0, 0}};
            int z22 = !FASTMODE && !LEVELSET;
            int cthr = nthreads();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(cthr)
#endif
            for (unsigned k = 0; k < NCAND + 2; k++) {
                if (k >= NCAND) {
                    if (k == NCAND + 1 && !z22) continue;
                    int lv = k == NCAND ? LEVEL : 22;
                    ZSTD_CCtx *sc = ZSTD_createCCtx();
                    size_t cap = ZSTD_compressBound(ni > nl ? ni : nl) + 1024;
                    uint8_t *zo = xmalloc(cap);
                    if (ni) zb[k - NCAND][0] = ZSTD_compressCCtx(sc, zo, cap, ids, ni, lv);
                    if (nl) zb[k - NCAND][1] = ZSTD_compressCCtx(sc, zo, cap, src, nl, lv);
                    free(zo); ZSTD_freeCCtx(sc);
                    continue;
                }
                /* A candidate the selection below cannot take is not worth
                   compressing. On a seekable track that is both bzip3
                   entries, whose bz3_new alone asks for its 16 or 64 MB. */
                if (bz3dup(CAND[k], ni > nl ? ni : nl, maxstream, CAND, NCAND)) { ci[k] = cl[k] = (size_t)-1 / 2; continue; }
                if (!((!seekable || fastenough(CAND[k]))
                      && (!big || (CAND[k] & (CODEC_BZIP3 | CODEC_LZMA))))) {
                    ci[k] = cl[k] = (size_t)-1 / 2; continue;
                }
                /* the sample sizes; the winner round trips below */
                uint8_t *t = NULL; size_t g;
                if (packdup(CAND[k], ids, ni)) ci[k] = (size_t)-1 / 2;
                else if (ni && (g = paycomp(ids, ni, CAND[k], &t))) ci[k] += g;
                else ci[k] = (size_t)-1 / 2;
                free(t); t = NULL;
                if (packdup(CAND[k], src, nl)) cl[k] = (size_t)-1 / 2;
                else if (nl && (g = paycomp(src, nl, CAND[k], &t))) cl[k] += g;
                else cl[k] = (size_t)-1 / 2;
                free(t);
            }
            if (ZSTD_isError(zb[0][0]) || ZSTD_isError(zb[0][1])
                || ZSTD_isError(zb[1][0]) || ZSTD_isError(zb[1][1]))
                die("sample compress failed");
            zi += zb[0][0]; zl += zb[0][1];
            z22i += zb[1][0]; z22l += zb[1][1];
        }

        int wi = 0, wl = 0;
        {   /* every mode: the smallest that round trips on what it was ranked on */
            size_t *ti = xcalloc(NCAND + 2, sizeof *ti), *tl = xcalloc(NCAND + 2, sizeof *tl);
            /* a candidate lives per stream: one that does not round-trip
               the ids may still be the smallest for the lengths */
            int *live_i = xcalloc(NCAND + 2, sizeof *live_i), *live_l = xcalloc(NCAND + 2, sizeof *live_l);
            for (unsigned k = 0; k < NCAND; k++) {
                int ok = !big || (CAND[k] & (CODEC_BZIP3 | CODEC_LZMA));
                live_i[k] = ci[k] < (size_t)-1 / 4 && ok;
                live_l[k] = cl[k] < (size_t)-1 / 4 && ok;
            }
            live_i[NCAND] = live_l[NCAND] = 1; live_i[NCAND + 1] = live_l[NCAND + 1] = !LEVELSET && !FASTMODE;
            /* one block whose streams fit the sample: the sample already
               compressed exactly this, so its sizes are the file's */
            int reuse = !SLOW || (nraw == 1 && bm[0].lens <= SCAP && bm[0].len - bm[0].lens <= SCAP);
            if (reuse) {
                for (unsigned k = 0; k < NCAND; k++) { ti[k] = ci[k]; tl[k] = cl[k]; }
                ti[NCAND] = zi; tl[NCAND] = zl;             /* the zstd baselines */
                ti[NCAND + 1] = z22i; tl[NCAND + 1] = z22l;
            }
            for (size_t b2 = 0; !reuse && b2 < nraw; b2++) {
                if (bm[b2].len > sbcap) { sbcap = bm[b2].len; sbuf = xrealloc(sbuf, sbcap); }
                rawget(rawfd, (off_t)bm[b2].off, sbuf, bm[b2].len);
                const uint8_t *src = sbuf;
                size_t nl = bm[b2].lens, ni = bm[b2].len - nl;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads())
#endif
                for (unsigned k = 0; k < NCAND + 2; k++) {
                    if (!live_i[k] && !live_l[k]) continue;
                    uint8_t *t = NULL; size_t g;
                    if (k < NCAND) {
                        if (ni && live_i[k]) { g = paycomp(src + nl, ni, CAND[k], &t); ti[k] += g ? g : (size_t)-1 / 4; free(t); t = NULL; }
                        if (nl && live_l[k]) { g = paycomp(src, nl, CAND[k], &t); tl[k] += g ? g : (size_t)-1 / 4; free(t); }
                    } else {
                        int lv = k == NCAND ? LEVEL : 22;
                        if (ni) { t = malloc(ZSTD_compressBound(ni));
                                  g = t ? ZSTD_compress(t, ZSTD_compressBound(ni), src + nl, ni, lv) : (size_t)-1;
                                  ti[k] += (!t || ZSTD_isError(g)) ? (size_t)-1 / 4 : g; free(t); t = NULL; }
                        if (nl) { t = malloc(ZSTD_compressBound(nl));
                                  g = t ? ZSTD_compress(t, ZSTD_compressBound(nl), src, nl, lv) : (size_t)-1;
                                  tl[k] += (!t || ZSTD_isError(g)) ? (size_t)-1 / 4 : g; free(t); }
                    }
                }
            }

            unsigned qi, ql;
            for (;;) {
                qi = NCAND; ql = NCAND;
                for (unsigned k = 0; k < NCAND + 2; k++) {
                    if (live_i[k] && ti[k] < ti[qi]) qi = k;
                    if (live_l[k] && tl[k] < tl[ql]) ql = k;
                }
                int bad_i = 0, bad_l = 0;
                for (size_t j2 = 0; j2 < (SLOW ? nraw : np) && !(bad_i && bad_l); j2++) {
                    size_t b2 = SLOW ? j2 : pick[j2];
                    if (bm[b2].len > sbcap) { sbcap = bm[b2].len; sbuf = xrealloc(sbuf, sbcap); }
                    rawget(rawfd, (off_t)bm[b2].off, sbuf, bm[b2].len);
                    const uint8_t *src = sbuf;
                    size_t nl = bm[b2].lens, ni = bm[b2].len - nl;
                    const uint8_t *ids = src + nl;
                    if (!SLOW && nl > SCAP) nl = SCAP;   /* the sample it was ranked on */
                    if (!SLOW && ni > SCAP) ni = SCAP;
                    uint8_t *t = NULL, *v = NULL; size_t g;
                    if (qi < NCAND && ni && !bad_i) {
                        g = paycomp(ids, ni, CAND[qi], &t); v = g ? malloc(ni) : NULL;
                        if (!g || !v || paydecomp(t, g, CAND[qi], v, ni) != ni || memcmp(v, ids, ni)) bad_i = 1;
                        free(t); free(v); t = v = NULL; }
                    if (ql < NCAND && nl && !bad_l) {
                        g = paycomp(src, nl, CAND[ql], &t); v = g ? malloc(nl) : NULL;
                        if (!g || !v || paydecomp(t, g, CAND[ql], v, nl) != nl || memcmp(v, src, nl)) bad_l = 1;
                        free(t); free(v); }
                }
                if (!bad_i && !bad_l) break;
                if (bad_i) live_i[qi] = 0;
                if (bad_l) live_l[ql] = 0;
            }
            wi = qi < NCAND ? CAND[qi] : 0;
            wl = ql < NCAND ? CAND[ql] : 0;
            if (!wi) LEVEL_I = qi == NCAND + 1 ? 22 : LEVEL;
            if (!wl) LEVEL_L = ql == NCAND + 1 ? 22 : LEVEL;
            fprintf(stderr, "codec: %s ranked %zu candidates on %s\n", FASTMODE ? "effort low" : SLOW ? "effort high" : "effort medium",
                    (size_t)NCAND + 2, SLOW ? "the whole file" : "the sample");
            free(ti); free(tl); free(live_i); free(live_l);
        }
        free(sbuf);
        IDCODEC = wi; LENCODEC = wl;
        fprintf(stderr, "codec: len=0x%x id=0x%x from %zu sample blocks, %s\n",
                LENCODEC, IDCODEC, np, seekable ? "seekable" : "one block");
}

static size_t complens(ZSTD_CCtx *tc, const uint8_t *src, size_t n, uint8_t **out, size_t bound) {
    if (LENCODEC) return paycomp(src, n, LENCODEC, out);
    *out = malloc(bound + 64);
    if (!*out) return 0;
    size_t a = ZSTD_compressCCtx(tc, *out, bound + 64, src, n, LEVEL_L < 0 ? LEVEL : LEVEL_L);
    if (ZSTD_isError(a)) { free(*out); *out = NULL; return 0; }
    return a;
}

/* A respelt block is a different stream from the one the codec was chosen
   on, so it counts only at the size it comes back at. */
/* Every stream by every codec of the mode in one pool, LZMA and bzip3
   first since they take longest. Per stream the smallest that round trips
   on its sample wins, ties to the earlier of zstd, the table, zstd-22.
   m[s] bytes are ranked of a stream n[s] bytes whole. Out: cd[s], lv[s]
   the zstd level or -1, sz[s] 0 if none came back. */
typedef struct { int s, c, lv, r, cls; size_t z, m; } rjob;
/* the slowest first, so no long job starts last: LZMA and zstd 17 and up,
   then bzip3, then the rest; within a class the larger stream first */
static int rjcmp(const void *a, const void *b) {
    const rjob *x = a, *y = b;
    return x->cls != y->cls ? x->cls - y->cls : x->m != y->m ? (x->m < y->m) - (x->m > y->m)
         : x->s != y->s ? x->s - y->s : x->r - y->r;
}
static void rankall(int ns, const uint8_t *const *src, const size_t *m, const size_t *n,
                    int *cd, int *lv, size_t *sz) {
    int cb[KCAND_N]; size_t nc = candlist(cb), nj = 0;
    rjob *jb = xmalloc((size_t)ns * (nc + 2) * sizeof *jb);
    for (int t = 0; t < ns; t++) {
        cd[t] = 0; lv[t] = -1; sz[t] = 0;
        if (!m[t]) continue;
        jb[nj++] = (rjob){ t, 0, LEVEL, 0, LEVEL >= 17 ? 0 : 2, 0, m[t] };
        for (size_t k = 0; k < nc; k++) {
            if (packdup(cb[k], src[t], m[t]) || bz3dup(cb[k], m[t], n[t], cb, nc)) continue;
            if (n[t] >= 0x7fffffffu && !(cb[k] & (CODEC_BZIP3 | CODEC_LZMA))) continue;   /* rANS, arith: below INT_MAX */
            jb[nj++] = (rjob){ t, cb[k], -1, (int)k + 1, cb[k] & CODEC_LZMA ? 0 : cb[k] & CODEC_BZIP3 ? 1 : 2, 0, m[t] };
        }
        if (!FASTMODE && !LEVELSET) jb[nj++] = (rjob){ t, 0, 22, (int)nc + 1, 0, 0, m[t] };
    }
    qsort(jb, nj, sizeof *jb, rjcmp);
    double *jt = xcalloc(nj, sizeof *jt), tp = nowsec();   /* per job, under --profile */
    uint8_t **bo = xcalloc((size_t)ns, sizeof *bo); long *bj = xmalloc((size_t)ns * sizeof *bj);
    for (int t = 0; t < ns; t++) bj[t] = -1;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads())
#endif
    for (size_t j = 0; j < nj; j++) {
        const uint8_t *p = src[jb[j].s]; size_t k = m[jb[j].s]; uint8_t *o = NULL;
        double t0 = PROF ? nowsec() : 0;
        if (jb[j].c) {   /* the stream's best so far is kept for the read-back */
            jb[j].z = paycomp(p, k, jb[j].c, &o);
            int t = jb[j].s;
#ifdef _OPENMP
#pragma omp critical(rbest)
#endif
            if (jb[j].z && (bj[t] < 0 || jb[j].z < jb[bj[t]].z || (jb[j].z == jb[bj[t]].z && jb[j].r < jb[bj[t]].r))) {
                free(bo[t]); bo[t] = o; o = NULL; bj[t] = (long)j; }
            free(o);
        } else {
            size_t cap = ZSTD_compressBound(k); o = malloc(cap);
            size_t z = o ? ZSTD_compress(o, cap, p, k, jb[j].lv) : 0;
            jb[j].z = !o || ZSTD_isError(z) ? 0 : z; free(o);
        }
        if (PROF) jt[j] = nowsec() - t0;
    }
    if (PROF) fprintf(stderr, "  rank pool: %zu jobs, %.3f s\n", nj, nowsec() - tp);
    tp = nowsec();
    if (PROF)   /* the slowest jobs set the pool's critical path */
        for (int top = 0; top < 5; top++) {
            size_t q = 0; for (size_t j = 1; j < nj; j++) if (jt[j] > jt[q]) q = j;
            if (jt[q] <= 0) break;
            fprintf(stderr, "  rank job: stream %d %s on %zu B, %.3f s\n", jb[q].s, cname(jb[q].c ? jb[q].c : 0), m[jb[q].s], jt[q]);
            jt[q] = -1;
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads() < ns ? nthreads() : ns)
#endif
    for (int t = 0; t < ns; t++)
        for (;;) {
            size_t q = nj;
            for (size_t j = 0; j < nj; j++)
                if (jb[j].s == t && jb[j].z && (q == nj || jb[j].z < jb[q].z || (jb[j].z == jb[q].z && jb[j].r < jb[q].r))) q = j;
            if (q == nj) break;
            int ok = 1;
            if (jb[q].c) {                        /* zstd is trusted, the rest read back */
                uint8_t *o = NULL, *v = NULL;
                size_t z = (long)q == bj[t] ? jb[q].z : paycomp(src[t], m[t], jb[q].c, &o);
                if ((long)q == bj[t]) { o = bo[t]; bo[t] = NULL; }
                ok = z && (v = malloc(m[t])) && paydecomp(o, z, jb[q].c, v, m[t]) == m[t] && !memcmp(v, src[t], m[t]);
                free(o); free(v);
            }
            if (ok) { cd[t] = jb[q].c; lv[t] = jb[q].lv; sz[t] = jb[q].z; break; }
            jb[q].z = 0;                          /* it does not come back: the next */
        }
    if (PROF) fprintf(stderr, "  rank read-back: %.3f s\n", nowsec() - tp);
    for (int t = 0; t < ns; t++) free(bo[t]);
    free(jb); free(jt); free(bo); free(bj);
}

static size_t lensok(ZSTD_CCtx *tc, const uint8_t *src, size_t n) {
    uint8_t *c = NULL; size_t z = complens(tc, src, n, &c, ZSTD_compressBound(n));
    size_t ok = 0;
    if (z) {
        uint8_t *v = malloc(n ? n : 1); size_t g = 0;
        if (v) {
            if (LENCODEC) g = paydecomp(c, z, LENCODEC, v, n);
            else { g = ZSTD_decompress(v, n, c, z);
                   if (ZSTD_isError(g)) g = 0; }
            if (g == n && !memcmp(v, src, n)) ok = z;
        }
        free(v);
    }
    free(c);
    return ok;
}

static int u64cmp(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b; return x < y ? -1 : x > y; }
static void namesunique(const buf *names, size_t n) {
    uint64_t *hs = xmalloc(mul(n, 8)); size_t k = 0;
    const char *p = (const char *)names->d, *e = p + names->n;
    while (p < e && k < n) { const char *nl = memchr(p, '\n', (size_t)(e - p)); if (!nl) break;
                             hs[k++] = fnvbuf(FNV0, p, (size_t)(nl - p)); p = nl + 1; }
    qsort(hs, k, 8, u64cmp);
    size_t dup = 0; for (size_t i = 1; i < k; i++) if (hs[i] == hs[i - 1] && (i == 1 || hs[i] != hs[i - 2])) dup++;
    if (dup) {   /* the repeats: compare the names themselves */
        uint64_t *rep = xmalloc(mul(dup, 8)); size_t nr = 0;
        for (size_t i = 1; i < k; i++) if (hs[i] == hs[i - 1] && (i == 1 || hs[i] != hs[i - 2])) rep[nr++] = hs[i];
        buf *seen = xcalloc(nr, sizeof *seen);   /* the distinct names of each hash, NUL-joined */
        for (p = (const char *)names->d; p < e; ) {
            const char *nl = memchr(p, '\n', (size_t)(e - p)); if (!nl) break;
            uint64_t h = fnvbuf(FNV0, p, (size_t)(nl - p));
            uint64_t *f = bsearch(&h, rep, nr, 8, u64cmp);
            if (f) { size_t j = (size_t)(f - rep), L = (size_t)(nl - p); buf *b = &seen[j];
                     for (size_t o = 0; o < b->n; o += strlen((char *)b->d + o) + 1)
                         if (strlen((char *)b->d + o) == L && !memcmp(b->d + o, p, L)) {
                             fprintf(stderr, "karyozip: sequence %.*s appears twice\n", (int)L, p); die("a sequence name repeats"); }
                     bput(b, p, L); bput(b, "", 1); }
            p = nl + 1;
        }
        for (size_t j = 0; j < nr; j++) free(seen[j].d);
        free(seen); free(rep);
    }
    free(hs);
}

/* the index blob: chunk bounds first, then every chunk on its own thread */
static uint64_t idxbuild(const buf *ib, buf *ichunks, buf *itab) {
    uint64_t nichunk = 0;
    const uint8_t *q = ib->d, *qe = ib->d + ib->n;
    uint64_t cpend = 0, seen = 0, ncap = ib->n / ICHUNK + 2;
    const uint8_t **cb = xmalloc(ncap * sizeof *cb); uint64_t *cpe = xmalloc(ncap * 8);
    cb[0] = q; cpe[0] = 0;
    while (q < qe) {
        uint64_t nruns = vgetb(&q, qe);
        cpend += nruns;
        if (++seen % ICHUNK == 0 || q >= qe) { nichunk++; cb[nichunk] = q; cpe[nichunk] = cpend; }
    }
    uint8_t **cz = xcalloc(nichunk ? nichunk : 1, sizeof *cz);
    size_t *cl = xcalloc(nichunk ? nichunk : 1, sizeof *cl);
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads())
    for (uint64_t i = 0; i < nichunk; i++) {
        buf one = { (uint8_t *)cb[i], (size_t)(cb[i + 1] - cb[i]), 0 };
        cl[i] = zcl(&one, &cz[i], METALEVEL);
    }
    for (uint64_t i = 0; i < nichunk; i++) {
        bput(ichunks, cz[i], cl[i]); free(cz[i]);
        bvar(itab, cl[i]); bvar(itab, cpe[i]);
    }
    free(cz); free(cl); free(cb); free(cpe);
    return nichunk;
}

static void ktw_close(ktw *w) {
    buf names = w->names, sizes = w->sizes;
    size_t nidx = w->nidx;
    uint64_t nblk = w->nblk, rows = w->rows;
    FILE *tmp = w->tmp;
    const blkmeta *bm = w->bm; size_t nraw = w->nraw;

    double _ts = nowsec();
    if (!(nraw == 1 && w->ibx.n)) pickcodecs(w->rawfd, bm, nraw);   /* else the layout search ranks it all */
    PH_ADD(8, _ts);

    size_t maxraw = 0;
    for (size_t i = 0; i < nraw; i++) {
        if (bm[i].len > maxraw) maxraw = bm[i].len;
        for (int v = 0; ALNPAF && v < 2; v++)
            if (bm[i].alen[v] + (bm[i].len - bm[i].lens) > maxraw) maxraw = bm[i].alen[v] + (bm[i].len - bm[i].lens);
    }
    size_t obound = ZSTD_compressBound(maxraw ? maxraw : 1);

    /* Four spellings of one block, all measured, the smallest kept: the two
       frames as they are; the run count heading each read's lengths, which
       empties the index; the two frames interleaved into one, where a
       repeated read pattern is a single string; and both. One block only,
       since a split read has no count to head it. */
    uint8_t *altbuf = NULL; size_t altlen = 0, alt_lens = 0;
    uint8_t *altv[3] = {0}; size_t altvn[3] = {0}, altvl[3] = {0};
    double _tl = nowsec();
    if (nraw == 1 && w->ibx.n) {
        size_t nl0 = bm[0].lens, ni0 = bm[0].len - nl0;
        uint8_t *rb0 = xmalloc(bm[0].len);
        rawget(w->rawfd, (off_t)bm[0].off, rb0, bm[0].len);
        ZSTD_CCtx *mc = ZSTD_createCCtx();
        /* The layout is decided on the mode's sample, as the codec is: whole
           reads up to 4 MB of lengths -e low, DEFSAMPLE MB the default, all
           of them -e high. Each layout gets its own best codec there; only the
           winner is built whole, compressed once and read back. */
        size_t cap = SLOW ? SIZE_MAX : (size_t)(FASTMODE ? 4 : DEFSAMPLE) << 20, pl = 0, pi = 0, pb = 0;
        {   const uint8_t *p = rb0, *e = rb0 + nl0, *q = rb0 + nl0, *qe = q + ni0, *r = w->ibx.d, *re = r + w->ibx.n;
            while (r < re) {
                uint64_t k = vgetb(&r, re);
                for (uint64_t j = 0; j < k && p < e && q < qe; j++) { vgetb(&p, e); vgetb(&q, qe); }
                if ((size_t)(p - rb0) > cap && pl) break;           /* at least one read */
                pl = (size_t)(p - rb0); pi = (size_t)(q - rb0 - nl0); pb = (size_t)(r - w->ibx.d);
            } }
        const uint8_t *ls = rb0, *is = rb0 + nl0, *bs = w->ibx.d;
        buf ibs = { (uint8_t *)w->ibx.d, pb, pb }, ic = {0}, it = {0};
        idxbuild(&ibs, &ic, &it);
        size_t sz[4] = {0}, n4[4] = {0};
        uint8_t *lay[4] = {0};                    /* plain, heads, grouped, interleaved */
        n4[1] = splice_heads(ls, pl, bs, pb, &lay[1]);
        if (!FASTMODE) n4[2] = grp_pack(ls, pl, is, pi, bs, pb, &lay[2]);   /* a sort a read back */
        n4[3] = ilv_pack(ls, pl, is, pi, bs, pb, &lay[3]);
        /* the labels and the four length spellings, ranked in one pool; a
           spelling's whole size is its sample's scaled by the lengths' */
        const uint8_t *src[5] = { is, ls, lay[1], lay[2], lay[3] };
        size_t m5[5] = { pi, pl, n4[1], n4[2], n4[3] }, n5[5] = { ni0, nl0 }, r5[5];
        for (int k = 2; k < 5; k++) n5[k] = pl ? (size_t)((unsigned __int128)m5[k] * nl0 / pl) : m5[k];
        int c5[5], l5[5];
        rankall(5, src, m5, n5, c5, l5, r5);
        int idok = !pi || r5[0];
        if (r5[1] && idok) sz[0] = r5[1] + r5[0] + ic.n;      /* plain pays its index */
        if (r5[2] && idok) sz[1] = r5[2] + r5[0];
        if (r5[3] && idok) sz[2] = r5[3] + r5[0];
        if (r5[4]) sz[3] = r5[4];                             /* interleaved: no ids frame */
        int L = 0;
        for (int k = 1; k < 4; k++) if (sz[k] && (!sz[L] || sz[k] < sz[L])) L = k;
        for (int k = 1; k < 4; k++) free(lay[k]);
        IDCODEC = c5[0]; LEVEL_I = l5[0];
        LENCODEC = c5[L + 1]; LEVEL_L = l5[L + 1];
        size_t hn = 0; uint8_t *hb = NULL;
        if (L) hn = splice_heads(rb0, nl0, w->ibx.d, w->ibx.n, &hb);
        if (L == 1) { altbuf = hb; altlen = alt_lens = hn; hb = NULL; }
        if (L == 2) { altlen = alt_lens = grp_pack(rb0, nl0, rb0 + nl0, ni0, w->ibx.d, w->ibx.n, &altbuf); GRP = 1; }
        if (L == 3) { altlen = ilv_pack(rb0, nl0, rb0 + nl0, ni0, w->ibx.d, w->ibx.n, &altbuf); alt_lens = hn; ILV = 1; }
        free(hb);
        RUNHEAD = L > 0;
        if (L && !lensok(mc, altbuf, altlen)) LENCODEC = 0;  /* it does not come back: zstd does */
        static const char *LN[4] = {"plain", "headed", "grouped", "interleaved"};
        fprintf(stderr, "layout: %s, lengths %s", LN[L], cname(LENCODEC));
        fprintf(stderr, "%s%s, decided on %zu of %zu B of lengths\n", L == 3 ? "" : ", labels ", L == 3 ? "" : cname(IDCODEC), pl, nl0);
        /* --aln keeps three length spellings: respell every one the same way */
        if (ILV || RUNHEAD) {
            altv[0] = altbuf; altvn[0] = altlen; altvl[0] = alt_lens;
            for (int v = 1; ALNPAF && v < 3; v++) {
                size_t an = bm[0].alen[v - 1];
                if (!an) continue;
                const uint8_t *alt = w->rawalt[v - 1].d + bm[0].aoff[v - 1];
                uint8_t *ah = NULL;
                size_t ahn = GRP ? grp_pack(alt, an, rb0 + nl0, ni0, w->ibx.d, w->ibx.n, &ah)
                                 : splice_heads(alt, an, w->ibx.d, w->ibx.n, &ah);
                if (!ahn) { free(ah); continue; }
                if (!ILV) { altv[v] = ah; altvn[v] = ahn; altvl[v] = ahn; continue; }
                uint8_t *ai = NULL;
                size_t ain = ilv_pack(alt, an, rb0 + nl0, ni0, w->ibx.d, w->ibx.n, &ai);
                free(ah);
                if (ain) { altv[v] = ai; altvn[v] = ain; altvl[v] = ahn; } else free(ai);
            }
        }
        free(rb0); free(ic.d); free(it.d); ZSTD_freeCCtx(mc);
    } else {                                      /* blocks: no layout search, the codec search's pick */
        fprintf(stderr, "layout: plain, lengths %s", cname(LENCODEC));
        fprintf(stderr, ", labels %s, decided on the codec search's sample\n", cname(IDCODEC));
    }
    PH_ADD(9, _tl);
    if (RUNHEAD) { free(w->ibx.d); w->ibx.d = NULL; w->ibx.n = 0; }

    {
        /* Parallel: compress apart, write in order. */
        double _tc = nowsec();
        int nthr = nthreads();
        _Atomic int failed = 0, wfail = 0;
#ifdef _OPENMP
#pragma omp parallel num_threads(nthr)
#endif
        {
            ZSTD_CCtx *tc = ZSTD_createCCtx();
            uint8_t *o = malloc(2 * obound + 64);        /* one per thread */
            uint8_t *rb = malloc(maxraw ? maxraw : 1);   /* the block, read back */
            if (!tc || !o || !rb) failed = 1;
#ifdef _OPENMP
#pragma omp for ordered schedule(dynamic, 1)
#endif
            for (size_t i = 0; i < nraw; i++) {
                if (failed) continue;
                rawget(w->rawfd, (off_t)bm[i].off, rb, bm[i].len);
                const uint8_t *src = rb;
                size_t nl = bm[i].lens, ni = bm[i].len - nl;
                const uint8_t *idsrc = src + nl;   /* fixed: nl may change below */
                size_t a = 0, b2; uint64_t nosent = 0;
                uint8_t *fa = NULL, *fb = NULL;

                const uint8_t *vs[3] = { src, NULL, NULL }; size_t vn[3] = { nl, 0, 0 };
                for (int v = 0; ALNPAF && v < 2; v++) {
                    size_t an = bm[i].alen[v]; const uint8_t *alt = w->rawalt[v].d + bm[i].aoff[v];
                    if (!an || (an == nl && !memcmp(alt, src, an))) continue;   /* same bytes as every-sentinel */
                    if (v == 1 && vs[1] && an == vn[1] && !memcmp(alt, vs[1], an)) continue;   /* plus equals plain */
                    vs[v + 1] = alt; vn[v + 1] = an;
                }
                if (RUNHEAD || ILV)
                    for (int v = 0; v < 3; v++) {
                        if (altv[v]) { vs[v] = altv[v]; vn[v] = altvn[v]; }
                        else if (v) vs[v] = NULL;
                    }
                static const int ORDER[3] = { 1, 2, 0 };
                size_t best = 0; int bv = -1;
                for (int k = 0; k < 3; k++) {
                    int v = ORDER[k]; if (!vs[v]) continue;
                    uint8_t *f2 = NULL; size_t a2;
                    a2 = complens(tc, vs[v], vn[v], &f2, ZSTD_compressBound(vn[v]));
                    if (!a2) { failed = 1; break; }
                    if (bv < 0 || a2 < best) { free(fa); fa = f2; best = a2; bv = v; }
                    else free(f2);
                }
                if (failed) { free(fa); continue; }
                a = best; nl = ILV ? altvl[bv] : vn[bv]; nosent = (uint64_t)bv;
                if (ILV) { b2 = 0; fb = NULL; o[0] = 0; }
                else if (IDCODEC) { b2 = paycomp(idsrc, ni, IDCODEC, &fb);
                               if (!b2) { failed = 1; free(fa); continue; } }
                else
                { b2 = ZSTD_compressCCtx(tc, o, obound, idsrc, ni, LEVEL_I < 0 ? LEVEL : LEVEL_I);
                  if (ZSTD_isError(b2)) { failed = 1; free(fa); continue; } }
#ifdef _OPENMP
#pragma omp atomic
#endif
                BYT_LENS += a;
#ifdef _OPENMP
#pragma omp atomic
#endif
                BYT_LABS += b2;
                const uint8_t *ids = fb ? fb : o;
#ifdef _OPENMP
#pragma omp ordered
#endif
                {
                    if (fwrite(fa, 1, a, tmp) != a || fwrite(ids, 1, b2, tmp) != b2) wfail = 1;
                    bvar(&sizes, a); bvar(&sizes, b2);
                    bvar(&sizes, ni); bvar(&sizes, ALNPAF ? nl << 2 | nosent : nl);
                    bvar(&sizes, bm[i].rows); bvar(&sizes, bm[i].bp);
                    bvar(&sizes, fnvbuf(fnvbuf(FNV0, fa, a), ids, b2));
                }
                free(fa); free(fb);
            }
            free(o); free(rb);
            ZSTD_freeCCtx(tc);
        }
        if (wfail) die("tmp write");
        if (failed) die("block compression failed");
        RAWHELD = nraw ? bm[nraw - 1].off + bm[nraw - 1].len : 0;   /* staged, not held */
        PH_ADD(1, _tc);
    }

    /* index blob */
    buf ib = w->ibx;
    buf labs = {0};
    for (int i = 0; i < nlabtab; i++) { bput(&labs, labtab[i], strlen(labtab[i]));
                                        bput(&labs, "\n", 1); }
    uint8_t *nz = NULL, *lz = NULL;

    int ordinal = w->ordok && nidx > 0;   /* the writer tracked this */
    if (NAMEMODE != 1 && nidx > 1) namesunique(PERMON ? &natbuf : &names, PERMON ? NATN : nidx);

    int dropnames = NAMEMODE == 1;
    if (dropnames && ALNPAF && !ordinal) {
        fprintf(stderr, "--names drop: the names are not dense ordinals and --aln "
                        "needs them to reach the PAF, so they are stored as "
                        "integer deltas\n");
        dropnames = 0;
    }

    if (dropnames && PERMON) { dropnames = 0; ordinal = 0; }
    else if (dropnames) ordinal = 1;                   /* drop: store none */

    if (!ordinal && w->ordok && !names.n && nidx) {
        for (uint64_t i2 = 0; i2 < nidx; i2++) {
            char ob[24];
            int on = sprintf(ob, "%llu\n", (unsigned long long)i2);
            bput(&names, ob, (size_t)on);
        }
    }

    if (ordinal && NAMES_OUT && NAMEMODE != 1)
        die("--names drop:FILE: the names are already the record ordinals, so a "
            "file beside the track would be redundant: the track rebuilds them");

    if (PERMON && !NATORD && ordinal && NAMEMODE != 1) {
        free(names.d); names = natbuf; natbuf.n = 0; natbuf.d = NULL;
        ordinal = 0; PERMON = 0;
    }
    if (ordinal) { names.n = 0; NAMES_OUT = NULL; }

    size_t nzl = 0; int tokused = 0, intused = 0, fldused = 0;
    double _tn = nowsec();
    if (!ordinal) {
        /* whole names up to 4 MB -e low, DEFSAMPLE MB the default, all -e high */
        buf smp = { names.d, 0, 0 }; uint64_t smpn = 0;   /* a view, not a copy */
        size_t ncap = SLOW ? SIZE_MAX : (size_t)(FASTMODE ? 4 : DEFSAMPLE) << 20;
        { const char *q = (const char *)names.d, *e2 = q + names.n;
          while (q < e2 && (size_t)(q - (const char *)names.d) < ncap) {
              const char *nl = memchr(q, '\n', (size_t)(e2 - q));
              if (!nl) break;
              q = nl + 1; smpn++; }
          smp.n = smp.cap = (size_t)((const uint8_t *)q - names.d); }
        uint8_t *smpz = NULL; int pick;               /* rank on the sample, run the winner */
        size_t smpl = bestnames(&smp, smpn, &smpz, &pick);
        if (smp.n == names.n) {                       /* the sample was the whole
                                                         buffer: the winner is done */
            nz = smpz; nzl = smpl;
            intused = pick == 1; tokused = pick == 2; fldused = pick == 3;
        } else {
            free(smpz);
            if (pick == 1) { nzl = intnames(&names, nidx, &nz); intused = nzl != 0; }
            else if (pick == 2) { nzl = tokc(&names, &nz, TOKPICK); tokused = nzl != 0; }
            else if (pick == 3) { nzl = fldnames(&names, nidx, &nz); fldused = nzl != 0; }
        }
        if (!nzl) { nzl = zc(&names, &nz); tokused = intused = fldused = 0; }  /* refused */
        { static const char *NK[4] = { "zstd", "integers", "tokeniser", "fields" };
          int k = intused ? 1 : tokused ? 2 : fldused ? 3 : 0;
          fprintf(stderr, "names: %s", NK[k]);
          if (k == 2) fprintf(stderr, " chunk %u", TOKPICK);
          fprintf(stderr, ", decided on %zu of %zu B\n", smp.n, names.n); }
        /* The rows carry a rank per record. */
        if (PERMON && NAMEMODE == 1) { free(natbuf.d); natbuf.n = 0; natbuf.d = NULL; }
        else if (PERMON && !NATORD && natbuf.n) {
            char **nt = xmalloc((size_t)nidx * sizeof *nt);
            { char *q3 = (char *)natbuf.d; uint64_t k3 = 0;
              char *e3 = q3 + natbuf.n;
              while (k3 < nidx && q3 < e3) {
                  char *nl3 = memchr(q3, 10, (size_t)(e3 - q3));
                  if (!nl3) break;
                  nt[k3++] = q3; *nl3 = 0; q3 = nl3 + 1; }
              if (k3 != nidx) die("input-order names do not match the records"); }
            buf srt = {0};
            { const char *q4 = (const char *)names.d, *e4 = q4 + names.n;
              while (q4 < e4) {
                  const char *nl4 = memchr(q4, 10, (size_t)(e4 - q4));
                  if (!nl4) break;
                  char *e5; unsigned long long v = strtoull(q4, &e5, 10);
                  if (e5 != nl4 || e5 == q4 || v >= nidx) die("a stored rank is not a record number");
                  bput(&srt, nt[v], strlen(nt[v])); bput(&srt, "\n", 1);
                  q4 = nl4 + 1; } }
            /* nt[] NUL-terminated every name in place and the tokeniser
               chunks on newlines: without this the whole buffer is one
               chunk and tok3 refuses it, which cost 6.7% on 20547946 reads */
            for (uint64_t k5 = 0; k5 < nidx; k5++) nt[k5][strlen(nt[k5])] = '\n';
            uint8_t *bp = NULL; int bpick;
            size_t best = bestnames(&srt, nidx, &bp, &bpick);
            uint8_t *na = NULL, *nb2 = NULL;
            size_t p1 = zc(&natbuf, &na), p2 = tokbest(&natbuf, &nb2);
            size_t pbest = p1; uint8_t *pp = na; int ptok = 0;
            if (p2 && p2 < pbest) { pbest = p2; pp = nb2; ptok = 1; }
            if (best && best <= nzl + pbest) {     /* sorted names win: no rank */
                free(nz); nz = bp; nzl = best;
                tokused = bpick == 2; intused = bpick == 1; fldused = bpick == 3;
                free(na); free(nb2); free(natbuf.d); natbuf.n = 0; natbuf.d = NULL;
                PERMON = 0;
            } else {
                PERMBLOB = pp; PERMBLOBLEN = pbest; PERMBLOBTOK = ptok;
                if (pp != na) free(na);
                if (pp != nb2) free(nb2);
                free(bp);                          /* the sorted-name form loses */
            }
            free(srt.d); free(nt);
        }
    }
    PH_ADD(2, _tn);
    double _ti = nowsec();
    buf ichunks = {0}, itab = {0};
    uint64_t nichunk = idxbuild(&ib, &ichunks, &itab);
    PH_ADD(3, _ti);
    uint8_t *iz = NULL, *itz = NULL;
    size_t izl = ichunks.n, itl = zcl(&itab, &itz, METALEVEL);
    iz = ichunks.d;
    uint8_t *sz; size_t szl = zcl(&sizes, &sz, METALEVEL);
    size_t lzl = labs.n ? zc(&labs, &lz) : 0;

#define WR(buf, n) do { size_t _n = (n); \
    if (_n && fwrite((buf), 1, _n, w->out) != _n) die("write failed"); } while (0)

    uint8_t *npz = PERMBLOB; size_t npl = PERMON ? PERMBLOBLEN : 0;
    int permtok = PERMBLOBTOK;
    WR(MAGIC, 8);
    uint64_t fh = CARTUSED ? HIERHASH : FEATHASH;
    if (tokused + intused + fldused > 1) die("internal: names flagged more than one way");
    uint64_t h[16] = { nidx, (NAMES_OUT ? 0 : nzl) | (tokused ? TOKFLAG : 0) | (CARTUSED ? CARTFLAG : 0)
                             | (intused ? INTFLAG : 0) | (fldused ? FLDFLAG : 0) | (dropnames && !(w->ordok && nidx > 0) ? DROPPED : 0)
                             | (ALNPAF ? LENFLAG : 0) | (TAILON > 0 ? TAILFLAG : 0)
                             | (RUNHEAD ? RUNFLAG : 0) | (ILV ? ILVFLAG : 0) | (GRP ? GRPFLAG : 0)
                             | (ALNPAF ? ALNFLAG : 0) | (DBFEATS ? DBFLAG : 0)
                             | (NAMES_OUT ? SIDEFLAG : 0)
                             | (npl || (PERMON && NAMEMODE == 1) ? PERMFLAG : 0)
                             | (permtok ? PERMTOK : 0), izl, szl, rows,
                      (uint64_t)(uint32_t)LENCODEC | (DBFEATS ? (fh & 0xffffffffu) << 32 : 0),
                      lzl, (uint64_t)nlabtab, BLOCK,
                      (uint64_t)(uint32_t)IDCODEC | (DBFEATS ? (fh >> 32) << 32 : 0), nichunk, itl,
                      TAILPRED, (uint64_t)TAILSKIP | ((ALNPAF ? ALNHASH : ROWHASH) << 8), npl, 0 };

    {   uint64_t mh = fnvbuf(FNV0, h, 128); mh = fnvbuf(mh, lz, lzl);
        mh = fnvbuf(mh, nz, nzl);
        mh = fnvbuf(mh, iz, izl); mh = fnvbuf(mh, sz, szl); mh = fnvbuf(mh, itz, itl);
        h[15] = fnvbuf(mh, npz, npl);
    }
    if (fwrite(h, 8, 16, w->out) != 16) die("write failed");
    WR(lz, lzl);
    if (NAMES_OUT) {
        FILE *nf = fopen(NAMES_OUT, "wb");
        if (!nf) die("cannot open the names file");
        if (fwrite(nz, 1, nzl, nf) != nzl || fclose(nf)) die("names write failed");
    } else {
        WR(nz, nzl);
    }
    WR(iz, izl);
    WR(sz, szl);
    WR(itz, itl);
    WR(npz, npl);
    double _tw = nowsec();
    rewind(tmp);
    char cp[1 << 20]; size_t r;
    while ((r = fread(cp, 1, sizeof cp, tmp)) > 0) WR(cp, r);
    if (ferror(tmp)) die("read error on the temporary file");
    if (fflush(w->out)) die("write failed");
    PH_ADD(7, _tw);
    fclose(tmp);
    if (ALNPAF)
        fprintf(stderr, "aln: %llu of %llu rows matched the reference (%.1f%%)%s",
                (unsigned long long)ALNHIT, (unsigned long long)rows,
                rows ? 100.0 * (double)ALNHIT / (double)rows : 0.0, "\n");
    fprintf(stderr, "encoded rows=%llu reads=%zu blocks=%llu names=%zuB%s index=%zuB offs=%zuB "
            "lens=%lluB labs=%lluB raw=%zuB\n",
            (unsigned long long)rows, nidx, (unsigned long long)nblk,
            nzl, dropnames ? " (dropped, needs --names FILE)" : ordinal ? " (ordinal, not stored)" : tokused ? " (tokenised)" : intused ? " (int deltas)" : fldused ? " (field deltas)" : "", izl, szl,
            (unsigned long long)BYT_LENS, (unsigned long long)BYT_LABS, RAWHELD);

    free(itab.d); free(ichunks.d); free(altbuf);
    for (int v = 1; v < 3; v++) free(altv[v]);
    free(lz); free(nz); free(sz); free(itz); free(npz);
    free(names.d); free(sizes.d); free(labs.d);
    free(w->bm);
    for (int v = 0; v < 2; v++) { free(w->rawalt[v].d); free(w->lalt[v].d); }
    close(w->rawfd);
    free(w->lens.d); free(w->ids.d); free(w->ibx.d);
    free(w);
}

/* ------------------------------------------------------- the library API */
/* karyozip.h documents the contract. These are thin: the CLI and a caller
   that skips the BED go through exactly the same writer, so a track built
   row by row is byte-identical to one encoded from the equivalent BED. */
kt_writer *kt_open(FILE *out) {
    /* the encoder keeps global state: a second writer in one process starts
       from nothing, and the codec flags as first given. Not in ktw_open: the
       CLI interns the reference's labels before it opens the writer. */
    static int idc0 = -2, lenc0, lvl0, lvlset0, blkset0; static uint64_t blk0;
    if (idc0 == -2) { idc0 = IDCODEC; lenc0 = LENCODEC; lvl0 = LEVEL; lvlset0 = LEVELSET; blk0 = BLOCK; blkset0 = BLOCKSET; }
    IDCODEC = idc0; LENCODEC = lenc0; LEVEL = lvl0; LEVELSET = lvlset0; BLOCK = blk0; BLOCKSET = blkset0;
    for (int i = 0; i < nlabtab; i++) free(labtab[i]);
    free(labhash); labhash = NULL;
    nlabtab = 0; lastid = -1; labhash_init = 0; lcache_init = 0;
    TAILSEEN = 0; TAILON = -1; TAILPRED = 0; ALNHIT = 0; RUNHEAD = ILV = GRP = 0;
    return ktw_open(out);
}
void kt_block(kt_writer *w, uint64_t rows) { (void)w; BLOCK = rows ? rows : UINT64_MAX; BLOCKSET = 1; }
void kt_level(kt_writer *w, int level) {
    (void)w;
    if (level < 1 || level > 22) die("kt_level: a zstd level is 1 to 22");
    LEVEL = level; LEVELSET = 1;
}
void kt_row(kt_writer *w, const char *name, size_t nlen,
            uint64_t run_len, const char *label, size_t llen, uint64_t start) {
    /* The CLI gets these for free: its parser split on tabs and newlines, so
       a name or label can never hold one, and a zero-length run cannot reach
       the writer. A library caller can produce all of them, and a track with
       a newline in a label does not decode: the label table is newline
       separated, so it comes back as "the label table has more lines than
       its count". Refusing here is what keeps encode lossless. */
    if (!nlen) die("empty sequence name");
    if (!run_len) die("end equals start: a run must cover at least one base");
    if (start > UINT64_MAX - run_len) die("a coordinate does not fit in 64 bits");
    int bad = -1;                          /* one pass per string, not five memchrs */
    for (size_t i = 0; i < nlen; i++) {
        unsigned char c = (unsigned char)name[i];
        if (c == 0 || c == '\n' || c == '\t') { bad = c; break; }
    }
    for (size_t i = 0; bad < 0 && i < llen; i++) {
        unsigned char c = (unsigned char)label[i];
        if (c == 0 || c == '\n') { bad = c; break; }
    }
    if (bad == 0) die("a NUL byte inside a row");
    if (bad == '\n') die("a newline inside a name or label");
    if (bad == '\t') die("a tab inside a sequence name");
    ktw_row(w, name, nlen, run_len, label, llen, start, 0);
}
void kt_close(kt_writer *w) {
    ktw_tail(w);        /* the held row was the last of its sequence */
    ktw_pushread(w);
    ktw_flush(w);
    ktw_close(w);
}
const char *kt_version(void) { return KARYOZIP_VERSION; }

/* --------------------------------------------------------------- decoding */
/* one block: file offset, both frame sizes, both decoded sizes, rows, bp,
   its FNV, and which of the three length spellings it used */
typedef struct { uint64_t off, size, l, raw, lraw, rows, bp, hash; uint8_t nos; } blk;
typedef struct {
    FILE *f; uint64_t nreads, rows, data0, block;
    char **name; char **lab; uint64_t nlab;
    blk *b; uint64_t nblk;
    int lenfile;                  /* per block: plain lengths; file: --aln lengths */
    const uint8_t *lp0, *ip0;
    uint32_t *htab; uint64_t hmask;
    ZSTD_DCtx *d; uint8_t *blkbuf; size_t blkcap;
    uint64_t cached;
    uint64_t curblk, curoff;      /* sequential cursor in the cached block */
    uint8_t *nz_keep; size_t nzl_keep; int names_ready, names_tok, names_int, names_fld, dropped;
    uint8_t *ibuf; int idcodec, lencodec;
    const uint8_t *cp, *cpe;   /* cpe: end of the decoded index chunk */
    uint64_t ci, cpend; int cvalid; rec crec;
    off_t nz_off; uint64_t mhash, feathash; int cart;
    uint64_t fsz; int hasidx;   /* data size, and whether an index trailer follows it */
    uint8_t *pz; size_t pzl; int perm_tok, perm_ext;  /* names in input order */
    char *namebuf, *natbufp;
    int name_borrowed;         /* the table belongs to another reader */
    int blk_borrowed;          /* so do the one block and its run counts */
    uint8_t *ichunk; uint64_t nichunk, *coff, *cpendv, curchunk;
    uint64_t *sst, *sen; int *sfeat; uint8_t *snov; uint64_t scap;
    const uint8_t *payend;        /* the cached block's payload end, per reader */
    const uint8_t *curlp, *curip;
    uint64_t rowhash, tailpred, tailskip; int taildec;   /* for a track read as --ref */
    int runhead, ilv, grp, rfilled; uint32_t *rtab;   /* run counts lifted out of the lengths */
} R;

static void loadichunk(R *r, uint64_t c) {
    if (c >= r->nichunk) die("index chunk out of range");
    size_t cl = (size_t)(r->coff[c + 1] - r->coff[c]);
    const uint8_t *src = r->ichunk + r->coff[c];
    size_t rb; free(r->ibuf); r->ibuf = (uint8_t *)unz(src, cl, &rb, "index chunk");
    r->curchunk = c;
    r->cp = r->ibuf; r->cpe = r->ibuf + (size_t)rb;
    r->ci = c * ICHUNK; r->cpend = r->cpendv[c];

    {   const uint8_t *q = r->cp; uint64_t want = c + 1 < r->nichunk ? ICHUNK : r->nreads - c * ICHUNK, sum = 0, k = 0;
        for (; k < want && q < r->cpe; k++) { uint64_t v = vgetb(&q, r->cpe); if (v > UINT32_MAX) die("corrupt index: a run count past 32 bits"); sum += v; }
        if (k != want || q != r->cpe || r->cpendv[c] + sum != r->cpendv[c + 1])
            die("corrupt index: an index chunk does not match its table entry"); }
}

static void loadblk(R *r, uint64_t b);

/* the counts live in the lengths, so the first block load tabulates them */
static rec runidx(R *r, uint64_t i) {
    if (!r->rfilled) loadblk(r, 0);
    if (r->cvalid && r->ci == i) return r->crec;
    uint64_t off = 0, k0 = 0;
    if (r->cvalid && r->ci < i) { k0 = r->ci + 1; off = (uint64_t)r->crec.off + r->crec.nruns; }
    for (uint64_t k = k0; k < i; k++) off += r->rtab[k];
    if (off > UINT32_MAX) die("corrupt index: a sequence offset past 32 bits");
    rec out = { 0, (uint32_t)off, r->rtab[i] };
    r->ci = i; r->crec = out; r->cvalid = 1;
    if (out.off > r->b[0].rows || (uint64_t)out.off + out.nruns > r->rows)
        die("corrupt index: a sequence points outside the block table");
    return out;
}

static rec getidx(R *r, uint64_t i) {
    if (r->runhead) { if (i >= r->nreads) die("corrupt index: a sequence past the end");
                      return runidx(r, i); }
    if (r->cvalid && r->ci == i) return r->crec;
    uint64_t want = i / ICHUNK;
    if (!r->cvalid || i < r->ci || want != r->curchunk) {
        loadichunk(r, want);
        r->cvalid = 0;
    } else {
        r->ci++;                     /* cursor sits on ci, step past it */
    }
    rec out = r->cvalid ? r->crec : (rec){0, 0, 0};
    for (uint64_t k = r->ci; k <= i; k++) {

        uint64_t nruns = vgetb(&r->cp, r->cpe);
        if (r->cpend / (uint64_t)r->block > UINT32_MAX) die("corrupt index: a block number past 32 bits");
        out.blk = (uint32_t)(r->cpend / (uint64_t)r->block);
        out.off = (uint32_t)(r->cpend % (uint64_t)r->block);
        out.nruns = (uint32_t)nruns;
        r->cpend += nruns;
    }
    r->ci = i; r->crec = out; r->cvalid = 1;

    if (out.blk >= r->nblk || out.off > r->b[out.blk].rows
        || out.nruns > r->rows || (uint64_t)out.off + out.nruns > r->rows)
        die("corrupt index: a sequence points outside the block table");
    return out;
}

static R *ropen(const char *path) {
    R *r = xcalloc(1, sizeof *r);
    r->f = fopen(path, "rb"); if (!r->f) die("open");
    char m[8];
    if (fread(m, 1, 8, r->f) != 8) die("not a KST file: shorter than the magic");
    if (memcmp(m, MAGIC, 8)) {                 /* name what was found, and the version */
        char sh[9]; for (int i = 0; i < 8; i++) sh[i] = (m[i] >= 32 && m[i] < 127) ? m[i] : '.';
        sh[8] = 0;
        if (!memcmp(m, MAGIC, 3))
            fprintf(stderr, "karyozip: this is KST version %.5s, this build writes %.5s\n", m + 3, MAGIC + 3);
        else fprintf(stderr, "karyozip: file starts with \"%s\", not \"%.8s\"\n", sh, MAGIC);
        die("not a KST file this build can read");
    }
    uint64_t h[16]; if (fread(h, 8, 16, r->f) != 16) die("short header");
    uint64_t h1_raw = h[1];             /* flags included, for the hash */

    if (fseeko(r->f, 0, SEEK_END)) die("seek");
    uint64_t fsz = (uint64_t)ftello(r->f);
    /* An index trailer, if any, sits past the data: take it off fsz so every
       size check below still bounds the data alone. Magic last, so a torn
       append is simply not seen. */
    if (fsz >= 40 + 8 + 128) {           /* n, bucket bits, binding, payload, magic */
        uint64_t t[5];
        if (fseeko(r->f, (off_t)fsz - 40, SEEK_SET) || fread(t, 8, 5, r->f) != 5) die("seek");
        if (t[4] == IDXMAGIC && t[0] == h[0] && t[1] < 40 && t[3] <= fsz - 40 - 8 - 128)
            { fsz -= 40 + t[3]; r->hasidx = 1; }
    }
    r->fsz = fsz;
    ALLOCMAX = (size_t)1 << 38;
    if (fseeko(r->f, 8 + 128, SEEK_SET)) die("seek");
    uint64_t nml = h[1] & ~(TOKFLAG | INTFLAG | FLDFLAG | DROPPED | ALNFLAG | DBFLAG
                            | SIDEFLAG | LENFLAG | TAILFLAG | PERMFLAG
                            | PERMTOK | CARTFLAG | RUNFLAG | ILVFLAG | GRPFLAG);   /* h[1] also holds flags */

    if (nml > fsz || h[2] > fsz || h[3] > fsz || h[6] > fsz || h[11] > fsz
        || h[14] > fsz || h[11] + h[14] > fsz
        || 8 + 128 + nml + h[2] + h[3] + h[6] + h[11] + h[14] > fsz)
        die("corrupt header: the stream sizes do not fit the file");
    r->block = h[8] ? h[8] : 512;
    r->nreads = h[0]; r->rows = h[4];
    size_t lzl = (size_t)h[6]; r->nlab = h[7]; r->lab = NULL;
    uint8_t *lzraw = NULL;              /* kept until the metadata hash */
    if (lzl) {
        uint8_t *lz = xmalloc(lzl);
        if (fread(lz, 1, lzl, r->f) != lzl) die("short labels");
        lzraw = lz;   /* kept for the hash, freed there */
        size_t lb; char *lbuf = unz(lz, lzl, &lb, "label table");

        if (r->nlab > lb) die("label count exceeds the label table");
        if (!r->nlab) die("corrupt file: a label table with no labels");
        r->lab = xmalloc(mul(r->nlab, sizeof(char *)));
        char *q = lbuf;
        for (uint64_t i = 0; i < r->nlab; i++) {
            r->lab[i] = q; char *nl = strchr(q, '\n'); if (!nl) die("bad label table");
            *nl = 0; q = nl + 1; }
        if (q != lbuf + lb) die("corrupt file: the label table has more lines than its count");
    } else if (r->nlab) die("corrupt file: a label count with no label table");

    int names_tok = (h[1] & TOKFLAG) != 0, names_int = (h[1] & INTFLAG) != 0;
    int names_fld = (h[1] & FLDFLAG) != 0;
    if (names_tok + names_int + names_fld > 1) die("corrupt file: names flagged more than one way");
    int was_aln = (h[1] & ALNFLAG) != 0, was_db = (h[1] & DBFLAG) != 0;
    int was_side = (h[1] & SIDEFLAG) != 0, was_len = (h[1] & LENFLAG) != 0;
    if (NAMESARG) { if (was_side) NAMESIN = NAMESARG; else NAMEFILE = NAMESARG; }
    int was_perm = (h[1] & PERMFLAG) != 0, perm_tok = (h[1] & PERMTOK) != 0;

    if (!RAW && !INDEXING) {   /* index reads names only: no companion to bind */
        if (was_aln && !ALNPAF)
            die("this track was written with --aln: decode it with --aln PAF --ref R, "
                "or --raw to see what is stored");
        if (!was_aln && ALNPAF)
            die("this track was not written with --aln");
        if (was_db && !DBFEATS)
            die("this track holds featureIDs: decode it with --db features.tsv --out-dir DIR, "
                "or --raw to see what is stored");
        if (!was_db && DBFEATS)
            die("this track does not hold featureIDs");
        if (was_len && !ALNPAF)
            die("this track leaves out every length the reference implies: "
                "decode it with --aln PAF --ref R");
        if (was_side && !NAMESIN) {
            char dflt[4096]; snprintf(dflt, sizeof dflt, "%s.names", path);
            FILE *tf = fopen(dflt, "rb");
            if (!tf) die("the names of this track are in a sidecar: give it with "
                         "--names FILE, or --raw to see the ordinals");
            fclose(tf);                 /* <track>.names is beside it: use it */
        }
    }
    int dropped = (h[1] & DROPPED) != 0;

    r->lenfile = LENDEC = was_len && !RAW;

    TAILDEC = (h[1] & TAILFLAG) != 0;
    int runhead = (h[1] & RUNFLAG) != 0, ilv = (h[1] & ILVFLAG) != 0;
    TAILPRED_R = h[12]; TAILSKIP_R = h[13] & 0xff;
    if (TAILDEC && TAILSKIP_R != TAILSKIP) die("corrupt file: tail skip is not the format constant");
    r->taildec = TAILDEC; r->tailpred = h[12]; r->tailskip = TAILSKIP_R; r->rowhash = was_aln ? 0 : h[13] >> 8;
    if (was_aln && ALNPAF && (h[13] >> 8) != (ALNHASH & ((1ull << 56) - 1)))
        die("this PAF, reference or name map is not the one the track was encoded with");
    r->cart = (h[1] & CARTFLAG) != 0;
    r->runhead = runhead; r->ilv = ilv; r->grp = (h[1] & GRPFLAG) != 0;
    if (r->grp && (ilv || !runhead)) die("corrupt file: grouped lengths without their run counts");
    if (runhead) r->rtab = xmalloc(mul(h[0] ? h[0] : 1, sizeof *r->rtab));
    if (r->cart && !was_db) die("corrupt file: vocabulary featureIDs without --db");
    h[1] &= ~(TOKFLAG | INTFLAG | DROPPED | ALNFLAG | DBFLAG | SIDEFLAG | LENFLAG
              | TAILFLAG | PERMFLAG | PERMTOK | CARTFLAG | RUNFLAG | ILVFLAG | FLDFLAG | GRPFLAG);
    uint8_t *iz = xmalloc(h[2]), *sz = xmalloc(h[3]);
    uint8_t *nz = NULL;
    r->nz_off = ftello(r->f);
    if (fseeko(r->f, (off_t)h[1], SEEK_CUR)) die("seek past names");
    if (fread(iz, 1, h[2], r->f) != h[2]) die("short index");
    if (fread(sz, 1, h[3], r->f) != h[3]) die("short sizes");
    size_t nzl_eff = (size_t)h[1]; uint8_t *sidebuf = NULL;
    if (!nzl_eff && was_side) {
        /* the sidecar: whatever --names drop:FILE named, else <track>.names */
        char side[4096];
        if (NAMESIN) snprintf(side, sizeof side, "%s", NAMESIN);
        else snprintf(side, sizeof side, "%s.names", path);
        FILE *nf = fopen(side, "rb");
        if (nf) {
            fseeko(nf, 0, SEEK_END); nzl_eff = (size_t)ftello(nf); rewind(nf);
            sidebuf = xmalloc(nzl_eff ? nzl_eff : 1);
            if (fread(sidebuf, 1, nzl_eff, nf) != nzl_eff) die("short names file");
            fclose(nf);
        } else if (!RAW) {

            fprintf(stderr, "cannot open the names sidecar %s\n", side);
            die("this track keeps its names in a sidecar and it is missing");
        } else fprintf(stderr, "sidecar missing: names are ordinals and the metadata hash is unchecked\n");
    }
    if (!(was_side && !sidebuf && RAW))
    {
        uint64_t want = h[15], mh = FNV0, hz[16];
        memcpy(hz, h, 128); hz[15] = 0;

        hz[1] = h1_raw;
        mh = fnvbuf(mh, hz, 128); mh = fnvbuf(mh, lzraw, lzl);
        if (h[1]) { off_t keep = ftello(r->f);
            uint8_t *nb = xmalloc(h[1]);
            if (fseeko(r->f, (off_t)r->nz_off, SEEK_SET) || fread(nb, 1, h[1], r->f) != h[1])
                die("short names");
            mh = fnvbuf(mh, nb, h[1]); free(nb);
            if (fseeko(r->f, keep, SEEK_SET)) die("seek"); }
        else if (sidebuf) mh = fnvbuf(mh, sidebuf, nzl_eff);
        mh = fnvbuf(mh, iz, h[2]); mh = fnvbuf(mh, sz, h[3]);
        {   off_t keep = ftello(r->f);
            uint8_t *tb = xmalloc(h[11] + h[14]);
            if (fread(tb, 1, h[11] + h[14], r->f) != h[11] + h[14]) die("short label table");
            mh = fnvbuf(mh, tb, h[11] + h[14]); free(tb);
            if (fseeko(r->f, keep, SEEK_SET)) die("seek"); }
        free(lzraw);
        if (mh != want) {
            if (sidebuf) die("this names file does not belong to this track: give the "
                             "file --names drop:FILE wrote, not a FASTA or a list");
            die("corrupt file: the metadata does not match its hash");
        }
    }
    r->mhash = h[15];
    r->feathash = (h[5] >> 32) | ((h[9] >> 32) << 32);
    if (was_perm) {                     /* it sits between itz and the rows */
        r->pzl = (size_t)h[14]; r->perm_tok = perm_tok;
        r->perm_ext = r->pzl == 0;          /* the text is in the FASTA */
        r->pz = r->pzl ? xmalloc(r->pzl) : NULL;
        off_t here = ftello(r->f);
        if (fseeko(r->f, (off_t)h[11], SEEK_CUR)) die("seek past the label table");
        if (r->pzl && fread(r->pz, 1, r->pzl, r->f) != r->pzl)
            die("short input-order names");
        if (fseeko(r->f, here, SEEK_SET)) die("seek back");
    }
    r->data0 = 8 + 128 + lzl + h[1] + h[2] + h[3] + h[11] + r->pzl;

    r->nz_keep = sidebuf ? sidebuf : nz; r->nzl_keep = nzl_eff; r->names_ready = 0; r->names_tok = names_tok; r->names_int = names_int; r->names_fld = names_fld; r->dropped = dropped;
    r->name = NULL;                    /* 160 MB, allocated by ensurenames */

    /* Walk varints with a cursor, allocate nothing. */

    r->nichunk = h[10];
    r->ichunk = iz;                      /* the concatenated chunk frames */
    r->curchunk = UINT64_MAX;
    {   uint8_t *itz = xmalloc(h[11]);
        if (fread(itz, 1, h[11], r->f) != h[11]) die("short index table");
        size_t tb; uint8_t *tbuf = (uint8_t *)unz(itz, h[11], &tb, "index table");
        free(itz);
        if (r->nichunk == SIZE_MAX) die("corrupt file: chunk count");
        r->coff = xmalloc(mul(r->nichunk + 1, 8));
        r->cpendv = xmalloc(mul(r->nichunk + 1, 8));
        const uint8_t *q = tbuf, *qe = tbuf + tb;
        uint64_t acc = 0;
        for (uint64_t c = 0; c < r->nichunk; c++) {
            uint64_t cl = vgetb(&q, qe);
            r->cpendv[c] = vgetb(&q, qe);
            r->coff[c] = acc; acc += cl;
        }
        r->coff[r->nichunk] = acc;
        if (acc != h[2]) die("index chunk table does not match the index");
        if (!r->runhead && r->nichunk != (r->nreads + ICHUNK - 1) / ICHUNK) die("corrupt file: index chunks do not match the sequence count");
        r->cpendv[r->nichunk] = r->rows;
        for (uint64_t c = 0; c < r->nichunk; c++)
            if (r->cpendv[c] > r->cpendv[c + 1] || (c == 0 && r->cpendv[0])) die("corrupt file: index chunk starts are not in order");
        free(tbuf);
        r->ibuf = NULL;
    }
    r->lencodec = (int)(uint32_t)h[5]; r->idcodec = (int)(uint32_t)h[9];
    if (nosz(r->lencodec) || nosz(r->idcodec)) die("corrupt file: a NOSZ codec cannot be size-checked");
    if ((uint32_t)h[5] >> 31 || (uint32_t)h[9] >> 31) die("corrupt file: codec id out of range");
    r->ci = 0; r->cpend = 0; r->cvalid = 0;

    size_t sb; uint8_t *sbuf = (uint8_t *)unz(sz, h[3], &sb, "block sizes");
    { const uint8_t *p = sbuf, *end = sbuf + sb; uint64_t n = 0, cap = 1 << 16;
      r->b = xmalloc(mul(cap, sizeof *r->b));
      uint64_t acc = 0;
      while (p < end) { if (n == cap) { cap *= 2; r->b = xrealloc(r->b, mul(cap, sizeof *r->b)); }
        uint64_t a = vgetb(&p, end), b2 = vgetb(&p, end),
                 c2 = vgetb(&p, end), d2 = vgetb(&p, end), e2 = vgetb(&p, end),
                 f2 = vgetb(&p, end), g2 = vgetb(&p, end);
        if (was_len) { r->b[n].nos = d2 & 3; d2 >>= 2; } else r->b[n].nos = 0;

        if (a > fsz || b2 > fsz || a + b2 > fsz || acc + a + b2 > fsz)
            die("corrupt index: a block runs past the end of the file");
        if (e2 > r->rows) die("corrupt index: a block claims more rows than the file has");
        if (e2 && (!c2 || !d2)) die("corrupt index: a block with rows but an empty stream");
        if (c2 > r->rows * 16 + 64 || d2 > r->rows * 16 + 64)
            die("corrupt index: a block claims an impossible decoded size");
        r->b[n] = (blk){ acc, a + b2, a, c2, d2, e2, f2, g2, r->b[n].nos };
        acc += a + b2; n++; }
      r->nblk = n;
      if (r->runhead && n != 1)
          die("corrupt file: run counts head the lengths but the track has more than one block");
      uint64_t tot = 0;
      for (uint64_t i = 0; i < n; i++) {
          tot += r->b[i].rows;
          if (i + 1 < n ? r->b[i].rows != r->block : !r->b[i].rows || r->b[i].rows > r->block)
              die("corrupt index: a block does not hold the block size of rows");
      }
      if (tot != r->rows || (!n && r->rows)) die("corrupt index: the blocks do not add up to the row count"); }
    free(sz); free(sbuf);   /* both are consumed by the walk above */

    r->htab = NULL; r->hmask = 0;
    r->d = ZSTD_createDCtx();
    { uint64_t mx = 0;
      for (uint64_t i = 0; i < r->nblk; i++) {
          uint64_t t = r->b[i].lraw + r->b[i].raw;
          if (t > mx) mx = t; }
      r->blkcap = mx ? mx + 64 : (1 << 16); }
    r->blkbuf = xmalloc(r->blkcap + VGUARD);
    r->cached = UINT64_MAX;
    r->curblk = UINT64_MAX; r->curoff = 0; r->curlp = r->curip = NULL;
    return r;
}

static void ensurenames(R *r);

static void buildhash(R *r) {
    if (r->htab) return;
    ensurenames(r);
    uint64_t cap = 1; while (cap < r->nreads * 2) cap <<= 1;
    r->hmask = cap - 1;
    r->htab = xmalloc(mul(cap, 4));
    memset(r->htab, 0xff, cap * 4);
    for (uint64_t i = 0; i < r->nreads; i++) {
        uint64_t s = fnvstr(r->name[i]) & r->hmask;
        while (r->htab[s] != 0xffffffffu) {
            if (!strcmp(r->name[r->htab[s]], r->name[i])) die("a sequence name repeats: the names are not usable");
            s = (s + 1) & r->hmask;
        }
        r->htab[s] = (uint32_t)i;
    }
}

static uint64_t NAMESRC_H;
static uint64_t stathash(const struct stat *st) {
    uint64_t h = fnvbuf(FNV0, NAMEFILE, strlen(NAMEFILE));
    h = fnvbuf(h, &st->st_size, sizeof st->st_size); h = fnvbuf(h, &st->st_ino, sizeof st->st_ino);
    h = fnvbuf(h, &st->st_dev, sizeof st->st_dev); h = fnvbuf(h, &st->st_mtim, sizeof st->st_mtim);
    return fnvbuf(h, &st->st_ctim, sizeof st->st_ctim);
}
static char *readnamefile(uint64_t nreads) {
    char *names;
        FILE *nf = zfopen(NAMEFILE);
        if (!nf) die("cannot open the file given to --names");
        /* stat the path: a cookie stream has no descriptor of its own */
        struct stat st0; if (stat(NAMEFILE, &st0)) die("cannot stat the file given to --names");
        size_t cap = 1 << 20, w = 0;
        names = xmalloc(cap);
        char *ln = NULL; size_t lcap = 0; uint64_t got = 0; int fq = 2, lineno = 0;

        char *hd[3] = { NULL, NULL, NULL }; size_t hc[3] = { 0, 0, 0 }; int nhd = 0;
        while (nhd < 3 && getline(&hd[nhd], &hc[nhd], nf) > 0) nhd++;
        if (nhd >= 2) {
            if (hd[0][0] == '@' && nhd >= 3 && hd[2][0] == '+') fq = 1;
            else if (hd[0][0] == '>' && hd[1][0] != '>') fq = 0;
        }
        for (int h = 0; ; h++) {
            if (h < nhd) { free(ln); ln = hd[h]; lcap = hc[h]; hd[h] = NULL; }
            else if (getline(&ln, &lcap, nf) <= 0) break;
            lineno++;
            char *q = ln;
            if (fq == 1) { if ((lineno - 1) % 4) continue; q++; }
            else if (fq == 0) { if (*q != '>') continue; q++; }
            size_t L = strcspn(q, fq == 2 ? "\r\n" : " \t\r\n");   /* a list keeps spaces */
            if (!L) continue;
            if (L >= 512) die("--names: a name is longer than any the track could hold");
            while (w + L + 2 > cap) { cap *= 2; names = xrealloc(names, cap); }
            memcpy(names + w, q, L); w += L; names[w++] = '\n'; got++;
        }
        if (ferror(nf)) die("read error on the file given to --names");
        for (int h = 0; h < 3; h++) free(hd[h]);   /* unused when the file is shorter */
        {   struct stat st1; if (stat(NAMEFILE, &st1)) die("cannot stat the file given to --names");
            NAMESRC_H = stathash(&st0);
            if (S_ISREG(st0.st_mode) && stathash(&st1) != NAMESRC_H)
                die("the file given to --names changed while it was read");
            if (!S_ISREG(st0.st_mode)) NAMESRC_H = 0; }   /* a pipe cannot bind an index */
        if (ferror(nf)) die("read error on a companion file");
        free(ln); fclose(nf);

        if (got != nreads) {
            fprintf(stderr, "--names: the file has %llu records, the track has %llu "
                            "sequences. Rank N is the Nth record only when every "
                            "record produced rows, so this file cannot name it.\n",
                    (unsigned long long)got, (unsigned long long)nreads);
            die("name file does not match the track");
        }
        names[w] = 0;
    return names;
}

/* expand the names blob on first use */
static void ensurenames(R *r) {
    if (r->names_ready) return;
    r->names_ready = 1;
    r->name = xmalloc(mul(r->nreads, sizeof(char *)));
    if (r->nzl_keep && !r->nz_keep) {   /* deferred by ropen */
        r->nz_keep = xmalloc(r->nzl_keep);
        if (fseeko(r->f, r->nz_off, SEEK_SET)) die("seek to names");
        if (fread(r->nz_keep, 1, r->nzl_keep, r->f) != r->nzl_keep)
            die("short names");
    }
    char *names;

    if (NAMEFILE && !r->perm_ext && !r->dropped && r->nzl_keep)
        die("this track stores its own names: --names FILE would replace "
            "them, so it is refused. Unset it to read the names in the file");
    if (NAMEFILE && !r->perm_ext) {                        /* rank N is the Nth record */
        names = readnamefile(r->nreads);
    } else if (r->dropped) {

        if (!RAW && !ORDNAMES)
            die("the names of this track were dropped at encode time: give "
                "them back with --names FILE, or pass --names ordinal to "
                "accept the ordinals that are stored instead");
        names = xmalloc(mul(r->nreads, 21) + 1);
        size_t w = 0;
        for (uint64_t i = 0; i < r->nreads; i++)
            w += (size_t)sprintf(names + w, "%llu\n", (unsigned long long)i);
        names[w] = 0;
    } else if (!r->nzl_keep) {             /* ordinal: the name is the row */
        names = xmalloc(mul(r->nreads, 21) + 1);
        size_t w = 0;
        for (uint64_t i = 0; i < r->nreads; i++)
            w += (size_t)sprintf(names + w, "%llu\n", (unsigned long long)i);
        names[w] = 0;
    } else
    if (r->names_fld) {
        size_t nb; names = fldnamesd(r->nz_keep, r->nzl_keep, r->nreads, &nb);
    } else if (r->names_int) {
        size_t nb; names = intnamesd(r->nz_keep, r->nzl_keep, r->nreads, &nb);
    } else if (r->names_tok) {
        size_t nb; names = tokd(r->nz_keep, r->nzl_keep, &nb);
    } else {
        size_t nb; names = unz(r->nz_keep, r->nzl_keep, &nb, "names");
    }
    r->namebuf = names;               /* r->name[] points into this */
    char *p = names;
    uint64_t i = 0;
    for (; i < r->nreads; i++) {
        char *nl = strchr(p, '\n');
        if (!nl) break;
        r->name[i] = p; *nl = 0; p = nl + 1;
    }

    if (i != r->nreads || *p)
        die("names stream does not match the track's sequence count");
    if (r->pzl || r->perm_ext) {  /* ranks: resolve through the input order */
        char *nat;
        if (r->perm_ext) {
            if (!NAMEFILE) {
                if (RAW || ORDNAMES) goto ranks_are_names;
                die("--sort --names drop keeps the rank of every record and no "
                    "name text: give the original FASTA/FASTQ in "
                    "--names FILE, or --names ordinal to see the ranks");
            }
            nat = readnamefile(r->nreads);
        }
        else if (r->perm_tok) { size_t nb; nat = tokd(r->pz, r->pzl, &nb); }
        else {
            size_t nb; nat = unz(r->pz, r->pzl, &nb, "input-order names");
        }
        r->natbufp = nat;             /* after the remap r->name[] points here */
        char **tab = xmalloc(mul(r->nreads, sizeof(char *)));
        char *q2 = nat; uint64_t k = 0;
        for (; k < r->nreads; k++) {
            char *nl2 = strchr(q2, 10);
            if (!nl2) break;
            tab[k] = q2; *nl2 = 0; q2 = nl2 + 1;
        }
        if (k != r->nreads || *q2)
            die("input-order names do not match the track's sequence count");
        for (uint64_t j = 0; j < r->nreads; j++) {
            char *e2; unsigned long long v = strtoull(r->name[j], &e2, 10);
            if (*e2 || v >= r->nreads) die("a stored rank is not a record number");
            r->name[j] = tab[v];
        }
        free(tab);
    }
ranks_are_names: ;
}


typedef struct { uint64_t h; uint32_t i; } idxrec;

static uint64_t h64b(const char *s) {
    uint64_t h = 0x9e3779b97f4a7c15ULL;
    for (; *s; s++) { h ^= (unsigned char)*s; h *= 0xff51afd7ed558ccdULL; h ^= h >> 32; }
    return h;
}
static uint64_t namesrc(void) {
    if (!NAMEFILE) return 0;
    if (NAMESRC_H) return NAMESRC_H;
    struct stat st; if (stat(NAMEFILE, &st)) die("cannot stat the file given to --names");
    return stathash(&st);
}
static int kcmp(const void *a, const void *b) {
    uint64_t x = ((const uint64_t *)a)[0], y = ((const uint64_t *)b)[0];
    return x < y ? -1 : x > y ? 1 : 0;
}

static void rclose(R *r);

/* The name index: reads in 2^k buckets by the top k bits of a 64-bit hash of
   the name, a byte a bucket for its count and a running total every 1024
   buckets, then per read its low 32 hash bits and its ordinal in w bits, in
   bucket order, bit-packed. The 32 bits are bound to the ordinal beside them,
   so a damaged entry fails to match and the query falls back to the names.
   32 + w + 2 bits a read where the old table spent 160: a lookup still reads
   one bucket and decodes no name, and an absent name passes for a present
   one 4 times in 2^32. */
#define IXFP(h, o) ((uint32_t)((h) ^ (o) * 0x9e3779b97f4a7c15ULL))
static int ixw(uint64_t n) { int w = 1; while (w < 32 && (n - 1) >> w) w++; return w; }
static uint64_t ixget(const uint8_t *p, uint64_t bit, int nb) {
    uint64_t v; memcpy(&v, p + (bit >> 3), 8);
    return (v >> (bit & 7)) & ((1ULL << nb) - 1);
}

static void buildindex(const char *path) {
    R *r = ropen(path);
    ensurenames(r);
    if (r->hasidx) { fprintf(stderr, "karyozip: %s already carries an index\n", path); rclose(r); return; }
    if (r->nreads > 0xffffffffu)
        die("more than 2^32 sequences: the index stores the ordinal in 32 bits");
    uint64_t n = r->nreads, *key = xmalloc(mul(n ? n : 1, 16));
    int w = ixw(n ? n : 1), k = 0;
    while (k < 30 && (1ULL << (k + 2)) < n) k++;          /* about 4 reads a bucket */
    for (;; k++) {                                        /* a count must fit its byte */
        for (uint64_t i = 0; i < n; i++) {
            uint64_t h = h64b(r->name[i]);
            key[2 * i] = (k ? (h >> (64 - k)) << 32 : 0) | IXFP(h, i); key[2 * i + 1] = i;
        }
        qsort(key, n, 16, kcmp);
        uint64_t run = 0, top = 0;
        for (uint64_t i = 0; i < n; i++) { run = i && key[2 * i] >> 32 == key[2 * i - 2] >> 32 ? run + 1 : 1; if (run > top) top = run; }
        if (top <= 255) break;
    }
    uint64_t nb = 1ULL << k, ns = (nb + 1023) / 1024, eb = (n * (32 + (uint64_t)w) + 7) / 8 + 8;
    uint64_t pay = nb + ns * 8 + eb, dup = 0;
    uint8_t *t = xcalloc(pay, 1), *cnt = t, *ent = t + nb + ns * 8;
    for (uint64_t i = 0; i < n; i++) {
        cnt[key[2 * i] >> 32]++;
        if (i && key[2 * i] == key[2 * i - 2]) dup++;
        uint64_t bit = i * (32 + (uint64_t)w), v = (uint32_t)key[2 * i] | (key[2 * i + 1] << 32);
        for (int j = 0; j < 32 + w; j++) if (v >> j & 1) ent[(bit + j) >> 3] |= (uint8_t)(1 << ((bit + j) & 7));
    }
    for (uint64_t g = 0, sum = 0; g < ns; g++) {
        memcpy(t + nb + g * 8, &sum, 8);
        for (uint64_t b = g * 1024; b < nb && b < (g + 1) * 1024; b++) sum += cnt[b];
    }
    free(key);
    /* Appended to the track, magic last: a torn write leaves no magic and
       the reader simply does not see an index. One file, not two. */
    uint64_t off0 = (uint64_t)r->fsz;
    FILE *f = fopen(path, "r+b");
    if (!f || fseeko(f, 0, SEEK_END)) die("cannot append the index to the track");
    uint64_t tail[5] = { n, (uint64_t)k, r->mhash ^ namesrc(), pay, IDXMAGIC };
    int bad = fwrite(t, 1, pay, f) != pay || fwrite(tail, 8, 5, f) != 5 || fflush(f);
    if (bad) { (void)!ftruncate(fileno(f), (off_t)off0); fclose(f); die("index write failed"); }
    if (fclose(f)) die("index write failed");
    free(t);
    fprintf(stderr, "index: %llu sequences, %llu B appended, %llu fingerprint collisions\n",
            (unsigned long long)n, (unsigned long long)(pay + 40), (unsigned long long)dup);
    rclose(r);
}

/* Returns the ordinal, or -1 when absent. */
static int64_t idxlookup(R *r, const char *path, const char *nm) {
    int fd = open(path, O_RDONLY);
    if (fd < 0) return -1;
    struct stat st;
    if (fstat(fd, &st) || (size_t)st.st_size < 40) { close(fd); return -1; }
    void *m = mmap(NULL, (size_t)st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    if (m == MAP_FAILED) return -1;
    /* every word is read with memcpy, never through a cast: nothing is aligned */
    const char *tail = (const char *)m + st.st_size - 40;
    uint64_t n = ld64(tail), k = ld64(tail + 8), pay = ld64(tail + 24);
    uint64_t nb = 1ULL << (k < 40 ? k : 0), ns = (nb + 1023) / 1024;
    int w = ixw(n ? n : 1);
    if (ld64(tail + 32) != IDXMAGIC || n != r->nreads || k >= 40 || pay > (uint64_t)st.st_size - 40
        || pay != nb + ns * 8 + (n * (32 + (uint64_t)w) + 7) / 8 + 8
        || ld64(tail + 16) != (r->mhash ^ namesrc())) { munmap(m, (size_t)st.st_size); return -1; }
    const uint8_t *cnt = (const uint8_t *)tail - pay, *ent = cnt + nb + ns * 8;
    uint64_t h = h64b(nm), b = k ? h >> (64 - k) : 0, at = ld64((const char *)cnt + nb + (b >> 10) * 8);
    for (uint64_t j = b & ~1023ULL; j < b; j++) at += cnt[j];
    uint64_t hit[256], nh = 0;
    for (uint64_t e = at; e < at + cnt[b] && e < n; e++) {
        uint64_t bit = e * (32 + (uint64_t)w), o = ixget(ent, bit + 32, w);
        if (o < r->nreads && ixget(ent, bit, 32) == IXFP(h, o)) hit[nh++] = o;
    }
    munmap(m, (size_t)st.st_size);
    if (nh == 1 && !r->names_ready) return (int64_t)hit[0];   /* one fingerprint: no name to decode */
    if (nh) ensurenames(r);
    for (uint64_t j = 0; j < nh; j++) if (!strcmp(r->name[hit[j]], nm)) return (int64_t)hit[j];
    return -1;
}

static int64_t rfind(R *r, const char *nm) {
    buildhash(r);
    uint64_t s = fnvstr(nm) & r->hmask;
    while (r->htab[s] != 0xffffffffu) {
        if (!strcmp(r->name[r->htab[s]], nm)) return r->htab[s];
        s = (s + 1) & r->hmask;
    }
    return -1;
}

/* The run count heads every read on disk. Lift them into rtab and close the
   gap, so every walk below sees the layout it always saw. */
static size_t rmheads(R *r, uint8_t *b, size_t n) {
    uint8_t *w = b; const uint8_t *p = b, *e = b + n;
    uint64_t i = 0, seen = 0;
    while (p < e) {
        if (i >= r->nreads) die("corrupt block: more reads than the header claims");
        uint64_t k = vget(&p);
        if (p > e) die("corrupt block: a run count runs past the lengths");
        if (!k || k > r->rows - seen) die("corrupt block: a run count is out of range");
        r->rtab[i++] = (uint32_t)k; seen += k;
        for (uint64_t j = 0; j < k; j++) {
            const uint8_t *s0 = p;
            if (p >= e) die("corrupt block: a read runs past its lengths");
            vget(&p);
            if (p > e) die("corrupt block: a length runs past the lengths");
            size_t len = (size_t)(p - s0);
            memmove(w, s0, len); w += len;
        }
    }
    if (i != r->nreads || seen != r->rows)
        die("corrupt block: the run counts do not add up to the row count");
    r->rfilled = 1;
    return (size_t)(w - b);
}

/* One stream a read, its run count then a label beside its own length. Split
   it back into the two regions the walks expect; the counts stay in front of
   the lengths for rmheads to lift out. */
static size_t ilv_split(R *r, uint8_t *b, size_t n, size_t nl, size_t ni) {
    uint8_t *tmp = xmalloc(ni ? ni : 1);
    const uint8_t *p = b, *e = b + n;
    uint8_t *wl = b, *wi = tmp;
    uint64_t seen = 0;
    #define TAKE(dst, cap, base) do {                                        \
        const uint8_t *s0 = p;                                               \
        if (p >= e) die("corrupt block: the interleaved frame is short");     \
        vget(&p);                                                            \
        if (p > e) die("corrupt block: a value runs past the frame");         \
        size_t g = (size_t)(p - s0);                                         \
        if ((size_t)((dst) - (base)) + g > (cap))                            \
            die("corrupt block: the interleaved frame does not fit its regions"); \
        memmove((dst), s0, g); (dst) += g;                                   \
    } while (0)
    while (p < e) {
        const uint8_t *h0 = p;
        uint64_t k = vget(&p);
        if (p > e) die("corrupt block: a run count runs past the frame");
        if (!k || k > r->rows - seen) die("corrupt block: a run count is out of range");
        size_t hg = (size_t)(p - h0);
        if ((size_t)(wl - b) + hg > nl) die("corrupt block: the length region overflows");
        memmove(wl, h0, hg); wl += hg;
        for (uint64_t j = 0; j < k; j++) {
            if (p >= e) die("corrupt block: the interleaved frame is short");
            uint64_t lb = vget(&p);
            if (p > e) die("corrupt block: a label runs past the frame");
            uint8_t lv[10], *lw = lv; vput(&lw, ilv_flip(lb));
            size_t lg = (size_t)(lw - lv);
            if ((size_t)(wi - tmp) + lg > ni) die("corrupt block: the label region overflows");
            memcpy(wi, lv, lg); wi += lg;
            TAKE(wl, nl, b);
        }
        seen += k;
    }
    #undef TAKE
    if (seen != r->rows || (size_t)(wl - b) != nl || (size_t)(wi - tmp) != ni)
        die("corrupt block: the interleaved frame does not split into its two regions");
    memcpy(b + nl, tmp, ni);
    free(tmp);
    return nl;
}

/* Undo grp_pack: the labels, still in their own order, say which group each
   length came from. Runs after rmheads, so rtab holds the run counts. */
static void grp_unpack(R *r, uint8_t *lens, size_t nl, const uint8_t *ids, size_t ni) {
    uint8_t *out = xmalloc(nl + VGUARD), *w = out;
    const uint8_t *p = lens, *e = lens + nl, *q = ids, *qe = ids + ni;
    size_t cap = 64;
    uint64_t *gv = xmalloc(mul(cap, 8)), *kv = xmalloc(mul(cap, 8)), *lv = xmalloc(mul(cap, 8));
    for (uint64_t i = 0; i < r->nreads; i++) {
        uint64_t k = r->rtab[i];
        if (k > cap) { cap = (size_t)k; gv = xrealloc(gv, mul(cap, 8));
                       kv = xrealloc(kv, mul(cap, 8)); lv = xrealloc(lv, mul(cap, 8)); }
        for (uint64_t j = 0; j < k; j++) {
            if (p >= e || q >= qe) die("corrupt block: a grouped read runs past its frames");
            gv[j] = vget(&p);
            uint64_t id = vget(&q);
            if (p > e || q > qe) die("corrupt block: a grouped value runs past its frame");
            if (id >> 32) die("corrupt block: a grouped label is out of range");
            kv[j] = (id << 32) | j;
        }
        qsort(kv, (size_t)k, 8, gcmp);
        for (uint64_t j = 0; j < k; j++) lv[kv[j] & 0xffffffffu] = gv[j];
        for (uint64_t j = 0; j < k; j++) vput(&w, lv[j]);
    }
    free(gv); free(kv); free(lv);
    if (p != e || q != qe) die("corrupt block: the grouped frames do not line up");
    if ((size_t)(w - out) != nl) die("corrupt block: the grouped lengths change size");
    memcpy(lens, out, nl); free(out);
}

static void loadblk(R *r, uint64_t b) {

    if (b >= r->nblk) die("corrupt index: a block number is past the end");
    /* PAYEND and LENDEC describe the block being read, and two readers
       interleave under --ref R.kst: re-establish them even when this
       reader's block is already decoded. */
    if (r->cached == b) { PAYEND = r->payend; LENDEC = r->lenfile && r->b[b].nos != 1; return; }
    r->curblk = UINT64_MAX;   /* a new block invalidates the cursor */
    uint8_t *cb = xmalloc(r->b[b].size);
    if (fseeko(r->f, (off_t)(r->data0 + r->b[b].off), SEEK_SET)) die("seek");
    if (fread(cb, 1, r->b[b].size, r->f) != r->b[b].size)
        die("file is truncated: a block is missing or short");
    if (fnvbuf(FNV0, cb, r->b[b].size) != r->b[b].hash)
        die("corrupt block: its bytes do not match the stored hash");
    LENDEC = r->lenfile && r->b[b].nos != 1;
    uint64_t la = r->b[b].l, lb = r->b[b].size - la;
    size_t g, h2;
    PAYEND = NULL;                    /* not valid until this block decodes */
    uint64_t want = r->ilv ? r->b[b].lraw + r->b[b].raw : r->b[b].lraw;
    if (r->lencodec) { g = paydecomp(cb, la, r->lencodec, r->blkbuf, want);
                       if (!g && want) die("frame decode failed"); }
    else
    { g = ZSTD_decompressDCtx(r->d, r->blkbuf, r->blkcap, cb, la);
      if (ZSTD_isError(g)) die(ZSTD_getErrorName(g)); }
    if (g != want) die("corrupt block: a stream decodes to a size the index does not claim");
    if (r->ilv) h2 = 0;                     /* the one frame carried both */
    else if (r->idcodec) { h2 = paydecomp(cb + la, lb, r->idcodec, r->blkbuf + g, r->b[b].raw);
                      if (!h2 && r->b[b].raw) die("frame decode failed"); }
    else
    { h2 = ZSTD_decompressDCtx(r->d, r->blkbuf + g, r->blkcap - g, cb + la, lb);
      if (ZSTD_isError(h2)) die(ZSTD_getErrorName(h2)); }
    free(cb);
    if (r->ilv) {
        if (h2) die("corrupt block: a stream decodes to a size the index does not claim");
        g = ilv_split(r, r->blkbuf, g, r->b[b].lraw, r->b[b].raw);
        h2 = r->b[b].raw;
    } else if (h2 != r->b[b].raw)
        die("corrupt block: a stream decodes to a size the index does not claim");
    if (r->runhead) {
        size_t g2 = rmheads(r, r->blkbuf, g);
        memmove(r->blkbuf + g2, r->blkbuf + g, h2);
        g = g2;
        if (r->grp) grp_unpack(r, r->blkbuf, g, r->blkbuf + g, h2);
    }
    r->lp0 = r->blkbuf; r->ip0 = r->blkbuf + g;
    r->payend = PAYEND = r->blkbuf + g + h2;
    r->cached = b;
}

typedef struct { uint64_t blk, left; } span;

static uint64_t rlen(const uint8_t **lp, const char *nm, uint64_t st, int last) {
    uint64_t v = vget(lp);
    if (last && TAILDEC) {          /* holds the sequence total, not a length */
        uint64_t t = TAILPRED_R + (uint64_t)unzz(v);
        if (t <= st) die("corrupt track: sequence total is not past its last run");
        return t - st;
    }
    if (!v && !LENDEC && !RAW) die("corrupt track: a zero length outside a sentinel block");
    if (v || !LENDEC) return v;
    uint64_t e = refend_at(nm, st);
    if (e <= st) die("a length sentinel has no reference answer: wrong PAF or --ref");
    return e - st;
}

static inline void spannext(R *r, span *sp, const uint8_t **lp, const uint8_t **ip) {
    if (sp->left) return;
    loadblk(r, ++sp->blk);
    *lp = r->lp0; *ip = r->ip0; sp->left = r->b[sp->blk].rows;
}

typedef struct { uint64_t blk, off; const uint8_t *lp, *ip; } cur_t;

static void seqptrs(R *r, uint64_t ri, cur_t *c, const uint8_t **lp,
                    const uint8_t **ip, span *sp) {
    loadblk(r, getidx(r, ri).blk);
    uint64_t skip;
    rec q = getidx(r, ri);
    if (c->blk == q.blk && c->off <= q.off) {
        *lp = c->lp; *ip = c->ip; skip = q.off - c->off;
    } else {
        *lp = r->lp0; *ip = r->ip0; skip = q.off;
    }
    for (uint64_t k = 0; k < skip; k++) { vget(lp); vget(ip); }
    sp->blk = q.blk; sp->left = r->b[q.blk].rows - q.off;
}

static void seqdone(R *r, uint64_t ri, cur_t *c, const uint8_t *lp,
                    const uint8_t *ip, span *sp) {
    getidx(r, ri);          /* advances the cursor; result unused */
    c->blk = sp->blk; c->off = r->b[sp->blk].rows - sp->left;
    c->lp = lp; c->ip = ip;
}

/* the digit count by compares, then one division a digit, in place: a
   count by division doubles the divisions and a scratch array costs a
   memcpy, and this writes every coordinate of every row */
static inline char *putu64(char *o, uint64_t v) {
    static const uint64_t P[19] = { 10ULL, 100ULL, 1000ULL, 10000ULL, 100000ULL,
        1000000ULL, 10000000ULL, 100000000ULL, 1000000000ULL, 10000000000ULL,
        100000000000ULL, 1000000000000ULL, 10000000000000ULL, 100000000000000ULL,
        1000000000000000ULL, 10000000000000000ULL, 100000000000000000ULL,
        1000000000000000000ULL, 10000000000000000000ULL };
    static const char D2[201] =
        "00010203040506070809101112131415161718192021222324"
        "25262728293031323334353637383940414243444546474849"
        "50515253545556575859606162636465666768697071727374"
        "75767778798081828384858687888990919293949596979899";
    int n = (int)((63 - __builtin_clzll(v | 1)) * 1233 >> 12) + 1;
    n += (n <= 19 && v >= P[n - 1]);          /* the estimate is off by one */
    int i = n;
    while (v >= 100) { unsigned q = (unsigned)(v % 100); v /= 100; i -= 2;
                       o[i] = D2[2 * q]; o[i + 1] = D2[2 * q + 1]; }
    if (v < 10) o[--i] = (char)('0' + v);
    else { i -= 2; o[i] = D2[2 * v]; o[i + 1] = D2[2 * v + 1]; }
    return o + n;
}

static void setcols(void);
/* --db and --aln companions, for every reader */
static void loadinformed(const char *who) {
    if (ALNPAF && !ALNREF) { fprintf(stderr, "karyozip: %s--aln needs --ref\n", who); exit(1); }
    if (ALNREF && !ALNPAF) die("--ref does nothing without --aln");
    if (ALNNAMES && !ALNPAF) die("--aln-names does nothing without --aln");
    if (SETNAME && !DBFEATS) die("--set needs --db");
    if (HIERTAB && !DBFEATS) die("--hier does nothing without --db");
    if (DBFEATS) { loadfeatures(DBFEATS); if (HIERTAB) loadvocab(HIERTAB); setcols(); }
    if (ALNPAF) { if (ALNNAMES) loadordinals(ALNNAMES); REFLAZY = 1; loadref(ALNREF); loadaln(ALNPAF); }
}
static R *ropen_informed(const char *path) {
    R *r = ropen(path);
    if (DBFEATS && r->cart && !CARTON) die("this track holds label tuples outside features.tsv: give --hier HIERARCHY");
    if (DBFEATS && (!r->feathash || r->feathash != (r->cart ? HIERHASH : FEATHASH)))
        die(r->cart ? "this features.tsv and hierarchy are not the ones the track was encoded with"
                    : "this features.tsv is not the one the track was encoded with");
    return r;
}
static size_t *SETK, NSETK;   /* --set as columns of features.tsv, default all */
static void setcols(void) {
    SETK = xmalloc((nfset ? nfset : 1) * sizeof *SETK); NSETK = 0;
    if (!SETNAME || !strcmp(SETNAME, "all")) { for (size_t k = 0; k < nfset; k++) SETK[NSETK++] = k; return; }
    char *cp = xstrdup(SETNAME);
    for (char *t = strtok(cp, ","); t; t = strtok(NULL, ",")) {
        if (!strcmp(t, "all")) die("--set all takes no other set");
        size_t k = 0; while (k < nfset && strcmp(fset[k], t)) k++;
        if (k == nfset) { fprintf(stderr, "karyozip: --set: %s is not a set in features.tsv\n", t); exit(1); }
        for (size_t j = 0; j < NSETK; j++) if (SETK[j] == k) die("--set names a set twice");
        SETK[NSETK++] = k;
    }
    free(cp);
    if (!NSETK) die("--set needs at least one set");
}
static uint32_t featid(const char *L) {   /* a --db label, checked */
    uint32_t id = (uint32_t)num_arg(L, CARTON ? CARTBASE - 1 + VPROD : CARTBASE - 1, "corrupt --db track: a label is not a featureID");
    if (id && id < CARTBASE && (id >= nfid || !fid2lab[id])) die("a featureID in the track is not in the features table");
    return id;
}
/* the label a reader prints: restored under --aln, a set tuple under --db */
static char *putlab(char *o, const R *r, uint64_t f, const char *seq, uint64_t st, uint64_t len) {
    if (r->lab && f >= r->nlab) die("corrupt block: a label id past the table");
    if (!r->lab) return putu64(o, f);
    const char *L = r->lab[f];
    if (ALNPAF && L[0] == '=' && !L[1]) {
        uint32_t rl = refrow_at(seq, st, len);
        if (rl == UINT32_MAX || rl >= r->nlab) die("a delta row has no reference answer: wrong PAF or reference track");
        L = r->lab[rl];
    }
    if (!DBFEATS) { size_t n = strlen(L); memcpy(o, L, n); return o + n; }
    uint32_t id = featid(L);
    for (size_t j = 0; j < NSETK; j++) {
        const char *Lb = flab(id, SETK[j]); size_t n = strlen(Lb);
        if (j) *o++ = '\t';
        memcpy(o, Lb, n); o += n;
    }
    return o;
}
static size_t maxlablen(const R *r) {   /* the longest label putlab writes */
    size_t m = 20;
    if (!DBFEATS) { if (r->lab) for (uint64_t j = 0; j < r->nlab; j++) { size_t L = strlen(r->lab[j]); if (L > m) m = L; } return m; }
    if (6 * NSETK > m) m = 6 * NSETK;   /* the novel tuple */
    for (uint32_t id = 0; id < nfid; id++) if (fid2lab[id]) {
        size_t L = 0; for (size_t j = 0; j < NSETK; j++) L += strlen(fid2lab[id][SETK[j]]) + 1;
        if (L > m) m = L; }
    if (CARTON) for (size_t j = 0; j < NSETK; j++) { size_t L = 0;
        for (size_t i = 0; i < NVOC[SETK[j]]; i++) if (strlen(VOC[SETK[j]][i]) > L) L = strlen(VOC[SETK[j]][i]);
        m += L + 1; }
    return m;
}

static void dump(const char *path) {
    R *r = ropen(path);
    ensurenames(r);
    size_t OCAP = 1 << 22;
    char *ob = xmalloc(OCAP), *o = ob;
    size_t *lablen = NULL, maxlab = 20;      /* 20 covers an integer label */
    if (r->lab) {
        lablen = xmalloc(mul(r->nlab, sizeof(size_t)));
        for (uint64_t j = 0; j < r->nlab; j++) {
            lablen[j] = strlen(r->lab[j]);
            if (lablen[j] > maxlab) maxlab = lablen[j];
        }
    }
    for (int j = 0; j < nlabtab; j++) {      /* a restored label may be longer */
        size_t n2 = strlen(labtab[j]);
        if (n2 > maxlab) maxlab = n2;
    }
    if (maxlab + 640 > OCAP) { OCAP = maxlab + 640; ob = xrealloc(ob, OCAP); o = ob; }
    cur_t c = { UINT64_MAX, 0, NULL, NULL };
    uint64_t emitted = 0, restored = 0;
    for (uint64_t i = 0; i < r->nreads; i++) {
        const uint8_t *lp, *ip;
        span sp;
        seqptrs(r, i, &c, &lp, &ip, &sp);
        const char *nm = r->name[i];
        size_t nlen = strlen(nm);
        const qseg *Q = ALNPAF ? qlookup(nm) : NULL;
        uint64_t st = 0, nr = getidx(r, i).nruns;
        for (uint64_t k = 0; k < nr; k++) {
            spannext(r, &sp, &lp, &ip); sp.left--;
            uint64_t e = st + rlen(&lp, nm, st, k + 1 == nr && i >= TAILSKIP_R), f = vget(&ip);
            if (r->lab && f >= r->nlab) die("corrupt block: a label id past the table");
            const char *L2 = (lablen && f < r->nlab) ? r->lab[f] : NULL;
            size_t Ln = L2 ? lablen[f] : 0;
            if (ALNPAF && L2 && L2[0] == '=' && !L2[1]) {
                uint32_t rl = refrow(Q, st, e - st);
                if (rl == UINT32_MAX) die("a delta row has no reference answer: "
                                          "wrong PAF or reference track");
                L2 = labtab[rl]; Ln = strlen(L2); restored++;
            }
            size_t need = nlen + (L2 ? Ln : 20) + 48;   /* this row: name, two numbers, label */
            if ((size_t)(o - ob) + need > OCAP) {
                if (o > ob && fwrite(ob, 1, (size_t)(o - ob), stdout) != (size_t)(o - ob))
                    die("write failed");
                if (need > OCAP) { OCAP = need; ob = xrealloc(ob, OCAP); }
                o = ob;
            }
            memcpy(o, nm, nlen); o += nlen; *o++ = '\t';
            o = putu64(o, st); *o++ = '\t';
            o = putu64(o, e); *o++ = '\t';
            if (L2) { memcpy(o, L2, Ln); o += Ln; } else o = putu64(o, f);
            *o++ = '\n';
            emitted++;
            st = e;
        }
        seqdone(r, i, &c, lp, ip, &sp);
    }
    if (o > ob && fwrite(ob, 1, (size_t)(o - ob), stdout) != (size_t)(o - ob))
        die("write failed");
    free(ob); free(lablen);
    if (ALNPAF)
        fprintf(stderr, "aln: %llu rows restored from the reference\n",
                (unsigned long long)restored);

    if (emitted != r->rows) {
        fprintf(stderr, "karyozip: emitted %llu rows, the header says %llu\n",
                (unsigned long long)emitted, (unsigned long long)r->rows);
        die("file is truncated");
    }
    rclose(r);
}

static char *putrow(char *o, const char *seq, size_t nlen, uint64_t st, uint64_t en, const char *lab, size_t ln) {
    memcpy(o, seq, nlen); o += nlen; *o++ = '\t'; o = putu64(o, st); *o++ = '\t';
    o = putu64(o, en); *o++ = '\t'; memcpy(o, lab, ln); o += ln; *o++ = '\n';
    return o;
}
/* --db: the sequence's rows over the selected sets, touching equal tuples
   merged, so one set gives what its own track gives; needs the whole sequence */
static void region_set(R *r, int64_t ri, const char *seq, uint64_t lo, uint64_t hi, uint64_t halo) {
    rec q = getidx(r, ri);
    loadblk(r, q.blk);
    const uint8_t *lp = r->lp0, *ip = r->ip0;
    for (uint64_t k = 0; k < q.off; k++) { vget(&lp); vget(&ip); }
    span sp = { q.blk, r->b[q.blk].rows - q.off };
    uint64_t n = 0, cap = 1024, *st = xmalloc(cap * 8), *en = xmalloc(cap * 8);
    char **lab = xmalloc(cap * sizeof *lab), *tmp = xmalloc(maxlablen(r) + 64);   /* a --db row joins every picked set */
    uint64_t s = 0;
    for (uint64_t k = 0; k < q.nruns; k++) {
        spannext(r, &sp, &lp, &ip); sp.left--;
        uint64_t e = s + rlen(&lp, seq, s, k + 1 == q.nruns && (uint64_t)ri >= TAILSKIP_R), f = vget(&ip);
        char *end = putlab(tmp, r, f, seq, s, e - s); *end = 0;
        if (n && en[n - 1] == s && !strcmp(lab[n - 1], tmp)) { en[n - 1] = e; s = e; continue; }
        if (n == cap) { cap *= 2; st = xrealloc(st, cap * 8); en = xrealloc(en, cap * 8); lab = xrealloc(lab, cap * sizeof *lab); }
        st[n] = s; en[n] = e; lab[n] = xstrdup(tmp); n++; s = e;
    }
    uint64_t a = 0, b = n;
    if (hi) {
        uint64_t f = n, l = 0;
        for (uint64_t k = 0; k < n; k++) if (en[k] > lo && st[k] < hi) { if (k < f) f = k; l = k; }
        if (f > l) { b = 0; } else { a = f > halo ? f - halo : 0; b = n - l - 1 <= halo ? n : l + halo + 1; }
    }
    size_t nlen = strlen(seq), cap2 = 1 << 20;
    if (cap2 < nlen + maxlablen(r) + 48) cap2 = nlen + maxlablen(r) + 48;   /* the longest row fits */
    char *ob = xmalloc(cap2), *o = ob;
    for (uint64_t k = a; k < b; k++) {
        if ((size_t)(o - ob) + nlen + strlen(lab[k]) + 48 > cap2) {
            if (fwrite(ob, 1, (size_t)(o - ob), stdout) != (size_t)(o - ob)) die("write failed");
            o = ob;
        }
        o = putrow(o, seq, nlen, st[k], en[k], lab[k], strlen(lab[k]));
    }
    if (o > ob && fwrite(ob, 1, (size_t)(o - ob), stdout) != (size_t)(o - ob)) die("write failed");
    for (uint64_t k = 0; k < n; k++) free(lab[k]);
    free(st); free(en); free(lab); free(ob); free(tmp); rclose(r);
}

static void region(const char *path, const char *seq, uint64_t lo, uint64_t hi,
                   uint64_t halo) {
    R *r = ropen_informed(path);

    int64_t ri = idxlookup(r, path, seq);
    if (ri < 0) ri = rfind(r, seq);
    if (ri < 0) die("sequence not found");
    if (DBFEATS) { region_set(r, ri, seq, lo, hi, halo); return; }
    rec q = getidx(r, ri);
    loadblk(r, q.blk);
    const uint8_t *lp = r->lp0, *ip = r->ip0;

    uint64_t start = (uint64_t)q.blk * (uint64_t)r->block + q.off;
    uint64_t endrun = start + q.nruns, jump = q.blk, kskip = 0;
    uint64_t st = 0, kfirst = 0, klast = q.nruns;
    if (lo) {
        while (jump + 1 < r->nblk && (jump + 1) * (uint64_t)r->block < endrun
               && r->b[jump + 1].bp && r->b[jump + 1].bp <= lo) jump++;
    }
    span sp;
    if (jump != q.blk) {           /* the block's first run is the sequence's */
        kskip = jump * (uint64_t)r->block - start;
        st = r->b[jump].bp;
        loadblk(r, jump); lp = r->lp0; ip = r->ip0;
        sp.blk = jump; sp.left = r->b[jump].rows;
    } else {
        for (uint64_t k = 0; k < q.off; k++) { vget(&lp); vget(&ip); }
        sp.blk = q.blk; sp.left = r->b[q.blk].rows - q.off;
    }
    if (halo && (lo || hi)) {

        const uint8_t *lp2 = lp, *ip2 = ip;
        span sp2 = sp;
        uint64_t s2 = st, first = q.nruns, last = 0;
        for (uint64_t k = kskip; k < q.nruns; k++) {
            if (hi && s2 >= hi) break;
            spannext(r, &sp2, &lp2, &ip2); sp2.left--;
            uint64_t e2 = s2 + rlen(&lp2, seq, s2, k + 1 == q.nruns && (uint64_t)ri >= TAILSKIP_R); vget(&ip2);
            if (e2 > lo && (hi == 0 || s2 < hi)) {
                if (k < first) first = k;
                last = k;
            }
            s2 = e2;
        }
        if (first > last) { rclose(r); return; }   /* atexit emitflush does the flush */
        kfirst = first > halo ? first - halo : 0;
        klast = q.nruns - last - 1 <= halo ? q.nruns : last + halo + 1;   /* no wrap */
        if (klast > q.nruns) klast = q.nruns;
        lo = 0; hi = 0;
        if (kfirst < kskip) { jump = q.blk; kskip = 0; }   /* the halo reaches back */
        loadblk(r, jump); lp = r->lp0; ip = r->ip0;   /* the pre-pass moved us */
        sp.blk = jump;
        if (jump != q.blk) sp.left = r->b[jump].rows;
        else { for (uint64_t k = 0; k < q.off; k++) { vget(&lp); vget(&ip); }
               sp.left = r->b[q.blk].rows - q.off; }
        st = jump != q.blk ? r->b[jump].bp : 0;
    }

    size_t nlen = strlen(seq), maxlab = maxlablen(r);
    size_t OCAP = 1 << 20;
    if (nlen + maxlab + 48 > OCAP) OCAP = nlen + maxlab + 48;
    char *ob = xmalloc(OCAP), *o = ob;
    for (uint64_t k = kskip; k < klast; k++) {
        spannext(r, &sp, &lp, &ip); sp.left--;
        uint64_t en = st + rlen(&lp, seq, st, k + 1 == q.nruns && (uint64_t)ri >= TAILSKIP_R), f = vget(&ip);
        if (k >= kfirst && en > lo && (hi == 0 || st < hi)) {
            if ((size_t)(o - ob) + nlen + maxlab + 48 > OCAP) {
                if (fwrite(ob, 1, (size_t)(o - ob), stdout) != (size_t)(o - ob))
                    die("write failed");
                o = ob;
            }
            memcpy(o, seq, nlen); o += nlen; *o++ = '\t';
            o = putu64(o, st); *o++ = '\t';
            o = putu64(o, en); *o++ = '\t';
            o = putlab(o, r, f, seq, st, en - st);
            *o++ = '\n';
        }
        st = en;
        if (hi && st >= hi) break;   /* nothing after this can overlap */
    }
    if (o > ob && fwrite(ob, 1, (size_t)(o - ob), stdout) != (size_t)(o - ob))
        die("write failed");
    free(ob);
    rclose(r);
}

#define MAXNODE 4096
static char *NODE[MAXNODE];
static int NPAR[MAXNODE];
static int NNODE = 0;
static size_t NLEN[MAXNODE], NLEN1[MAXNODE];   /* name and --v1 name lengths */
static uint8_t *ANC = NULL;
static int16_t *LCA = NULL;
static int ROOTID = -1;
static uint64_t MAXGAP = 1000;
/* runs of context around a requested region */

static int nodeid(const char *nm, size_t n) {
    for (int i = 0; i < NNODE; i++)
        if (strlen(NODE[i]) == n && !memcmp(NODE[i], nm, n)) return i;
    if (NNODE >= MAXNODE) die("too many hierarchy nodes");
    NODE[NNODE] = strndup(nm, n);
    NPAR[NNODE] = -1;
    return NNODE++;
}

static void loadhier(const char *path, const char *fs) {
    FILE *f = zfopen(path);
    if (!f) die("cannot open the hierarchy file");
    char *line = NULL; size_t cap = 0; ssize_t len;
    while ((len = getline(&line, &cap, f)) > 0) {
        chomp(line, &len);
        if (!len) continue;
        char *t1 = memchr(line, '\t', (size_t)len); if (!t1) continue;
        char *t2 = memchr(t1 + 1, '\t', (size_t)(len - (t1 + 1 - line))); if (!t2) continue;
        if (memchr(t2 + 1, '\t', (size_t)(len - (t2 + 1 - line)))) die("hierarchy: a row has more than three columns");
        if (!strncmp(line, "feature_set\t", 12)) continue;   /* v2 header */
        if ((size_t)(t1 - line) != strlen(fs) || memcmp(line, fs, strlen(fs))) continue;
        int c = nodeid(t1 + 1, (size_t)(t2 - t1 - 1));
        int pa = nodeid(t2 + 1, (size_t)(len - (t2 + 1 - line)));
        if (NPAR[c] >= 0) die("hierarchy: a child is listed twice");
        NPAR[c] = pa;
    }
    if (ferror(f)) die("read error on a companion file");
    free(line); fclose(f);
    if (!NNODE) die("feature set not found in the hierarchy");
    for (int a = 0; a < NNODE; a++)   /* every node but categorized has a parent */
        if (NPAR[a] < 0 && strcmp(NODE[a], "categorized")) die("hierarchy: a node without a parent that is not categorized");
    for (int a = 0; a < NNODE; a++) { int n = 0;
        for (int c = a; c >= 0; c = NPAR[c]) if (++n > NNODE) die("the hierarchy has a cycle"); }
    ROOTID = nodeid("categorized", 11);
}

static void buildanc(void) {
    for (int a = 0; a < NNODE; a++) {   /* emit needs no strlen a row */
        size_t l = strlen(NODE[a]);
        NLEN[a] = l;
        NLEN1[a] = a != ROOTID && l > 9 && !strcmp(NODE[a] + l - 9, "_specific") ? l - 9 : l;
    }
    ANC = xcalloc((size_t)NNODE * NNODE, 1);
    LCA = xmalloc(sizeof(int16_t) * (size_t)NNODE * NNODE);
    for (int a = 0; a < NNODE; a++)
        for (int c = a; c >= 0; c = NPAR[c]) ANC[(size_t)a * NNODE + c] = 1;
    for (int a = 0; a < NNODE; a++)
        for (int b = 0; b < NNODE; b++) {
            int got = -1;
            for (int c = a; c >= 0; c = NPAR[c])
                if (ANC[(size_t)b * NNODE + c]) { got = c; break; }
            LCA[(size_t)a * NNODE + b] = (int16_t)got;
        }
}
#define ISANC(node, pa) (ANC[(size_t)(node) * NNODE + (pa)])

static void smooth_one(uint64_t *st, uint64_t *en, int *feat, uint64_t n) {
    if (!n) return;
    uint8_t *was = xmalloc(n);

    uint32_t *dis = xmalloc(sizeof(uint32_t) * (size_t)NNODE);
    uint32_t gen = 0;
    int64_t *rel = xmalloc(sizeof(int64_t) * n);
    memset(dis, 0, sizeof(uint32_t) * (size_t)NNODE);
    for (;;) {
        int changes = 0;
        memset(was, 0, n);
        uint64_t i = 0;
        while (i < n) {
            uint64_t ws = i;
            uint64_t nrel = 0;
            rel[nrel++] = (int64_t)i;
            int64_t firstun = -1;
            gen++;
            int lrf = feat[i];
            uint64_t lri = i, lre = en[i];

            uint64_t j = i;
            while (j + 1 < n) {
                int nf = feat[j + 1];
                if (st[j + 1] - lre > MAXGAP) break;   /* runs tile: no wrap */
                if (ISANC(lrf, nf)) {
                    if (dis[nf] == gen) break;
                    rel[nrel++] = (int64_t)(j + 1);
                    lri = j + 1; lrf = nf; lre = en[j + 1]; j++;
                } else if (!ISANC(nf, lrf)) {
                    for (int c = nf; c >= 0; c = NPAR[c]) dis[c] = gen;
                    if (firstun < 0 && !was[j + 1]) firstun = (int64_t)(j + 1);
                    j++;
                } else break;
            }
            uint64_t we = lri;

            int peak = feat[we];
            uint64_t k = we;
            while (k + 1 < n) {
                int nf = feat[k + 1];
                if (st[k + 1] - lre > MAXGAP) break;
                if (ISANC(nf, lrf)) {
                    rel[nrel++] = (int64_t)(k + 1);
                    lri = k + 1; lrf = nf; lre = en[k + 1]; k++;
                } else if (!ISANC(lrf, nf)) {
                    if (ISANC(nf, peak)) break;
                    if (firstun < 0 && !was[k + 1]) firstun = (int64_t)(k + 1);
                    k++;
                } else break;
            }
            we = lri;

            if (nrel >= 2) nrel--;            /* archive quirk, kept */
            for (uint64_t t = 0; t < nrel; t++) was[rel[t]] = 1;

            if (we > ws) {
                int l = LCA[(size_t)feat[ws] * NNODE + feat[we]];
                if (l >= 0)
                    for (uint64_t w = ws + 1; w < we; w++)
                        if (feat[w] != l && ISANC(l, feat[w])) { feat[w] = l; changes = 1; }
            }

            if (firstun >= 0) i = ((uint64_t)firstun > we && we > ws) ? we : (uint64_t)firstun;
            else if (we > i) i = we;
            else while (i < n && was[i]) i++;
        }
        if (!changes) break;
    }
    free(was); free(dis); free(rel);
}

static uint8_t *isnovel = NULL;
static int *labnodes(R *r, const char *fs) {   /* label id to hierarchy node */
    size_t k = 0;
    if (DBFEATS) { while (k < nfset && strcmp(fset[k], fs)) k++;
                   if (k == nfset) die("smooth: the set is not in features.tsv"); }
    uint64_t n = r->nlab ? r->nlab : 1;
    isnovel = xmalloc(n); int *m = xmalloc(mul(n, sizeof(int)));
    for (uint64_t i = 0; i < r->nlab; i++) {
        const char *L = r->lab[i];
        if (ALNPAF && L[0] == '=' && !L[1]) { m[i] = INT_MIN; isnovel[i] = 0; continue; }
        if (DBFEATS) {                      /* a reference label need not be a featureID */
            size_t n = strspn(L, "0123456789");
            if (!n || L[n]) { m[i] = INT_MIN + 1; isnovel[i] = 0; continue; }
            L = flab(featid(L), k); }
        isnovel[i] = (uint8_t)(strcmp(L, "novel") == 0);
        m[i] = nodeid(L, strlen(L));
    }
    return m;
}

static uint64_t loadseq(R *r, int64_t ri, uint64_t **st, uint64_t **en, int **feat,
                        uint8_t **novel, int *labmap) {
    loadblk(r, getidx(r, ri).blk);
    const uint8_t *lp, *ip;
    uint64_t skip;

    rec q = getidx(r, ri);
    if (r->curblk == q.blk && r->curoff <= q.off) {
        lp = r->curlp; ip = r->curip; skip = q.off - r->curoff;
    } else {
        lp = r->lp0; ip = r->ip0; skip = q.off;
    }
    for (uint64_t k = 0; k < skip; k++) { vget(&lp); vget(&ip); }
    span sp = { q.blk, r->b[q.blk].rows - q.off };
    uint64_t n = q.nruns;
    if (n > r->scap) {                 /* grown, never freed per sequence */
        uint64_t c = r->scap ? r->scap : 64;
        while (c < n) c *= 2;
        r->sst = xrealloc(r->sst, 8 * c); r->sen = xrealloc(r->sen, 8 * c);
        r->sfeat = xrealloc(r->sfeat, 4 * c); r->snov = xrealloc(r->snov, c);
        r->scap = c;
    }
    *st = r->sst; *en = r->sen; *feat = r->sfeat; *novel = r->snov;
    uint64_t s = 0;
    for (uint64_t k = 0; k < n; k++) {
        spannext(r, &sp, &lp, &ip); sp.left--;
        uint64_t e = s + rlen(&lp, r->name[ri], s, k + 1 == n && (uint64_t)ri >= TAILSKIP_R), f = vget(&ip);
        (*st)[k] = s; (*en)[k] = e;
        if (f >= r->nlab) die("corrupt block: a label id past the table");
        int m = labmap[f];
        if (m == INT_MIN) {                 /* stored as =: the reference row's label */
            uint32_t rl = refrow_at(r->name[ri], s, e - s);
            if (rl == UINT32_MAX || rl >= r->nlab) die("a delta row has no reference answer: wrong PAF or reference track");
            f = rl; m = labmap[f];
        }
        if (m == INT_MIN + 1) die("corrupt --db track: a label is not a featureID");
        (*novel)[k] = isnovel[f];
        (*feat)[k] = isnovel[f] ? ROOTID : m;
        s = e;
    }
    r->curblk = q.blk; r->curoff = q.off + n;
    r->curlp = lp; r->curip = ip;
    return n;
}

/* --ref R.kst: the reference as a track. Its contigs are decoded when
   first asked, so a query pays one contig, not a gigabyte of text. Its
   labels are interned in table order, the order the text would have given. */
static R *REFR; static uint32_t *REFLAB;
static int ONESEQ;                   /* region: the answer needs one contig, not the reference */
static void reffill_r(R *r, rcon *rc);
static R *refopen(const char *path) {   /* a plain track, opened as no companion and under no --names */
    const char *na = NAMESARG; int ix = INDEXING, ld = LENDEC, td = TAILDEC; uint64_t tp = TAILPRED_R, ts = TAILSKIP_R;
    NAMESARG = NULL; INDEXING = 1;
    R *r = ropen(path);
    NAMESARG = na; INDEXING = ix; LENDEC = ld; TAILDEC = td; TAILPRED_R = tp; TAILSKIP_R = ts;
    return r;
}
static void loadref_kst(const char *path) {
    uint64_t hh[16];
    {   FILE *f = fopen(path, "rb");
        if (!f || fseek(f, 8, SEEK_SET) || fread(hh, 8, 16, f) != 16) die("cannot read the reference track");
        fclose(f); }
    if (hh[1] & (ALNFLAG | DBFLAG | SIDEFLAG | DROPPED | PERMFLAG))
        die("the reference .kst must be a plain track: not --aln, not --db, and its names in the file, neither dropped nor in a sidecar");
    R *r = refopen(path);
    if (!r->rowhash) die("the reference .kst carries no row hash: encode it again with this version");
    ensurenames(r);
    REFPERBLK = REFLAZY && ONESEQ && !r->taildec;
    REFLAB = xmalloc(mul(r->nlab ? r->nlab : 1, 4));
    for (uint64_t j = 0; j < r->nlab; j++) REFLAB[j] = (uint32_t)intern(r->lab[j], strlen(r->lab[j]));
    uint64_t tot = 0;
    for (uint64_t i = 0; i < r->nreads; i++) {
        if (rcon_get(r->name[i], 0)) die("reference track: a contig's rows are not together");
        rcon *rc = rcon_get(r->name[i], 1);
        rc->ri = (int64_t)i; rc->filled = 0; rc->blkin = NULL;
        tot += getidx(r, i).nruns;
    }
    ALNHASH ^= r->rowhash;
    REFR = r;
    fprintf(stderr, "reference: %zu contigs, %llu runs\n", nrc, (unsigned long long)tot);
    if (REFLAZY) return;
    int nthr = nthreads(); if (nthr > (int)nrc) nthr = (int)nrc; if (nthr < 1) nthr = 1;
#ifdef _OPENMP
#pragma omp parallel num_threads(nthr)
#endif
    {   R *rr;                              /* every contig at once: one reader per thread */
#ifdef _OPENMP
#pragma omp critical(refopen)
#endif
        rr = refopen(path);                 /* refopen toggles plain globals: one at a time */
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (size_t i = 0; i < nrc; i++) reffill_r(rr, RC + i);
        rclose(rr);
    }
}
static void reffill(rcon *rc) {
    if (rc->filled) return;
#ifdef _OPENMP
#pragma omp critical(reffill)
#endif
    if (!rc->filled) reffill_r(REFR, rc);
}
static void reffill_r(R *r, rcon *rc) {
    {   int ld = LENDEC, td = TAILDEC; uint64_t tp = TAILPRED_R, ts = TAILSKIP_R; const uint8_t *pe = PAYEND;
        TAILDEC = r->taildec; TAILPRED_R = r->tailpred; TAILSKIP_R = r->tailskip;
        rec q = getidx(r, (uint64_t)rc->ri);
        loadblk(r, q.blk);                          /* a plain track: LENDEC is 0 in it */
        const uint8_t *lp = r->lp0, *ip = r->ip0;
        for (uint32_t k = 0; k < q.off; k++) { vget(&lp); vget(&ip); }
        span sp = { q.blk, r->b[q.blk].rows - q.off };
        size_t n = q.nruns ? q.nruns : 1;
        /* starts into 4 bytes from the first row; a start past 32 bits, which
           no chromosome reaches, moves what is there into 8 and goes on */
        rc->st32 = xmalloc(mul(n, 4)); rc->lab = xmalloc(mul(n, 4)); rc->n = rc->cap = q.nruns;
        uint64_t s = 0;
        for (uint64_t k = 0; k < q.nruns; k++) {
            spannext(r, &sp, &lp, &ip); sp.left--;
            uint64_t e = s + rlen(&lp, rc->name, s, k + 1 == q.nruns && (uint64_t)rc->ri >= TAILSKIP_R), f = vget(&ip);
            if (f >= r->nlab) die("corrupt block: a label id past the table");
            if (e <= s) die("reference track: a run of no length");
            if (rc->st32 && s >> 32) {
                rc->st = xmalloc(mul(n, 8));
                for (uint64_t i = 0; i < k; i++) rc->st[i] = rc->st32[i];
                free(rc->st32); rc->st32 = NULL;
            }
            if (rc->st32) rc->st32[k] = (uint32_t)s; else rc->st[k] = s;
            rc->lab[k] = REFLAB[f]; s = e;
        }
        rc->end = s; r->curblk = UINT64_MAX;
        LENDEC = ld; TAILDEC = td; TAILPRED_R = tp; TAILSKIP_R = ts; PAYEND = pe;
        rc->filled = 1;
    }
}

/* Per block. Without a tail-encoded row every length is a plain varint, and
   a block that opens inside a sequence records the position it opens at, so
   the block a target position falls in can be found and decoded on its own:
   1 of chr1's 19, not 19. */
static void refprep(R *r, rcon *rc) {
    rec q = getidx(r, (uint64_t)rc->ri);
    rc->row0 = (uint64_t)q.blk * r->block + q.off;
    rc->n = rc->cap = q.nruns;
    rc->st = xmalloc(mul(rc->n ? rc->n : 1, 8));
    rc->lab = xmalloc(mul(rc->n ? rc->n : 1, 4));
    rc->bend = xcalloc(r->nblk, 8);
    rc->firstblk = q.blk;
    rc->lastblk = (uint32_t)((rc->row0 + (rc->n ? rc->n - 1 : 0)) / r->block);
    rc->blkin = xcalloc(r->nblk, 1);            /* last: it is the ready flag */
}
static void refblk(R *r, rcon *rc, uint64_t b) {
    loadblk(r, b);
    const uint8_t *lp = r->lp0, *ip = r->ip0;
    uint64_t row = b * r->block, pos = r->b[b].bp, end = row + r->b[b].rows;
    if (b == rc->firstblk) {                    /* skip the rows before it */
        for (; row < rc->row0; row++) { vget(&lp); vget(&ip); }
        pos = 0;
    }
    if (end > rc->row0 + rc->n) end = rc->row0 + rc->n;
    for (; row < end; row++) {
        uint64_t len = vget(&lp), f = vget(&ip);
        if (!len) die("corrupt track: a zero length outside a sentinel block");
        if (f >= r->nlab) die("corrupt block: a label id past the table");
        rc->st[row - rc->row0] = pos;
        rc->lab[row - rc->row0] = REFLAB[f];
        pos += len;
    }
    rc->bend[b] = pos;
    rc->blkin[b] = 1;
    /* every block in: the contig is whole, and the window checks stop */
    if (++rc->nblkin == rc->lastblk - rc->firstblk + 1) {
        rc->end = rc->bend[rc->lastblk];
        rc->filled = 1;
    }
}
/* One reference read, with the caller's block state put back. */
static void refload(rcon *rc, uint64_t b, int prep) {
    R *r = REFR;
    int ld = LENDEC, td = TAILDEC; uint64_t tp = TAILPRED_R, ts = TAILSKIP_R;
    const uint8_t *pe = PAYEND;
    TAILDEC = r->taildec; TAILPRED_R = r->tailpred; TAILSKIP_R = r->tailskip;
    if (prep) refprep(r, rc); else refblk(r, rc, b);
    r->curblk = UINT64_MAX;
    LENDEC = ld; TAILDEC = td; TAILPRED_R = tp; TAILSKIP_R = ts; PAYEND = pe;
}
/* The decoded rows around tpos: [lo, hi) and the position after them. */
static void refwin(rcon *rc, uint64_t tpos, size_t *lo, size_t *hi, uint64_t *wend) {
    if (rc->filled) { *lo = 0; *hi = rc->n; *wend = rc->end; return; }
    static __thread rcon *WRC; static __thread size_t WLO, WHI;
    static __thread uint64_t WST, WEN;
    if (rc == WRC && tpos >= WST && tpos < WEN) {    /* the window last answered */
        *lo = WLO; *hi = WHI; *wend = WEN;
        return;
    }
    R *r = REFR;
    if (!rc->blkin) {
#ifdef _OPENMP
#pragma omp critical(reffill)
#endif
        if (!rc->blkin) refload(rc, 0, 1);
    }
    uint64_t b = rc->firstblk, e = rc->lastblk;      /* the block holding tpos */
    while (b < e) { uint64_t m = (b + e + 1) / 2;
                    if (r->b[m].bp <= tpos) b = m; else e = m - 1; }
    if (!rc->blkin[b]) {
#ifdef _OPENMP
#pragma omp critical(reffill)
#endif
        if (!rc->blkin[b]) refload(rc, b, 0);
    }
    *lo = b == rc->firstblk ? 0 : (size_t)(b * r->block - rc->row0);
    *hi = rc->row0 + rc->n <= (b + 1) * r->block ? rc->n : (size_t)((b + 1) * r->block - rc->row0);
    *wend = rc->bend[b];
    WRC = rc; WLO = *lo; WHI = *hi; WST = rc->st[*lo]; WEN = *wend;   /* per block: st stays 8 bytes */
}

static char out[1 << 22];
static size_t eo = 0;
static void emitflush(void) {   /* runs from atexit: die would recurse */
    if (eo && fwrite(out, 1, eo, stdout) != eo) { fprintf(stderr, "karyozip: write failed\n"); _exit(1); }
    eo = 0;
    if (fflush(stdout)) { fprintf(stderr, "karyozip: write failed\n"); _exit(1); }
}

/* One row emitter. With a FILE it flushes a fixed buffer and can write a
   row larger than that buffer straight through; without one it grows the
   caller's buffer. The run merging and the label choice are the same either
   way, which is why they used to be written out twice. */
static void emit_to(char **ob, size_t *cap, size_t *len, FILE *f,
                    const char *seq, uint64_t *st, uint64_t *en, int *feat,
                    uint8_t *novel, uint64_t n, uint64_t lo, uint64_t hi, int v1) {
    if (!n) return;
    size_t sl = strlen(seq);
    uint64_t cs = st[0], ce = en[0];
    int cf = feat[0]; uint8_t cn = novel[0];
    for (uint64_t k = 1; k <= n; k++) {
        int same = k < n && feat[k] == cf && st[k] == ce
                   && (cf != ROOTID || novel[k] == cn);
        if (same) { ce = en[k]; continue; }
        if ((hi == 0 && lo == 0) || (ce > lo && (hi == 0 || cs < hi))) {
            int root = cf == ROOTID && cn;
            const char *nm = root ? "novel" : NODE[cf];
            size_t ml = root ? 5 : (v1 ? NLEN1[cf] : NLEN[cf]);
            size_t need = sl + ml + 48;
            if (*len + need > *cap) {
                if (f) {
                    if (*len && fwrite(*ob, 1, *len, f) != *len) die("write failed");
                    *len = 0;
                } else {
                    *cap = *cap * 2 + sl + ml + 64;
                    *ob = xrealloc(*ob, *cap);
                }
            }
            if (f && need > *cap) {                /* one row is bigger than the buffer */
                char nb[48]; char *q2 = nb; q2 = putu64(q2, cs); *q2++ = '\t'; q2 = putu64(q2, ce);
                if (fwrite(seq, 1, sl, f) != sl || fputc('\t', f) < 0
                    || fwrite(nb, 1, (size_t)(q2 - nb), f) != (size_t)(q2 - nb) || fputc('\t', f) < 0
                    || fwrite(nm, 1, ml, f) != ml || fputc('\n', f) < 0) die("write failed");
            } else {
                char *q = *ob + *len;
                memcpy(q, seq, sl); q += sl; *q++ = '\t';
                q = putu64(q, cs); *q++ = '\t';
                q = putu64(q, ce); *q++ = '\t';
                memcpy(q, nm, ml); q += ml; *q++ = '\n';
                *len = (size_t)(q - *ob);
            }
        }
        if (k == n) break;
        cs = st[k]; ce = en[k]; cf = feat[k]; cn = novel[k];
    }
}

static void emit(const char *seq, uint64_t *st, uint64_t *en, int *feat,
                 uint8_t *novel, uint64_t n, uint64_t lo, uint64_t hi, int v1) {
    char *o = out; size_t cap = sizeof out;      /* f is set, so never grown */
    emit_to(&o, &cap, &eo, stdout, seq, st, en, feat, novel, n, lo, hi, v1);
}

static void emit_buf(char **ob, size_t *cap, size_t *len, const char *seq,
                     uint64_t *st, uint64_t *en, int *feat, uint8_t *novel,
                     uint64_t n, uint64_t lo, uint64_t hi, int v1) {
    emit_to(ob, cap, len, NULL, seq, st, en, feat, novel, n, lo, hi, v1);
}

static void rclose(R *r) {
    if (!r) return;
    if (r->f) fclose(r->f);
    if (r->d) ZSTD_freeDCtx(r->d);
    if (!r->name_borrowed) { free(r->namebuf); free(r->natbufp); free(r->name); }
    if (r->lab) { free(r->lab[0]); free(r->lab); }
    if (!r->blk_borrowed) { free(r->blkbuf); free(r->rtab); }
    free(r->b);
    free(r->nz_keep); free(r->pz); free(r->ibuf); free(r->ichunk);
    free(r->coff); free(r->cpendv); free(r->htab);
    free(r->sst); free(r->sen); free(r->sfeat); free(r->snov);
    free(r);
}

static void smooth_mt(const char *path, R *r0, int64_t first, int64_t last,
                      int *labmap, uint64_t lo, uint64_t hi, int v1, int nthr) {

    const int64_t CHUNK = 4096;
    int64_t nchunk = (last - first + CHUNK - 1) / CHUNK;
    _Atomic int wfail = 0;
    if (r0->nblk == 1) loadblk(r0, 0);   /* one block: decode it once, for every thread */
#ifdef _OPENMP
#pragma omp parallel num_threads(nthr)
#endif
    {
        R *rr = ropen(path);
        rr->name = r0->name; rr->names_ready = 1; rr->name_borrowed = 1;
        if (r0->nblk == 1) {
            free(rr->blkbuf); free(rr->rtab); rr->blk_borrowed = 1;
            rr->blkbuf = r0->blkbuf; rr->rtab = r0->rtab; rr->rfilled = r0->rfilled;
            rr->lp0 = r0->lp0; rr->ip0 = r0->ip0; rr->payend = r0->payend; rr->cached = 0;
        }
        size_t cap = 1 << 20, len = 0;
        char *ob = xmalloc(cap);
#ifdef _OPENMP
#pragma omp for ordered schedule(dynamic, 1)
#endif
        for (int64_t c = 0; c < nchunk; c++) {
            int64_t a = first + c * CHUNK, b = a + CHUNK;
            if (b > last) b = last;
            len = 0;
            for (int64_t ri = a; ri < b; ri++) {
                uint64_t *st, *en; int *feat; uint8_t *novel;
                uint64_t n = loadseq(rr, ri, &st, &en, &feat, &novel, labmap);
                if (hi && (!n || lo >= en[n - 1])) continue;   /* the range misses it */
                smooth_one(st, en, feat, n);   /* the whole sequence: a range only cuts the output */
                emit_buf(&ob, &cap, &len, rr->name[ri], st, en, feat, novel, n, lo, hi, v1);
            }
#ifdef _OPENMP
#pragma omp ordered
#endif
            {
                if (len && fwrite(ob, 1, len, stdout) != len) wfail = 1;
            }
        }
        free(ob);
        rclose(rr);
    }
    if (wfail) die("write failed");
}

/* smooth TRACK HIER FS [SEQ [START END]] */
static void smoothcmd(const char *path, const char *hier, const char *fs,
                      const char *seq, uint64_t lo, uint64_t hi, int v1) {
    loadhier(hier, fs);
    R *r = ropen_informed(path);
    ensurenames(r);   /* emit prints names */
    int *labmap = labnodes(r, fs);
    buildanc();   /* after labnodes: it adds the unknown labels */

    int64_t first = 0, last = (int64_t)r->nreads;
    if (seq) { first = idxlookup(r, path, seq); if (first < 0) first = rfind(r, seq);
               if (first < 0) die("sequence not found");
               last = first + 1; }

    int nthr = nthreads();
#ifndef _OPENMP
    nthr = 1;                             /* no threaded path without OpenMP */
#endif
    if (last - first < 64) nthr = 1;      /* not worth a team */
    if (nthr > 1) {
        smooth_mt(path, r, first, last, labmap, lo, hi, v1, nthr);
        rclose(r);            /* the threaded path returned without this */
        free(labmap);
        return;
    }
    for (int64_t ri = first; ri < last; ri++) {
        uint64_t *st, *en; int *feat; uint8_t *novel;
        double _t = nowsec();
        uint64_t n = loadseq(r, ri, &st, &en, &feat, &novel, labmap);
        PH_ADD(6, _t);

        if (hi && (!n || lo >= en[n - 1])) continue;
        _t = nowsec();
        smooth_one(st, en, feat, n);
        PH_ADD(4, _t);
        _t = nowsec();
        emit(r->name[ri], st, en, feat, novel, n, lo, hi, v1);
        PH_ADD(5, _t);
    }
    free(labmap); rclose(r);
}

typedef struct { FILE *f; char *ln; size_t cap;
                 char seq[512]; uint64_t st, en; char *lab; size_t labcap; int eof;
                 char *labc; size_t labccap; uint32_t labid; int set;
                 char pseq[512]; uint64_t pen; char *plab; size_t plabcap;   /* the run before */
                 char *buf; size_t bp, bn, bcap; } brd;
static void sset(char **d, size_t *cap, const char *s, size_t n) {   /* a string on the heap */
    if (n + 1 > *cap) { *cap = n + 64; *d = xrealloc(*d, *cap); }
    memcpy(*d, s, n); (*d)[n] = 0;
}

#define BRDBUF (8u << 20)
static ssize_t brd_line(brd *b, char **lnp) {
    for (;;) {
        if (b->bp < b->bn) {
            char *nl = memchr(b->buf + b->bp, '\n', b->bn - b->bp);
            if (nl) {
                *lnp = b->buf + b->bp;
                size_t len = (size_t)(nl - (b->buf + b->bp));
                b->bp += len + 1;
                return (ssize_t)len + 1;         /* count the newline, as getline does */
            }
        }
        size_t rest = b->bn - b->bp;             /* carry the partial line down */
        if (rest && b->bp) memmove(b->buf, b->buf + b->bp, rest);
        b->bn = rest; b->bp = 0;
        if (b->bn + BRDBUF > b->bcap) {
            b->bcap = b->bn + BRDBUF + 1;
            b->buf = xrealloc(b->buf, b->bcap);
        }
        size_t got = fread(b->buf + b->bn, 1, BRDBUF, b->f);
        if (!got) {
            if (ferror(b->f)) die("read error on an input BED");
            if (!b->bn) return -1;               /* real end of file */
            die("input does not end with a newline");
        }
        b->bn += got;
    }
}

static void bnext(brd *b) {
    char *line = NULL;
    ssize_t n = brd_line(b, &line);
    if (n <= 0) { b->eof = 1; return; }
    b->ln = line;
    if (b->ln[n - 1] == '\n') b->ln[--n] = 0;
    if (memchr(b->ln, 0, (size_t)n)) die("a NUL byte inside a row");
    char *t1 = memchr(b->ln, '\t', (size_t)n);            if (!t1) die("bad line");
    char *t2 = memchr(t1 + 1, '\t', (size_t)(n - (t1 + 1 - b->ln))); if (!t2) die("bad line");
    char *t3 = memchr(t2 + 1, '\t', (size_t)(n - (t2 + 1 - b->ln))); if (!t3) die("bad line");
    size_t sl = (size_t)(t1 - b->ln), ll = (size_t)(n - (t3 + 1 - b->ln));
    if (sl >= sizeof b->seq) die("name too long");
    if (b->lab) { memcpy(b->pseq, b->seq, strlen(b->seq) + 1); b->pen = b->en; sset(&b->plab, &b->plabcap, b->lab, strlen(b->lab)); }
    memcpy(b->seq, b->ln, sl); b->seq[sl] = 0;
    sset(&b->lab, &b->labcap, t3 + 1, ll);
    if (!sl) die("empty sequence name");
    b->st = bedcoord(t1 + 1, t2, "start is not a number");
    b->en = bedcoord(t2 + 1, t3, "end is not a number");
    if (b->en < b->st) die("end is before start");
    if (b->en == b->st) die("end equals start: a run must cover at least one base");
    if (b->plab && b->st == b->pen && !strcmp(b->seq, b->pseq) && !strcmp(b->lab, b->plab))
        die("two touching runs carry the same label: merge them first");
    if (DBFAST) {                       /* the same label spans several rows */
        if (!b->labc || strcmp(b->labc, b->lab)) {
            sset(&b->labc, &b->labccap, b->lab, ll);
            b->labid = si_get(&SI[b->set], b->lab, 0);
        }
    }
}

static void encode_db(int nfiles, char **files) {
    if (HIERTAB && !DBFEATS) die("--hier does nothing without --db");
    loadfeatures(DBFEATS);
    if (HIERTAB) loadvocab(HIERTAB);
    /* --db takes either form: */
    if ((size_t)nfiles != nfset) {
        if (nfiles != 1)
            die("--db: give one BED per feature set in features.tsv column "
                "order, or one merged featureID track");
        DBIDS = 1;
    }
    brd *b = xcalloc((size_t)nfiles, sizeof *b);
    for (int i = 0; i < nfiles; i++) {
        b[i].set = i;
        b[i].f = zfopen(strcmp(files[i], "-") ? files[i] : NULL);
        if (!b[i].f) { fprintf(stderr, "cannot open %s: %s\n", files[i], strerror(errno));
                       die("cannot open an input BED"); }
        bnext(&b[i]);
        if (b[i].eof) die("an input BED is empty");
    }
    ktw *w = ktw_open(stdout);
    buf keyb = {0}; char idbuf[24];
    uint64_t rows = 0, missing = 0;
    size_t seqlen = 0;
    uint32_t prevfid = UINT32_MAX; uint64_t lastend = 0;
    char lastseq[512] = ""; int haveseq = 0;
    while (!b[0].eof) {
        for (int i = 1; i < nfiles; i++)
            if (b[i].eof || strcmp(b[i].seq, b[0].seq))
                die("the input BEDs disagree on sequence order");
        uint64_t lo = b[0].st, hi = b[0].en;
        for (int i = 1; i < nfiles; i++) {
            if (b[i].st != lo) die("the input BEDs disagree on a run start");
            if (b[i].en < hi) hi = b[i].en;          /* cut at the first change */
        }

        if (!haveseq || strcmp(lastseq, b[0].seq)) {
            haveseq = 1; snprintf(lastseq, sizeof lastseq, "%s", b[0].seq);
            seqlen = strlen(lastseq);
            prevfid = UINT32_MAX;               /* a new sequence starts fresh */
        }
        uint32_t id = UINT32_MAX;
        if (DBIDS) {                       /* the label column IS the id */
            if (b[0].lab[0] == '0' && b[0].lab[1]) die("--db: a featureID with a leading zero");
            uint64_t v = num_arg(b[0].lab, CARTON ? CARTBASE - 1 + VPROD : CARTBASE - 1, "--db: a label is not a featureID; give one BED per feature set instead");
            if (v < CARTBASE && (v >= nfid || (v && !fid2lab[v]))) die("--db: featureID not in the features table");
            if (v >= CARTBASE) CARTUSED = 1;
            id = (uint32_t)v;
        } else if (DBFAST) {
            uint32_t tv[DBMAXSET];
            for (int i = 0; i < nfiles; i++) tv[i] = b[i].labid;
            id = tup_get(tv, nfset);
        }
        if (id == UINT32_MAX) {          /* fallback keeps the old behaviour */
            keyb.n = 0;
            for (int i = 0; i < nfiles; i++) {
                if (i) bput(&keyb, "\t", 1);
                bput(&keyb, b[i].lab, strlen(b[i].lab));
            }
            bput(&keyb, "", 1);
            id = fget((char *)keyb.d);
        }
        if (id == UINT32_MAX) {

            int alln = 1;
            for (int i = 0; i < nfiles && alln; i++) alln = !strcmp(b[i].lab, "novel");
            if (alln) { id = 0; missing++; }
            else if (!CARTON) die("a label tuple is not in the features table");
            else {                       /* its rank in the vocabulary */
                uint64_t rank = 0;
                for (int i = 0; i < nfiles; i++) { size_t x = vocidx((size_t)i, b[i].lab);
                    if (x == SIZE_MAX) die("a label is neither in features.tsv nor in the hierarchy of its set");
                    rank += x * VSTRIDE[i]; }
                id = CARTBASE + (uint32_t)rank; CARTUSED = 1;
            }
        }

        if (prevfid != UINT32_MAX && id == prevfid && lastend == lo)
            die("two touching runs carry the same featureID: merge them first");
        prevfid = id; lastend = hi;
        char *e2 = putu64(idbuf, id);              /* snprintf per row is slow */
        const char *lab = idbuf; size_t llen = (size_t)(e2 - idbuf);

        if (ALNPAF && llen == 1 && *lab == '=') die("a label \"=\" is reserved under --aln");
        if (SORT) sortpush(b[0].seq, seqlen, lo, hi - lo, lab, llen);
        else ktw_row(w, b[0].seq, seqlen, hi - lo, lab, llen, lo, 0);
        rows++;
        for (int i = 0; i < nfiles; i++) {
            b[i].st = hi;                            /* consume up to the cut */
            if (b[i].st == b[i].en) bnext(&b[i]);
        }
    }
    for (int i = 1; i < nfiles; i++)
        if (!b[i].eof) die("the input BEDs have different lengths");
    fprintf(stderr, "db: %llu rows from %d tracks, %llu all-novel\n",
            (unsigned long long)rows, nfiles, (unsigned long long)missing);
    if (ALNPAF)
        fprintf(stderr, "aln: %llu of %llu rows matched the reference (%.1f%%)%s",
                (unsigned long long)ALNHIT, (unsigned long long)rows,
                rows ? 100.0 * (double)ALNHIT / (double)rows : 0.0, "\n");
    sortfinish(w);
    ktw_tail(w);
    ktw_pushread(w);
    ktw_flush(w);
    ktw_close(w); free(keyb.d);
    for (int i = 0; i < nfiles; i++) {
        if (b[i].f && b[i].f != stdin) fclose(b[i].f);
        free(b[i].buf); free(b[i].lab); free(b[i].labc); free(b[i].plab);   /* ln points into buf */
    }
    free(b);
}

typedef struct { FILE *f; char *b; size_t n, cap; } ow;

static void owput(ow *w, const char *nm, size_t nl, uint64_t a, uint64_t b,
                  const char *lab, size_t ll) {
    if (w->n + nl + ll + 48 > w->cap) {
        if (fwrite(w->b, 1, w->n, w->f) != w->n) die("write failed");
        w->n = 0;
        if (nl + ll + 48 > w->cap) { w->cap = nl + ll + 64; w->b = xrealloc(w->b, w->cap); }   /* one row that is bigger */
    }
    char *o = w->b + w->n;
    memcpy(o, nm, nl); o += nl; *o++ = '\t';
    o = putu64(o, a); *o++ = '\t';
    o = putu64(o, b); *o++ = '\t';
    memcpy(o, lab, ll); o += ll; *o++ = '\n';
    w->n = (size_t)(o - w->b);
}

static char **DBTMP; static size_t ndbtmp;
static void dbtmpclean(void) { for (size_t i = 0; i < ndbtmp; i++) if (DBTMP[i]) remove(DBTMP[i]); }
static void dbtmpsig(int s) { dbtmpclean(); (void)s; _exit(1); }   /* atexit does not run on a signal */

static void decode_db(const char *path, const char *dir) {
    /* --out-dir DIR: one BED per selected set, staged then swapped in together.
       No --out-dir: the selected columns, tab-joined, merged where they all
       agree, on stdout; one set gives that set's own track. */
    int tostd = !dir;
    size_t nout = tostd ? 1 : NSETK;
    ow *o = xcalloc(nout, sizeof *o);
    uint64_t *c0 = xcalloc(nout, sizeof *c0), *c1 = xcalloc(nout, sizeof *c1);
    const char **cl = xcalloc(nout, sizeof *cl);
    size_t *cll = xcalloc(nout, sizeof *cll);
    const char *cseq = NULL; size_t cseqn = 0;
    R *r = ropen_informed(path);
    ensurenames(r);
    char pb[4096];
    if (tostd) { o[0].f = stdout; o[0].cap = 1 << 20; o[0].b = xmalloc(o[0].cap); }
    else {
        char *d = xstrdup(dir);                     /* mkdir -p, as the Rust does */
        for (char *q = d + 1; ; q++) if (*q == '/' || !*q) {
            char c = *q; *q = 0;
            if (mkdir(d, 0777) && errno != EEXIST) die("cannot create the output directory");
            *q = c; if (!c) break;
        }
        free(d);
        for (size_t i = 0; i < NSETK; i++) {
            struct stat ts; int pn = snprintf(pb, sizeof pb, "%s/%s.bed", dir, fset[SETK[i]]);
            if (pn < 0 || (size_t)pn >= sizeof pb - 16) die("--out-dir: the output path is too long");
            if (!lstat(pb, &ts) && !S_ISREG(ts.st_mode)) die("--out-dir: a set's BED exists and is not a regular file");
        }
        DBTMP = xcalloc(NSETK, sizeof *DBTMP); ndbtmp = NSETK; atexit(dbtmpclean);
        signal(SIGINT, dbtmpsig); signal(SIGTERM, dbtmpsig); signal(SIGHUP, dbtmpsig);
        for (size_t i = 0; i < NSETK; i++) {
            int pn = snprintf(pb, sizeof pb, "%s/%s.bed.tmp.%d", dir, fset[SETK[i]], (int)getpid());
            if (pn < 0 || (size_t)pn >= sizeof pb - 16) die("--out-dir: the output path is too long");
            DBTMP[i] = xstrdup(pb);
            o[i].f = fopen(pb, "w");
            if (!o[i].f) die("cannot write a derived BED");
            o[i].cap = 1 << 20;
            o[i].b = xmalloc(o[i].cap);
        }
    }
    char *jb[2] = { NULL, NULL }; size_t jcap[2] = { 0, 0 }; int jw = 0;   /* joined labels, pending and current */
    cur_t c = { UINT64_MAX, 0, NULL, NULL };
    uint64_t rows = 0;
    for (uint64_t i = 0; i < r->nreads; i++) {
        const uint8_t *lp, *ip;
        span sp;
        seqptrs(r, i, &c, &lp, &ip, &sp);
        const char *nm = r->name[i];
        const qseg *Q = ALNPAF ? qlookup(nm) : NULL;

        if (cseq)
            for (size_t j = 0; j < nout; j++)
                owput(&o[j], cseq, cseqn, c0[j], c1[j], cl[j], cll[j]);
        cseq = NULL;
        uint64_t s = 0, nr = getidx(r, i).nruns;
        size_t nml = strlen(nm);
        for (uint64_t k = 0; k < nr; k++) {
            spannext(r, &sp, &lp, &ip); sp.left--;
            uint64_t len = rlen(&lp, nm, s, k + 1 == nr && i >= TAILSKIP_R), f = vget(&ip);
            if (r->lab && f >= r->nlab) die("corrupt block: a label id past the table");
            const char *idtxt = (r->lab && f < r->nlab) ? r->lab[f] : NULL;

            if (ALNPAF && idtxt && idtxt[0] == '=' && !idtxt[1]) {
                uint32_t rl = refrow(Q, s, len);
                if (rl == UINT32_MAX)
                    die("a delta row has no reference answer: wrong PAF or --ref");
                idtxt = labtab[rl];
            }
            /* the label text IS the featureID: parse it as strictly as encode did */
            uint32_t id = idtxt ? featid(idtxt) : (uint32_t)f;
            for (size_t j = 0; j < nout; j++) {
                const char *Lb; size_t Lbn;
                if (tostd) {                        /* the selected columns, one tuple */
                    size_t need = 1; for (size_t t = 0; t < NSETK; t++) need += strlen(flab(id, SETK[t])) + 1;
                    if (need > jcap[jw]) { jcap[jw] = need * 2; jb[jw] = xrealloc(jb[jw], jcap[jw]); }
                    char *q = jb[jw];
                    for (size_t t = 0; t < NSETK; t++) {
                        const char *L = flab(id, SETK[t]); size_t L2 = strlen(L);
                        if (t) *q++ = '\t';
                        memcpy(q, L, L2); q += L2;
                    }
                    *q = 0; Lb = jb[jw]; Lbn = (size_t)(q - jb[jw]);
                } else { Lb = flab(id, SETK[j]); Lbn = strlen(Lb); }
                if (cseq && c1[j] == s && cll[j] == Lbn && !memcmp(cl[j], Lb, Lbn)) {
                    c1[j] = s + len;                       /* extend the run */
                } else {
                    if (cseq) owput(&o[j], cseq, cseqn, c0[j], c1[j], cl[j], cll[j]);
                    c0[j] = s; c1[j] = s + len; cl[j] = Lb; cll[j] = Lbn;
                    if (tostd) jw ^= 1;                    /* the pending label lives on */
                }
            }
            cseq = nm; cseqn = nml;
            s += len; rows++;
        }
        seqdone(r, i, &c, lp, ip, &sp);
    }
    if (cseq)
        for (size_t j = 0; j < nout; j++)
            owput(&o[j], cseq, cseqn, c0[j], c1[j], cl[j], cll[j]);
    for (size_t j = 0; j < nout; j++) {
        if ((o[j].n && fwrite(o[j].b, 1, o[j].n, o[j].f) != o[j].n) || (tostd ? fflush(o[j].f) : fclose(o[j].f))) die("write failed");
        free(o[j].b);
    }
    free(jb[0]); free(jb[1]);
    if (tostd) fprintf(stderr, "db: %llu rows, %zu sets on stdout\n", (unsigned long long)rows, NSETK);
    else {
        char **old = xcalloc(NSETK, sizeof *old); size_t nmv = 0; int bad = 0;
        for (size_t j = 0; j < NSETK && !bad; j++) {
            snprintf(pb, sizeof pb, "%s/%s.bed", dir, fset[SETK[j]]);
            char *o2; if (asprintf(&o2, "%s/%s.bed.old.%d", dir, fset[SETK[j]], (int)getpid()) < 0) die("oom");
            if (!rename(pb, o2)) { old[j] = o2; } else { free(o2); if (errno != ENOENT) bad = 1; }
        }
        for (size_t j = 0; j < NSETK && !bad; j++) {
            snprintf(pb, sizeof pb, "%s/%s.bed", dir, fset[SETK[j]]);
            if (rename(DBTMP[j], pb)) bad = 1; else { free(DBTMP[j]); DBTMP[j] = NULL; nmv++; }
        }
        if (bad) {
            for (size_t j = 0; j < nmv; j++) {   /* take the new ones back out */
                snprintf(pb, sizeof pb, "%s/%s.bed", dir, fset[SETK[j]]);
                if (asprintf(&DBTMP[j], "%s/%s.bed.tmp.%d", dir, fset[SETK[j]], (int)getpid()) >= 0)
                    (void)!rename(pb, DBTMP[j]);
            }
            for (size_t j = 0; j < NSETK; j++) if (old[j]) {
                snprintf(pb, sizeof pb, "%s/%s.bed", dir, fset[SETK[j]]);
                (void)!rename(old[j], pb); free(old[j]);
            }
            free(old); die("cannot rename a derived BED into place");
        }
        for (size_t j = 0; j < NSETK; j++) if (old[j]) { remove(old[j]); free(old[j]); }
        free(old);
        fprintf(stderr, "db: %llu rows into %zu tracks under %s\n", (unsigned long long)rows, NSETK, dir);
    }
    free(o); free(c0); free(c1); free(cl); free(cll); rclose(r);
}

#ifndef KARYOZIP_LIB
/* -r NAME or NAME:START-END, the range BED's: 0-based, END excluded. A
   name may hold ':' itself, so only a digits-digits tail is a range. */
static const char *regparse(char *r, uint64_t *lo, uint64_t *hi) {
    *lo = *hi = 0;
    char *c = strrchr(r, ':'), *d = c ? strchr(c + 1, '-') : NULL;
    if (!d || d == c + 1 || !d[1] || strspn(c + 1, "0123456789") != (size_t)(d - c - 1)
        || strspn(d + 1, "0123456789") != strlen(d + 1)) return r;
    *c = *d = 0;
    *lo = num_arg(c + 1, UINT64_MAX - 1, "-r: START is not a number");
    *hi = num_arg(d + 1, UINT64_MAX, "-r: END is not a number");
    if (*lo >= *hi) die("-r: END must be past START");
    return r;
}

/* -i and -o: the file takes the place of stdin or stdout */
static void redirect(const char *path, int fd) {
    int f = fd ? open(path, O_WRONLY | O_CREAT | O_TRUNC, 0644) : open(path, O_RDONLY);
    if (f < 0 || dup2(f, fd) < 0) { fprintf(stderr, "karyozip: cannot open %s\n", path); exit(1); }
    close(f);
}

int main(int argc, char **argv) {
    /* glibc raises its mmap threshold to the first large mapped chunk it
       frees, after which every big buffer comes from a per-thread heap and
       stays resident: 113 MB against 57 on one -e high encode. Pin it. */
    mallopt(M_MMAP_THRESHOLD, 128 << 10);
    signal(SIGXFSZ, SIG_IGN);
    for (int fd = 0; fd < 3; fd++) if (fcntl(fd, F_GETFD) < 0) die("stdin, stdout or stderr is closed");
    atexit(profdump);
    atexit(emitflush);
    static const char *const USAGE = "usage:\n"
            "  karyozip encode [-i IN.bed[.gz]] [-o OUT.kst] [-e low|medium|high]\n"
            "                  [-n drop|drop:FILE] [--sort[=MEM]]\n"
            "                  [-d features.tsv [--hier H] -i SET.bed ...] [-a P.paf --ref R.bed|R.kst [--aln-names F]]\n"
            "  karyozip decode -i F [-o OUT.bed] [-r NAME[:START-END] [--halo N]] [-n FILE|ordinal] [--raw]\n"
            "                  [-d features.tsv [--hier H] [--set A,B] [--out-dir DIR]] [-a P.paf --ref R]\n"
            "  karyozip smooth -i F --hier H -f FEATURE_SET [-r NAME[:START-END]] [-o OUT]\n"
            "  karyozip index  -i F                  appends a name index to the track itself\n"
            "  karyozip aln-pack -a P.paf --ref REF.bed|REF.kst -o OUT   only what -a reads\n"
            "  encode reads stdin without -i; every command but index writes stdout without -o\n"
            "  any command: -t N, --threads N  --profile  --tmp-dir DIR  -h, --help\n"
            "               --version prints the zstd, libbz3 and htscodecs versions\n"
            "  encode:      -e, --effort: low is a quick encode with blocks and an index for fast\n"
            "               queries, medium the default, high the smallest file\n"
            "               -b N|all sets the blocks by hand\n"
            "  -r NAME:START-END is BED, 0-based, END excluded; -r NAME alone is the whole sequence;\n"
            "               decode -r prints the runs overlapping it, --halo N more each side\n"
            "  readers:     -n, -d, --hier, -a, --ref as encode was given; -n FILE takes the names from\n"
            "               drop:FILE, a FASTA, FASTQ or list; -n ordinal accepts dropped names\n"
            "               --set A,B picks feature sets of a -d track, default all; without\n"
            "               --out-dir decode prints them as columns, merged where they all agree\n"
            "               --raw shows what is stored, reversing no transform\n"
            "  with -a:     --aln-names FILE, on encode and decode, when the track names rows by ordinal\n"
            "  smooth:      --max-gap N  --v1\n"
            "  long forms:  -i --input  -o --output  -r --region  -f --feature-set  -b --block\n"
            "               -e --effort  -n --names  -d --db  -a --aln\n"
            "  Two variables are read, and neither changes the stored bytes:\n"
            "    TMPDIR         where temp files go, when --tmp-dir is not given\n"
            "    OMP_NUM_THREADS  the thread count, when --threads is not given";
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) { puts(USAGE); return 0; }
        /* The codec versions decide the stored bytes, so print them. */
        if (!strcmp(argv[i], "--version")) {
            printf("karyozip %s, format %.8s\n", KARYOZIP_VERSION, MAGIC);
            printf("  zstd       %s\n", ZSTD_versionString());
            printf("  libbz3     %s\n", bz3_version());
            printf("  htscodecs  %s\n", HTSCODECS_VERSION_TEXT);
            return 0;
        }
    }
    if (argc < 2) die(USAGE);
    const char *INS[64]; int NIN = 0;
    const char *OUTPUT = NULL, *HALOA = NULL, *FSA = NULL, *HIERA = NULL;
    char *REGA = NULL;
    static const struct { const char *f; int nv; } FLAGS[] = {
        {"--version",0},{"--threads",1},{"-t",1},{"--profile",0},{"--tmp-dir",1},{"--names",1},
        {"--max-gap",1},{"--v1",0},{"--aln-names",1},{"--sort",0},
        {"--effort",1},
        {"--block",1},{"--db",1},{"--aln",1},
        {"--ref",1},{"--out-dir",1},{"--set",1},{"--raw",0},{"--hier",1},
        {"--input",1},{"--output",1},{"--region",1},{"--halo",1},{"--feature-set",1},
    };
    static const char *const SHORT[][2] = {   /* each short form is its long one */
        {"-i","--input"},{"-o","--output"},{"-r","--region"},{"-f","--feature-set"},
        {"-b","--block"},{"-e","--effort"},{"-n","--names"},{"-d","--db"},{"-a","--aln"},
    };
    const char *cmd = "";   /* the first word that is not a flag or its value */
    for (int i = 1; i < argc; i++) {
        for (size_t k = 0; k < sizeof SHORT / sizeof *SHORT; k++)
            if (!strcmp(argv[i], SHORT[k][0])) { argv[i] = (char *)SHORT[k][1]; break; }
        const char *a = argv[i];
        if (a[0] != '-' || !a[1]) { cmd = a; break; }
        int nv = -1;
        if (!strncmp(a, "--sort=", 7)) nv = 0;        /* the value is in the flag */
        else for (size_t k = 0; k < sizeof FLAGS / sizeof *FLAGS; k++)
            if (!strcmp(a, FLAGS[k].f)) { nv = FLAGS[k].nv; break; }
        if (nv < 0) { fprintf(stderr, "karyozip: unknown flag %s\n", a); exit(1); }
        i += nv;
    }
    for (int i = 1; i < argc; i++)   /* short forms past the command too */
        for (size_t k = 0; k < sizeof SHORT / sizeof *SHORT; k++)
            if (!strcmp(argv[i], SHORT[k][0])) { argv[i] = (char *)SHORT[k][1]; break; }
    int enc = !strcmp(cmd, "encode");
    int sm = !strcmp(cmd, "smooth");
    int dec = !strcmp(cmd, "decode");
    int rd = sm || dec || !strcmp(cmd, "index");
    int inf = enc || rd;                   /* --db and --aln: every reader but index */
    if (!strcmp(cmd, "index")) inf = 0;
    int pack = !strcmp(cmd, "aln-pack");
    for (int i = 1; i < argc; i++) {
        int take = 0, ok = 1; const char *a = argv[i], *v = i + 1 < argc ? argv[i + 1] : NULL;
        for (size_t k = 0; k < sizeof FLAGS / sizeof *FLAGS; k++)   /* a value is required */
            if (!strcmp(a, FLAGS[k].f) && FLAGS[k].nv && !v) {
                fprintf(stderr, "karyozip: %s needs a value\n", a); exit(1); }
        if (!strcmp(a, "--threads") || !strcmp(a, "-t")) {
            NTHREADS = (int)int_arg(v, 1, 1024, "--threads is a number of 1 to 1024"); take = 2; }
        else if (!strcmp(a, "--profile")) { PROF = 1; take = 1; }
        else if (!strcmp(a, "--tmp-dir")) { TMPD = v; take = 2; }
        else if (!strcmp(a, "--names")) {
            if (!strcmp(v, "drop")) { ok = enc; NAMEMODE = 1; }
            else if (!strncmp(v, "drop:", 5)) { ok = enc; NAMES_OUT = v + 5;
                                                if (!*NAMES_OUT) die("--names drop:FILE needs a file"); }
            else if (!strcmp(v, "ordinal")) {
                if (enc) die("encode takes --names drop or drop:FILE");
                ok = rd; ORDNAMES = 1; }
            else if (enc) die("encode takes --names drop or drop:FILE; a name file is --aln-names");
            else { ok = rd; NAMESARG = v; }   /* the blob drop:FILE wrote, or a FASTA, FASTQ or list */
            take = 2; }
        else if (!strcmp(a, "--max-gap")) { ok = sm; MAXGAP = num_arg(v, UINT64_MAX, "--max-gap needs a number"); take = 2; }
        else if (!strcmp(a, "--v1")) { ok = sm; V1 = 1; take = 1; }
        else if (!strcmp(a, "--aln-names")) {   /* the alignment's own table: ordinal N is line N+1 */
            ok = inf; ALNNAMES = v; take = 2; }
        else if (!strcmp(a, "--sort")) { ok = enc; SORT = 1; take = 1; }
        else if (!strncmp(a, "--sort=", 7)) { ok = enc; SORT = 1; SORTMEM = mem_arg(a + 7); take = 1; }
        else if (!strcmp(a, "--effort")) { ok = enc;
            if (strcmp(v, "low") && strcmp(v, "medium") && strcmp(v, "high")) die("--effort is low, medium or high");
            FASTMODE = !strcmp(v, "low"); SLOW = !strcmp(v, "high"); take = 2; }
        else if (!strcmp(a, "--block")) { ok = enc;
            BLOCK = !strcmp(v, "all") ? UINT64_MAX      /* one block, no random access */
                  : (uint64_t)int_arg(v, 1, 0x7fffffffffffffffL,
                                      "--block is a number of runs, 1 or more, or \"all\"");
            BLOCKSET = 1; take = 2; }
        else if (!strcmp(a, "--db")) { ok = inf; DBFEATS = v; take = 2; }
        else if (!strcmp(a, "--hier")) { ok = inf; HIERA = v; take = 2; }
        else if (!strcmp(a, "--aln")) { ok = inf || pack; ALNPAF = v; take = 2; }
        else if (!strcmp(a, "--ref")) { ok = inf || pack; ALNREF = v; take = 2; }
        else if (!strcmp(a, "--out-dir")) { ok = dec; OUTDIR = v; take = 2; }
        else if (!strcmp(a, "--input")) { ok = !pack;
            if (NIN == 64) die("too many input files");
            INS[NIN++] = v; take = 2; }
        else if (!strcmp(a, "--output")) { ok = strcmp(cmd, "index") != 0; OUTPUT = v; take = 2; }
        else if (!strcmp(a, "--region")) { ok = sm || dec; REGA = (char *)v; take = 2; }
        else if (!strcmp(a, "--halo")) { ok = dec; HALOA = v; take = 2; }
        else if (!strcmp(a, "--feature-set")) { ok = sm; FSA = v; take = 2; }
        else if (!strcmp(a, "--set")) { ok = dec; SETNAME = v; take = 2; }
        else if (!strcmp(a, "--raw")) { ok = dec; RAW = 1; take = 1; }
        if (take && !ok) { if (!*cmd) die(USAGE);   /* flags but no command */
            fprintf(stderr, "karyozip: %s does not apply to %s\n", a, cmd); exit(1); }
        if (take) { for (int j = i; j + take <= argc; j++) argv[j] = argv[j + take]; argc -= take; i--; }
    }
    if (!*cmd) die(USAGE);   /* argv[1] is now past the end */
    if (!enc && !rd && !pack) die("bad command or missing arguments");
    if (FASTMODE && !BLOCKSET) BLOCK = 262144;   /* -e low: blocks and an index, for queries */
    for (int i = 2; i < argc; i++) {   /* every argument is a flag now */
        if (argv[i][0] == '-' && argv[i][1]) fprintf(stderr, "karyozip: unknown flag %s\n", argv[i]);
        else fprintf(stderr, "karyozip: unexpected argument %s\n", argv[i]);
        exit(1); }
    if (NTHREADS > 1024) NTHREADS = 1024;   /* libgomp keeps team state on the stack */
#ifdef _OPENMP
    if (omp_get_max_threads() > 1024) omp_set_num_threads(1024);   /* OMP_NUM_THREADS too */
#endif


    WALL0 = nowsec();

    if (sm) { HIERTAB = DBFEATS ? HIERA : NULL; } else HIERTAB = HIERA;   /* smooth's own, and --db's */
    if (NIN > 1 && !(enc && DBFEATS)) die("give one -i");
    if (OUTPUT && !pack) redirect(OUTPUT, 1);
    if (pack) {
        if (!ALNPAF || !ALNREF || !OUTPUT) die("usage: karyozip aln-pack -a P.paf --ref REF.bed|REF.kst -o OUT");
        loadref(ALNREF);
        alnpack(ALNPAF, OUTPUT);
        return 0;
    }
    if (enc) {
        if (NAMEMODE == 1 && NAMES_OUT) die("--names drop and --names drop:FILE are exclusive");
        if (ALNPAF && !ALNREF) die("--aln needs --ref");
        if (ALNREF && !ALNPAF) die("--ref does nothing without --aln");

        if (ALNPAF && NAMEMODE == 1 && !ALNNAMES)
            die("--aln with --names drop needs --aln-names: without the names "
                "decode cannot reach the alignment");

        if (ALNPAF && SORT && NAMEMODE == 1)
            die("--aln with --sort and --names drop cannot be decoded: sorting "
                "reorders the records, so the ordinal name map is invalid");
        if (ALNNAMES && !ALNPAF) die("--aln-names does nothing without --aln");

        PERMON = SORT && !NAMES_OUT;
        if (ALNPAF) { if (ALNNAMES) loadordinals(ALNNAMES);
                      loadref(ALNREF); loadaln(ALNPAF); }

        if (DBFEATS) {
            if (!NIN) die("--db: give one -i BED per feature set");
            encode_db(NIN, (char **)INS);
        } else { if (NIN) redirect(INS[0], 0); encode(); }
    }
    else if (!NIN) { fprintf(stderr, "karyozip: %s needs -i TRACK\n", cmd); exit(1); }
    else if (dec && !REGA) {
        if (HALOA) die("--halo needs -r NAME:START-END");
        loadinformed("decode ");
        if (DBFEATS) decode_db(INS[0], OUTDIR); else dump(INS[0]);
    }
    /* one sequence: -r NAME is all of it, NAME:START-END the runs
       overlapping it, plus an optional halo */
    else if (dec) {
        if (OUTDIR) die("-r and --out-dir: give one");
        if (RAW) die("-r and --raw: give one");
        uint64_t rlo, rhi; const char *nm = regparse(REGA, &rlo, &rhi);
        if (HALOA && !rhi) die("--halo needs -r NAME:START-END");
        ONESEQ = 1;   /* the answer needs one contig, not the whole reference */
        loadinformed("decode ");
        region(INS[0], nm, rlo, rhi, HALOA ? num_arg(HALOA, UINT64_MAX, "--halo is not a number") : 0);
    }
    else if (!strcmp(cmd, "index")) { INDEXING = 1; buildindex(INS[0]); }
    else if (sm) {
        if (!HIERA || !FSA) die("smooth needs --hier HIER and -f FEATURE_SET");
        uint64_t slo = 0, shi = 0; const char *nm = REGA ? regparse(REGA, &slo, &shi) : NULL;
        loadinformed("smooth ");
        smoothcmd(INS[0], HIERA, FSA, nm, slo, shi, V1);
    }
    else die("bad command or missing arguments");
    return 0;
}
#endif   /* KARYOZIP_LIB: the library has no main */
