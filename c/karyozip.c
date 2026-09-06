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
#include "karyozip.h"
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>
#include <fcntl.h>
#include <signal.h>
#include <unistd.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static uint64_t BLOCK = 262144;
static int BLOCKSET = 0;                 /* --block was given */
/* 0 means use every core; --threads overrides */
static int NTHREADS = 0;
static int LEVEL = 17;
static int LEVELSET = 0;              /* --level was given: do not override */
static int LEVEL_L = -1, LEVEL_I = -1;   /* per stream, chosen by the search */

static const char MAGIC[8] = "KST00001";   /* KaryoScope Track, format version 1 */
#define IDXMAGIC 0x313058444954534bULL   /* "KSTIDX01", the index trailer */
#define KARYOZIP_VERSION "0.2.0"

static int PROF = 0;
static double PH[9];
static const char *PHN[9] = {"parse", "compress", "names", "index",
                             "smooth", "parse:lines", "parse:append", "write", "codec search"};
static inline double nowsec(void) {
    struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}
#define PH_RSS(i) do { if (PROF) PHRSS[i] = rss_now(); } while (0)
/* the high-water mark is read from /proc at most every 10 ms: a span per
   sequence over 3,000,000 sequences made a 1.7 s smooth take 104 s */
#define PH_ADD(i, t0) do { if (PROF) { double _n = nowsec(); PH[i] += _n - (t0); \
    if (_n - PHLAST > 0.01) { PHRSS[i] = rss_now(); PHLAST = _n; } } } while (0)
static size_t PHRSS[9]; static double PHLAST;
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
    for (int i = 0; i < 9; i++) tot += PH[i];
    for (int i = 0; i < 9; i++)
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
    int g = gzread((gzFile)c, b, (unsigned)n);
    return g < 0 ? -1 : g;
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
static uint32_t ld32(const void *p) { uint32_t v; memcpy(&v, p, 4); return v; }

static long int_arg(const char *p, long lo, long hi, const char *what) {   /* digits only, in range */
    if (!p || !*p || strspn(p, "0123456789") != strlen(p) || strlen(p) > 18) die(what);
    long v = strtol(p, NULL, 10);
    if (v < lo || v > hi) die(what);
    return v;
}
static long codec_arg(const char *p, const char *what) {   /* a codec id: digits, or 0x and hex digits */
    int hex = p && p[0] == '0' && (p[1] == 'x' || p[1] == 'X'); if (hex) p += 2;
    if (!p || !*p || strlen(p) > 8 || strspn(p, hex ? "0123456789abcdefABCDEF" : "0123456789") != strlen(p)) die(what);
    long v = strtol(p, NULL, hex ? 16 : 10);
    if (v > 0x7fffffff) die(what);
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
    uint64_t v = 0; int s = 0; uint8_t b;

    do { if (PAYEND && *p >= PAYEND) die("corrupt block: a run reads past the payload");
         b = *(*p)++;
         if (s < 64) v |= (uint64_t)(b & 0x7f) << s;
         s += 7;
         if (s > 70) die("corrupt stream: a varint does not terminate");
    } while (b & 0x80);
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

#define LABHASH (1 << 16)
static char **labtab = NULL;
static int nlabtab = 0, labcap = 0;
static int32_t labhash[LABHASH];
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
    ck ^= (uint32_t)(uint8_t)s[0] << 7;
    ck ^= (uint32_t)(uint8_t)s[n - 1] << 15;
    ck ^= (uint32_t)(uint8_t)s[n >> 1] << 23;
    ck = (ck ^ (ck >> 13)) & (LCACHE - 1);
    int cid = lcache[ck];
    if (cid >= 0 && labnlen[cid] == n && !memcmp(labtab[cid], s, n)) return lastid = cid;
    if (!labhash_init) { memset(labhash, 0xff, sizeof labhash); labhash_init = 1; }
    uint64_t k = fnvbuf(FNV0, s, n) & (LABHASH - 1);
    uint64_t probes = 0;
    while (labhash[k] >= 0) {
        if (++probes >= LABHASH) die("more than LABHASH distinct labels");
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

#define METALEVEL (LEVEL < 19 ? LEVEL : 19)

#define ICHUNK 65536

#define TOKFLAG (1ULL << 63)

#define CODEC_ARITH 0x10000
#define CODEC_BZIP3 0x20000           /* 0x20000|MB: bzip3, block size in MB */
#define DROPPED (1ULL << 60)          /* names deliberately not stored */
#define ALNFLAG (1ULL << 59)          /* rows deltaed against a reference */
#define DBFLAG  (1ULL << 58)          /* labels are featureIDs, not labels */
#define CARTFLAG (1ULL << 52)         /* some featureIDs are 2^26 + a vocabulary rank */
#define PERMFLAG (1ULL << 54)         /* names in input order + a permutation */
#define PERMTOK  (1ULL << 53)         /* that blob went through the tokeniser */
#define SIDEFLAG (1ULL << 57)
#define LENFLAG  (1ULL << 56)         /* lengths the reference implies are 0 */
#define TAILFLAG (1ULL << 55)

static int NAMEMODE = 0;              /* 0 auto, 1 drop, 2 keep */
static const char *NAMEFILE = NULL;
static const char *NAMESIN = NULL;
/* One flag, two kinds of file: the blob "--names drop:F" wrote, or the
 * sequences themselves. The track says which, so the file is never sniffed
 * (a pipe cannot be reopened). Resolved in ropen, see NAMESARG. */
static const char *NAMESARG = NULL;
/* the default trades size for decode speed; --slow does not */
static int SLOW = 0;                  /* --slow: every codec over every block */
static int FASTMODE = 1;              /* the default: quick-decoding codecs */
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
static const char *OUTDIR = NULL;   /* --out DIR, for decode --db */
static const char *SETNAME = NULL;  /* --set NAME: one set of a --db track */
static int INDEXING = 0;            /* the index command: names only, no companions */
static int ORDNAMES = 0;              /* --names ordinal: accept dropped names */
static const char *ALNREF = NULL;     /* --ref reference track */
static const char *ALNNAMES = NULL;   /* --names FILE: ordinal N is line N+1 */
static char **ORD; static size_t nord;

static uint64_t BYT_LENS = 0, BYT_LABS = 0;
static int IDCODEC = -1, LENCODEC = -1;   /* -1 picks the winner on a sample */
#define INTFLAG (1ULL << 62)

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
        if (p == end) { free(v.d); return 0; }
        p++;
        int64_t d = x - prev; prev = x;

        bvar(&v, ((uint64_t)d << 1) ^ (uint64_t)(d >> 63));
    }
    if (p != end) { free(v.d); return 0; }
    size_t z = zc(&v, out); free(v.d);
    return z;
}

static char *intnamesd(const uint8_t *b, size_t n, uint64_t cnt, size_t *outn) {
    unsigned long long rb = ZSTD_getFrameContentSize(b, n);
    if (rb == ZSTD_CONTENTSIZE_UNKNOWN || rb == ZSTD_CONTENTSIZE_ERROR)
        die("int names have no usable frame size");
    uint8_t *v = xmalloc((size_t)rb);
    if (ZSTD_isError(ZSTD_decompress(v, rb, b, n))) die("int names decompress failed");
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

#define TOKMAX 9000000u

static uint32_t TOKCHUNK = 3000000;     /* --tok-chunk overrides */
static int TOKSET = 0, V1 = 0;
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
    int nthr = nthreads();
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
static size_t tokbest(buf *in, uint8_t **out) {
    if (TOKSET) return tokc(in, out, TOKCHUNK);
    size_t n = sizeof TCAND / sizeof *TCAND, best = 0;
    for (size_t i = 0; i < n; i++) {
        if (FASTMODE && i) continue;    /* fast: the smallest chunk only; --slow tries all */
        uint8_t *o2 = NULL;
        size_t b = tokc(in, &o2, TCAND[i]);
        /* a larger chunk must save 0.1% to win: chunks decode in parallel */
        if (b && (!best || b + b / 1000 < best)) { free(*out); *out = o2; best = b; }
        else free(o2);
    }
    return best;
}

/* smallest of zstd, integer deltas, tokeniser; losers freed; kind 0/1/2 */
static size_t bestnames(buf *b, uint64_t n, uint8_t **out, int *kind) {
    uint8_t *a1 = NULL, *a2 = NULL, *a3 = NULL;
    size_t s1 = 0, s2 = 0, s3 = 0;
    double _t0 = nowsec();
    /* zstd and the tokeniser are independent and neither gains from threads
       on a buffer this size, so run them side by side. tok3 overwrites every
       separator with NUL and puts it back, so it cannot share the buffer with
       a concurrent reader: it gets a copy. b is at most the 200,000-name
       sample here, one tokeniser chunk, so nothing is lost by not nesting. */
    buf tb = { xmalloc(b->n ? b->n : 1), b->n, b->n };
    memcpy(tb.d, b->d, b->n);
#ifdef _OPENMP
#pragma omp parallel sections num_threads(2)
#endif
    {
#ifdef _OPENMP
#pragma omp section
#endif
        { s1 = zc(b, &a1); s2 = intnames(b, n, &a2); }
#ifdef _OPENMP
#pragma omp section
#endif
        { s3 = tokbest(&tb, &a3); }
    }
    free(tb.d);
    if (PROF) fprintf(stderr, "  names %zu B: candidates %.3f s\n", b->n, nowsec() - _t0);
    size_t best = s1; *kind = 0; *out = a1;
    if (s2 && s2 < best) { best = s2; *kind = 1; *out = a2; }
    if (s3 && s3 < best) { best = s3; *kind = 2; *out = a3; }
    if (*out != a1) free(a1);
    if (*out != a2) free(a2);
    if (*out != a3) free(a3);
    return best;
}


static int nosz(int c) { unsigned u = (unsigned)c; return u && !(u & CODEC_BZIP3) && (u & 0x10); }
static int badcodec(int c) {   /* bits no codec defines, or a bzip3 block size out of range */
    unsigned u = (unsigned)c;
    if (u & CODEC_BZIP3) { unsigned mb = u & 0xffff;
                           return (u & ~(unsigned)(CODEC_BZIP3 | 0xffff)) || mb < 1 || mb > 511; }
    return (u & ~(unsigned)(CODEC_ARITH | 0x3ff)) != 0;
}
static int noext(int c) { unsigned u = (unsigned)c; return (u & CODEC_ARITH) && !(u & CODEC_BZIP3) && (u & 0x04); }

static size_t paycomp(const uint8_t *src, size_t n, int codec, uint8_t **out) {
    unsigned int ol;
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
static int TAILFORCE = -1;       /* --tail overrides the decision */
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
            if (TAILFORCE >= 0) TAILON = TAILFORCE;
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
typedef struct { char *name; uint64_t *st, end; uint32_t *lab; size_t n, cap; int64_t ri; _Atomic int filled; } rcon;

static qseg *QS; static size_t nqs, cqs;
static rcon *RC; static size_t nrc, crc;

static uint64_t ALNHASH;
static int REFLAZY;                  /* readers: decode a reference contig when first asked */
static void loadref_kst(const char *path);
static void reffill(rcon *rc);
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
    RC[nrc] = (rcon){ xstrdup(nm), NULL, 0, NULL, 0, 0, -1, 1 };
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
        if (ln[n - 1] == '\n') ln[--n] = 0;
        if (n && ln[n - 1] == '\r') ln[--n] = 0;
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
    free(ln); fclose(f); ALNHASH ^= rh;
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
        if (ln[n - 1] == '\n') ln[--n] = 0;
        if (n && ln[n - 1] == '\r') ln[--n] = 0;

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
#define ALNMAGIC "KSTAL1\n"   /* carries the reference contig-order hash */
#define ALNMAGLEN 6
#define ALNTOK 0x7fff0000   /* the name stream went through the tokeniser */

static uint64_t refhash(void) {   /* the reference contig names, in order */
    uint64_t h = FNV0;
    for (size_t i = 0; i < nrc; i++) { h = fnvbuf(h, RC[i].name, strlen(RC[i].name)); h = fnvbuf(h, "\n", 1); }
    return h;
}

static void alnpack(const char *inpaf, const char *out) {
    loadpaf(inpaf);
    buf names = {0}, cnt = {0}, pos = {0}, len = {0}, seg_ = {0};
    for (size_t i = 0; i < nqs; i++) {
        const qseg *Q = QS + i;
        if (!Q->n) continue;   /* no aligned base: nothing to pack */
        bput(&names, Q->name, strlen(Q->name)); bput(&names, "\n", 1);
        bvar(&cnt, Q->n);

        uint64_t pq = 0, pt = 0; uint32_t px = 0; int prev_rev = -1; uint64_t run = 0;
        for (size_t j = 0; j < Q->n; j++) {
            const seg *g = Q->s + j;
            if (!j) { bvar(&pos, g->q); bvar(&pos, g->tp); }
            else if (g->t != px || g->rev != prev_rev) {
                bvar(&pos, g->q); bvar(&pos, g->tp);       /* absolute again */
            } else {
                bvar(&pos, zz((int64_t)(g->q - pq)));
                bvar(&pos, zz((int64_t)(g->tp - pt)));
            }
            if (g->t != px || g->rev != prev_rev) {
                if (j) bvar(&seg_, run);
                bvar(&seg_, (uint64_t)g->t * 2 + (g->rev ? 1 : 0));
                px = g->t; prev_rev = g->rev; run = 0;
            }
            run++;
            bvar(&len, g->n);
            pq = g->q + g->n; pt = g->tp + g->n;
        }
        bvar(&seg_, run);                       /* close the last run */
    }
    /* each stream gets the codec search */
    buf all = {0};
    buf *st[5] = { &names, &cnt, &pos, &len, &seg_ };
    static const char *SNAME[5] = { "names", "counts", "positions", "lengths", "targets" };
    for (int k = 0; k < 5; k++) {
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
    free(names.d); free(cnt.d); free(pos.d); free(len.d); free(seg_.d); free(all.d);
}

static int alnload(const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f) return 0;
    char mg[ALNMAGLEN];
    if (fread(mg, 1, ALNMAGLEN, f) != ALNMAGLEN || memcmp(mg, ALNMAGIC, ALNMAGLEN)) {
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
    fclose(f);
    const uint8_t *p = d, *e = d + n;
    const uint8_t *sec[5]; size_t slen[5]; uint8_t *raw[5] = {0};
    for (int k = 0; k < 5; k++) {
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
    const char *nm = (const char *)sec[0], *nme = nm + slen[0];
    while (nm < nme) {
        const char *nl = memchr(nm, '\n', (size_t)(nme - nm));
        if (!nl) die("corrupt packed alignment: a name has no newline");
        char *name = xmalloc((size_t)(nl - nm) + 1);
        memcpy(name, nm, (size_t)(nl - nm)); name[nl - nm] = 0;
        qseg *Q = qseg_get(name, 1);
        free(name);
        uint64_t cn = vgetb(&cp, cpe);
        if (cn > (uint64_t)(lpe - lp)) die("corrupt packed alignment: more blocks than lengths");
        Q->s = xrealloc(Q->s, mul((size_t)cn, sizeof *Q->s)); Q->cap = cn; Q->n = 0;
        uint64_t pq = 0, pt = 0; uint32_t tx = 0; int rev = 0, fresh = 1;
        uint64_t left = 0;
        for (uint64_t j = 0; j < cn; j++) {
            uint64_t q, tp;
            if (!left) {                          /* next target/strand run */
                uint64_t xr = vgetb(&sp, spe);
                rev = (int)(xr & 1); tx = (uint32_t)(xr >> 1);
                if (tx >= nrc) die("corrupt packed alignment: a target is not in the reference");
                left = vgetb(&sp, spe);
                if (!left) die("corrupt packed alignment: an empty target run");
                fresh = 1;
            }
            if (fresh) { q = vgetb(&pp, ppe); tp = vgetb(&pp, ppe); fresh = 0; }
            else {
                q  = pq + (uint64_t)unzz(vgetb(&pp, ppe));
                tp = pt + (uint64_t)unzz(vgetb(&pp, ppe));
            }
            uint64_t ln = vgetb(&lp, lpe);
            Q->s[Q->n++] = (seg){ q, tx, (uint32_t)ln, tp, rev };
            pq = q + ln; pt = tp + ln; left--;
        }
        nm = nl + 1;
    }
    free(d);
    for (int k = 0; k < 5; k++) free(raw[k]);
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
    if (!f) die("cannot open the file given to --names");
    char *ln = NULL; size_t cap = 0; ssize_t n; size_t c = 0;
    while ((n = getline(&ln, &cap, f)) > 0) {
        if (ln[n - 1] == '\n') ln[--n] = 0;
        if (n && ln[n - 1] == '\r') ln[--n] = 0;
        if (!c && (ln[0] == '>' || ln[0] == '@')) die("--names on encode takes a plain list, one name a line");
        if (c == nord) { nord = nord ? nord * 2 : 1 << 16;
                         ORD = xrealloc(ORD, nord * sizeof *ORD); }
        ORD[c++] = xstrdup(ln);
    }
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
static size_t runat(const rcon *rc, uint64_t tpos) {   /* count of starts <= tpos */
    if (tpos >= rc->end) return 0;      /* past the last run: no answer, not the last run's */
    size_t a;
    if (rc == CR) {
        a = CRI;
        int k = 0;
        while (a < rc->n && rc->st[a] <= tpos && k++ < 8) a++;
        while (a > 0 && rc->st[a - 1] > tpos && k++ < 16) a--;
        if (k >= 16 || (a < rc->n && rc->st[a] <= tpos)) a = rc->n + 1;
    } else a = rc->n + 1;
    if (a == rc->n + 1) {
        size_t b = rc->n; a = 0;
        while (a < b) { size_t m = (a + b) / 2;
                        if (rc->st[m] <= tpos) a = m + 1; else b = m; }
    }
    CR = rc; CRI = a;
    return a;
}
static uint32_t reflabel(const qseg *Q, uint64_t pos) {
    const seg *g = segat(Q, pos);
    if (!g) return UINT32_MAX;
    uint64_t tpos = g->rev ? g->tp + (g->q + g->n - 1 - pos) : g->tp + (pos - g->q);
    rcon *rc = RC + g->t; reffill(rc);
    size_t a = runat(rc, tpos);
    return a ? rc->lab[a - 1] : UINT32_MAX;
}

static uint64_t refend(const qseg *Q, uint64_t pos) {
    const seg *g = segat(Q, pos);
    LASTREV = g ? g->rev : 0;
    if (!g) return 0;
    uint64_t tpos = g->rev ? g->tp + (g->q + g->n - 1 - pos) : g->tp + (pos - g->q);
    rcon *rc = RC + g->t; reffill(rc);
    size_t a = runat(rc, tpos);
    if (!a) return 0;
    uint64_t tend = a < rc->n ? rc->st[a] : rc->end;   /* target end of that run */
    if (g->rev) { uint64_t ts = rc->st[a - 1];
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
        if (n && ln[n - 1] == '\n') ln[--n] = 0;
        if (n && ln[n - 1] == '\r') ln[--n] = 0;
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
        if (len && line[len - 1] == '\n') line[--len] = 0;
        if (len && line[len - 1] == '\r') line[--len] = 0;
        if (!len) continue;
        char *t1 = memchr(line, '\t', (size_t)len); if (!t1) continue;
        char *t2 = memchr(t1 + 1, '\t', (size_t)(len - (t1 + 1 - line))); if (!t2) continue;
        if (memchr(t2 + 1, '\t', (size_t)(len - (t2 + 1 - line)))) die("hierarchy: a row has more than three columns");
        if (!strncmp(line, "feature_set\t", 12)) continue;
        for (size_t k = 0; k < nfset; k++)
            if ((size_t)(t1 - line) == strlen(fset[k]) && !memcmp(line, fset[k], strlen(fset[k]))) {
                vocadd(k, t1 + 1, (size_t)(t2 - t1 - 1)); vocadd(k, t2 + 1, (size_t)(len - (t2 + 1 - line))); }
    }
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

static int srec_cmp(const void *a, const void *b) {
    const srec *x = a, *y = b;
    uint32_t n = x->nrun < y->nrun ? x->nrun : y->nrun;
    const srun *p = rrun + x->off, *q = rrun + y->off;
    for (uint32_t i = 0; i < n; i++)                 /* labels lead */
        if (p[i].lab != q[i].lab) return p[i].lab < q[i].lab ? -1 : 1;
    if (x->nrun != y->nrun) return x->nrun < y->nrun ? -1 : 1;
    for (uint32_t i = 0; i < n; i++)                 /* then lengths */
        if (p[i].len != q[i].len) return p[i].len < q[i].len ? -1 : 1;

    if (x->nam != y->nam) return x->nam < y->nam ? -1 : 1;
    return 0;
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
    for (size_t i = 0; i < rrecn; i++) {
        char ib[24]; const char *nm = rbuf + rrec[i].nam; size_t nl2;
        ALNNAME = nm;
        if (PERMON) { nl2 = (size_t)sprintf(ib, "%u", rrec[i].idx); nm = ib; }
        else nl2 = strlen(nm);
        uint64_t pos = 0;
        for (uint32_t k = 0; k < rrec[i].nrun; k++) {
            const srun *u = rrun + rrec[i].off + k;
            const char *lb = labtab[u->lab];
            ktw_row(w, nm, nl2, u->len, lb, strlen(lb), pos, 0);
            pos += u->len;
        }
    }
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
    if (r->nrun > r->cap) { r->cap = r->nrun; r->rs = xrealloc(r->rs, r->cap * sizeof *r->rs);
                            if (!r->rs) die("oom"); }
    if (fread(r->rs, sizeof(srun), r->nrun, r->f) != r->nrun) die("sort temp is corrupt");
    return 1;
}

static int srd_cmp(const srun_rd *a, const srun_rd *b) {
    uint32_t n = a->nrun < b->nrun ? a->nrun : b->nrun;
    for (uint32_t i = 0; i < n; i++)
        if (a->rs[i].lab != b->rs[i].lab) return a->rs[i].lab < b->rs[i].lab ? -1 : 1;
    if (a->nrun != b->nrun) return a->nrun < b->nrun ? -1 : 1;
    for (uint32_t i = 0; i < n; i++)
        if (a->rs[i].len != b->rs[i].len) return a->rs[i].len < b->rs[i].len ? -1 : 1;

    if (a->idx != b->idx) return a->idx < b->idx ? -1 : 1;
    return 0;
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
        char ib[24]; const char *nmp = r->nm; size_t nl;
        ALNNAME = r->nm;
        if (PERMON) { nl = (size_t)sprintf(ib, "%u", r->idx); nmp = ib; }
        else nl = strlen(r->nm);
        uint64_t pos = 0;
        for (uint32_t k = 0; k < r->nrun; k++) {
            const char *lb = labtab[r->rs[k].lab];
            ktw_row(w, nmp, nl, r->rs[k].len, lb, strlen(lb), pos, 0);
            pos += r->rs[k].len;
        }
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


static int fastenough(int codec) {
    return (codec & CODEC_BZIP3) == 0;
}

static void pickcodecs(int rawfd, const blkmeta *bm, size_t nraw) {

    if (IDCODEC >= 0 && LENCODEC >= 0) return;

        static const int CAND_FAST[] = { 0x100, 1, 65, 521, CODEC_ARITH, CODEC_ARITH | 1,
                                         CODEC_BZIP3 | 16, CODEC_BZIP3 | 64 };

        /* rANS order 0/1 x {PACK,RLE,X32,STRIPE}, arith order 0/1 x {PACK,RLE},
           bzip3 at six block sizes. Was built by three loops at every call. */
        static const int CAND_SLOW[] = {
            0x100, 0x80, 0x40, 0xc0, 0x8, 0x88, 0x48, 0xc8,
            0x4, 0x84, 0x44, 0xc4, 0xc, 0x8c, 0x4c, 0xcc,
            0x1, 0x81, 0x41, 0xc1, 0x9, 0x89, 0x49, 0xc9,
            0x5, 0x85, 0x45, 0xc5, 0xd, 0x8d, 0x4d, 0xcd,
            CODEC_ARITH, CODEC_ARITH | 0x80, CODEC_ARITH | 0x40, CODEC_ARITH | 0xc0,
            CODEC_ARITH | 1, CODEC_ARITH | 0x81, CODEC_ARITH | 0x41, CODEC_ARITH | 0xc1,
            CODEC_BZIP3 | 8, CODEC_BZIP3 | 16, CODEC_BZIP3 | 32,
            CODEC_BZIP3 | 64, CODEC_BZIP3 | 128, CODEC_BZIP3 | 256,
        };
        const int *CAND = FASTMODE ? CAND_FAST : CAND_SLOW;
        const size_t NCAND = (FASTMODE ? sizeof CAND_FAST : sizeof CAND_SLOW) / sizeof *CAND;

        size_t rawtot = 0, maxstream = 0;
        for (size_t i2 = 0; i2 < nraw; i2++) { rawtot += bm[i2].len;
            if (bm[i2].lens > maxstream) maxstream = bm[i2].lens;
            if (bm[i2].len - bm[i2].lens > maxstream) maxstream = bm[i2].len - bm[i2].lens; }
        int big = maxstream >= 0x7fffffffu;   /* rANS and arith decode below INT_MAX */
        int seekable = nraw > 1 || (uint64_t)rawtot > BLOCK;
        size_t SCAP = seekable ? (4u << 20) : (96u << 20);

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
           and changes what --slow means. Measured and declined. */
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
            int z22 = SLOW && !LEVELSET;
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
                if (!((SLOW || !seekable || !FASTMODE || fastenough(CAND[k]))
                      && (!big || (CAND[k] & CODEC_BZIP3)))) {
                    ci[k] = cl[k] = (size_t)-1 / 2; continue;
                }
                /* Size alone does not qualify a codec: */
                uint8_t *t = NULL; size_t g;
                if (ni && (g = paycomp(ids, ni, CAND[k], &t))) {
                    uint8_t *v = malloc(ni);
                    if (v && paydecomp(t, g, CAND[k], v, ni) == ni
                        && !memcmp(v, ids, ni)) ci[k] += g;
                    else ci[k] = (size_t)-1 / 2;
                    free(v);
                } else ci[k] = (size_t)-1 / 2;
                free(t); t = NULL;
                if (nl && (g = paycomp(src, nl, CAND[k], &t))) {
                    uint8_t *v = malloc(nl);
                    if (v && paydecomp(t, g, CAND[k], v, nl) == nl
                        && !memcmp(v, src, nl)) cl[k] += g;
                    else cl[k] = (size_t)-1 / 2;
                    free(v);
                } else cl[k] = (size_t)-1 / 2;
                free(t);
            }
            if (ZSTD_isError(zb[0][0]) || ZSTD_isError(zb[0][1])
                || ZSTD_isError(zb[1][0]) || ZSTD_isError(zb[1][1]))
                die("sample compress failed");
            zi += zb[0][0]; zl += zb[0][1];
            z22i += zb[1][0]; z22l += zb[1][1];
        }

        int wi = 0, wl = 0;
        size_t bi = zi, bl2 = zl;
        if (SLOW && !LEVELSET && z22i < bi) { bi = z22i; }
        if (SLOW && !LEVELSET && z22l < bl2) { bl2 = z22l; }
        int lvl_i = (SLOW && !LEVELSET && z22i < zi) ? 22 : LEVEL;
        int lvl_l = (SLOW && !LEVELSET && z22l < zl) ? 22 : LEVEL;
        for (unsigned k = 0; k < NCAND; k++) {
            int ok = (SLOW || !seekable || !FASTMODE || fastenough(CAND[k])) && (!big || (CAND[k] & CODEC_BZIP3));
            size_t need_i = SLOW ? bi : (size_t)(0.995 * (double)zi);
            size_t need_l = SLOW ? bl2 : (size_t)(0.995 * (double)zl);
            if (zi && ci[k] < bi && ci[k] < need_i + (SLOW ? 1 : 0) && ok)
                { bi = ci[k]; wi = CAND[k]; }
            if (zl && cl[k] < bl2 && cl[k] < need_l + (SLOW ? 1 : 0) && ok)
                { bl2 = cl[k]; wl = CAND[k]; }
        }

        if (SLOW && !LEVELSET && !wi && !wl) LEVEL = (lvl_i == 22 || lvl_l == 22) ? 22 : LEVEL;
        else if (SLOW && !LEVELSET && (!wi || !wl))
            LEVEL = (!wi ? lvl_i : lvl_l);
        /* --slow decides on the whole file: */
        if (SLOW) {
            size_t *ti = xcalloc(NCAND + 2, sizeof *ti), *tl = xcalloc(NCAND + 2, sizeof *tl);
            int *live = xcalloc(NCAND + 2, sizeof *live);
            for (unsigned k = 0; k < NCAND; k++)
                live[k] = ((ci[k] < (size_t)-1 / 4) || (cl[k] < (size_t)-1 / 4)) && (!big || (CAND[k] & CODEC_BZIP3));
            live[NCAND] = 1; live[NCAND + 1] = !LEVELSET;
            for (size_t b2 = 0; b2 < nraw; b2++) {
                if (bm[b2].len > sbcap) { sbcap = bm[b2].len; sbuf = xrealloc(sbuf, sbcap); }
                rawget(rawfd, (off_t)bm[b2].off, sbuf, bm[b2].len);
                const uint8_t *src = sbuf;
                size_t nl = bm[b2].lens, ni = bm[b2].len - nl;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads())
#endif
                for (unsigned k = 0; k < NCAND + 2; k++) {
                    if (!live[k]) continue;
                    uint8_t *t = NULL; size_t g;
                    if (k < NCAND) {
                        if (ni) { g = paycomp(src + nl, ni, CAND[k], &t); ti[k] += g ? g : (size_t)-1 / 4; free(t); t = NULL; }
                        if (nl) { g = paycomp(src, nl, CAND[k], &t); tl[k] += g ? g : (size_t)-1 / 4; free(t); }
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
                    if (!live[k]) continue;
                    if (ti[k] < ti[qi]) qi = k;
                    if (tl[k] < tl[ql]) ql = k;
                }
                int bad_i = 0, bad_l = 0;
                for (size_t b2 = 0; b2 < nraw && !(bad_i && bad_l); b2++) {
                    if (bm[b2].len > sbcap) { sbcap = bm[b2].len; sbuf = xrealloc(sbuf, sbcap); }
                    rawget(rawfd, (off_t)bm[b2].off, sbuf, bm[b2].len);
                    const uint8_t *src = sbuf;
                    size_t nl = bm[b2].lens, ni = bm[b2].len - nl;
                    uint8_t *t = NULL, *v = NULL; size_t g;
                    if (qi < NCAND && ni && !bad_i) {
                        g = paycomp(src + nl, ni, CAND[qi], &t); v = g ? malloc(ni) : NULL;
                        if (!g || !v || paydecomp(t, g, CAND[qi], v, ni) != ni || memcmp(v, src + nl, ni)) bad_i = 1;
                        free(t); free(v); t = v = NULL; }
                    if (ql < NCAND && nl && !bad_l) {
                        g = paycomp(src, nl, CAND[ql], &t); v = g ? malloc(nl) : NULL;
                        if (!g || !v || paydecomp(t, g, CAND[ql], v, nl) != nl || memcmp(v, src, nl)) bad_l = 1;
                        free(t); free(v); }
                }
                if (!bad_i && !bad_l) break;
                if (bad_i) live[qi] = 0;
                if (bad_l) live[ql] = 0;
            }
            wi = qi < NCAND ? CAND[qi] : 0;
            wl = ql < NCAND ? CAND[ql] : 0;
            if (!wi) LEVEL_I = qi == NCAND + 1 ? 22 : LEVEL;
            if (!wl) LEVEL_L = ql == NCAND + 1 ? 22 : LEVEL;
            fprintf(stderr, "codec: --slow ran %zu live candidates over %zu blocks whole\n",
                    (size_t)NCAND + 2, nraw);
            free(ti); free(tl); free(live);
        }
        free(sbuf);
        if (IDCODEC < 0) IDCODEC = wi;
        if (LENCODEC < 0) LENCODEC = wl;
        fprintf(stderr, "codec: len=0x%x id=0x%x from %zu sample blocks, %s\n",
                LENCODEC, IDCODEC, np, seekable ? "seekable" : "one block");
    if (IDCODEC < 0) IDCODEC = 0;
    if (LENCODEC < 0) LENCODEC = 0;
}

static size_t complens(ZSTD_CCtx *tc, const uint8_t *src, size_t n, uint8_t **out, size_t bound) {
    if (LENCODEC) return paycomp(src, n, LENCODEC, out);
    *out = malloc(bound + 64);
    if (!*out) return 0;
    size_t a = ZSTD_compressCCtx(tc, *out, bound + 64, src, n, LEVEL_L < 0 ? LEVEL : LEVEL_L);
    if (ZSTD_isError(a)) { free(*out); *out = NULL; return 0; }
    return a;
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

static void ktw_close(ktw *w) {
    buf names = w->names, sizes = w->sizes;
    size_t nidx = w->nidx;
    uint64_t nblk = w->nblk, rows = w->rows;
    FILE *tmp = w->tmp;
    const blkmeta *bm = w->bm; size_t nraw = w->nraw;

    double _ts = nowsec();
    pickcodecs(w->rawfd, bm, nraw);
    PH_ADD(8, _ts);

    size_t maxraw = 0;
    for (size_t i = 0; i < nraw; i++) {
        if (bm[i].len > maxraw) maxraw = bm[i].len;
        for (int v = 0; ALNPAF && v < 2; v++)
            if (bm[i].alen[v] + (bm[i].len - bm[i].lens) > maxraw) maxraw = bm[i].alen[v] + (bm[i].len - bm[i].lens);
    }
    size_t obound = ZSTD_compressBound(maxraw ? maxraw : 1);

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
                static const int ORDER[3] = { 1, 2, 0 };
                size_t best = 0; int bv = -1;
                for (int k = 0; k < 3; k++) {
                    int v = ORDER[k]; if (!vs[v]) continue;
                    uint8_t *f2 = NULL; size_t a2;
                    if (FASTMODE && ALNPAF) {
                        a2 = ZSTD_compressCCtx(tc, o, 2 * obound + 64, vs[v], vn[v], 1);
                        if (ZSTD_isError(a2)) { failed = 1; break; }
                    } else {
                        a2 = complens(tc, vs[v], vn[v], &f2, ZSTD_compressBound(vn[v]));
                        if (!a2) { failed = 1; break; }
                    }
                    if (bv < 0 || a2 + (FASTMODE ? a2 / 100 : 0) < best) { free(fa); fa = f2; best = a2; bv = v; }
                    else free(f2);
                }
                if (!failed && FASTMODE && ALNPAF) { best = complens(tc, vs[bv], vn[bv], &fa, ZSTD_compressBound(vn[bv])); if (!best) failed = 1; }
                if (failed) { free(o); free(fa); continue; }
                a = best; nl = vn[bv]; nosent = (uint64_t)bv;
                if (IDCODEC) { b2 = paycomp(idsrc, ni, IDCODEC, &fb);
                               if (!b2) { failed = 1; free(fa); free(o); continue; } }
                else
                { b2 = ZSTD_compressCCtx(tc, o, obound, idsrc, ni, LEVEL_I < 0 ? LEVEL : LEVEL_I);
                  if (ZSTD_isError(b2)) { failed = 1; free(fa); free(o); continue; } }
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
    if (NAMEMODE == 2) ordinal = 0;                    /* keep: store them */
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
            "file beside the track would be redundant; use --names keep to store them");

    if (PERMON && !NATORD && ordinal && NAMEMODE != 1) {
        free(names.d); names = natbuf; natbuf.n = 0; natbuf.d = NULL;
        ordinal = 0; PERMON = 0;
    }
    if (ordinal) { names.n = 0; NAMES_OUT = NULL; }

    size_t nzl = 0; int tokused = 0, intused = 0;
    double _tn = nowsec();
    if (!ordinal) {
        buf smp = { names.d, 0, 0 }; uint64_t smpn = 0;   /* a view, not a copy */
        { const char *q = (const char *)names.d, *e2 = q + names.n;
          while (q < e2 && smpn < 200000) {
              const char *nl = memchr(q, '\n', (size_t)(e2 - q));
              if (!nl) break;
              q = nl + 1; smpn++; }
          smp.n = smp.cap = (size_t)((const uint8_t *)q - names.d); }
        uint8_t *smpz = NULL; int pick;               /* rank on the sample, run the winner */
        size_t smpl = bestnames(&smp, smpn, &smpz, &pick);
        if (smp.n == names.n) {                       /* the sample was the whole
                                                         buffer: the winner is done */
            nz = smpz; nzl = smpl;
            intused = pick == 1; tokused = pick == 2;
        } else {
            free(smpz);
            if (pick == 1) { nzl = intnames(&names, nidx, &nz); intused = nzl != 0; }
            else if (pick == 2) { nzl = tokbest(&names, &nz); tokused = nzl != 0; }
        }
        if (!nzl) { nzl = zc(&names, &nz); tokused = intused = 0; }  /* refused */
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
            uint8_t *bp = NULL; int bpick;
            size_t best = bestnames(&srt, nidx, &bp, &bpick);
            uint8_t *na = NULL, *nb2 = NULL;
            size_t p1 = zc(&natbuf, &na), p2 = tokbest(&natbuf, &nb2);
            size_t pbest = p1; uint8_t *pp = na; int ptok = 0;
            if (p2 && p2 < pbest) { pbest = p2; pp = nb2; ptok = 1; }
            if (best && best <= nzl + pbest) {     /* sorted names win: no rank */
                free(nz); nz = bp; nzl = best;
                tokused = bpick == 2; intused = bpick == 1;
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
    uint64_t nichunk = 0;
    {
        const uint8_t *q = ib.d, *qe = ib.d + ib.n;
        uint64_t cpend = 0, seen = 0;
        const uint8_t *cstart = q;
        uint64_t spend = 0;
        while (q < qe) {
            uint64_t nruns = vget(&q);
            cpend += nruns;
            if (++seen % ICHUNK == 0 || q >= qe) {
                buf one = { (uint8_t *)cstart, (size_t)(q - cstart), 0 };
                uint8_t *cz = NULL;
                size_t cl = zcl(&one, &cz, METALEVEL);
                bput(&ichunks, cz, cl); free(cz);
                bvar(&itab, cl); bvar(&itab, spend);
                nichunk++;
                cstart = q; spend = cpend;
            }
        }
    }
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
    uint64_t h[16] = { nidx, (NAMES_OUT ? 0 : nzl) | (tokused ? TOKFLAG : 0) | (CARTUSED ? CARTFLAG : 0)
                             | (intused ? INTFLAG : 0) | (dropnames && !(w->ordok && nidx > 0) ? DROPPED : 0)
                             | (ALNPAF ? LENFLAG : 0) | (TAILON > 0 ? TAILFLAG : 0)
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
            nzl, dropnames ? " (dropped, needs --names FILE)" : ordinal ? " (ordinal, not stored)" : tokused ? " (tokenised)" : intused ? " (int deltas)" : "", izl, szl,
            (unsigned long long)BYT_LENS, (unsigned long long)BYT_LABS, RAWHELD);

    free(itab.d); free(ichunks.d);
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
kt_writer *kt_open(FILE *out) { return ktw_open(out); }
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
    uint8_t *nz_keep; size_t nzl_keep; int names_ready, names_tok, names_int, dropped;
    uint8_t *ibuf; int idcodec, lencodec;
    const uint8_t *cp, *cpe;   /* cpe: end of the decoded index chunk */
    uint64_t ci, cpend; int cvalid; rec crec;
    off_t nz_off; uint64_t mhash, feathash; int cart;
    uint64_t fsz; int hasidx;   /* data size, and whether an index trailer follows it */
    uint8_t *pz; size_t pzl; int perm_tok, perm_ext;  /* names in input order */
    char *namebuf, *natbufp;
    int name_borrowed;         /* the table belongs to another reader */
    uint8_t *ichunk; uint64_t nichunk, *coff, *cpendv, curchunk, ibcap;
    uint64_t *sst, *sen; int *sfeat; uint8_t *snov; uint64_t scap;
    const uint8_t *curlp, *curip;
    uint64_t rowhash, tailpred, tailskip; int taildec;   /* for a track read as --ref */
} R;

static void loadichunk(R *r, uint64_t c) {
    if (c >= r->nichunk) die("index chunk out of range");
    size_t cl = (size_t)(r->coff[c + 1] - r->coff[c]);
    const uint8_t *src = r->ichunk + r->coff[c];
    unsigned long long rb = ZSTD_getFrameContentSize(src, cl);
    if (rb == ZSTD_CONTENTSIZE_UNKNOWN || rb == ZSTD_CONTENTSIZE_ERROR)
        die("index chunk has no usable frame size");
    if (rb > r->ibcap) {
        r->ibuf = xrealloc(r->ibuf, (size_t)rb);
        r->ibcap = rb;
    }
    if (ZSTD_isError(ZSTD_decompress(r->ibuf, rb, src, cl)))
        die("index chunk decompress failed");
    r->curchunk = c;
    r->cp = r->ibuf; r->cpe = r->ibuf + (size_t)rb;
    r->ci = c * ICHUNK; r->cpend = r->cpendv[c];

    {   const uint8_t *q = r->cp; uint64_t want = c + 1 < r->nichunk ? ICHUNK : r->nreads - c * ICHUNK, sum = 0, k = 0;
        for (; k < want && q < r->cpe; k++) { uint64_t v = vgetb(&q, r->cpe); if (v > UINT32_MAX) die("corrupt index: a run count past 32 bits"); sum += v; }
        if (k != want || q != r->cpe || r->cpendv[c] + sum != r->cpendv[c + 1])
            die("corrupt index: an index chunk does not match its table entry"); }
}

static rec getidx(R *r, uint64_t i) {
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
    if (fsz >= 32 + 8 + 128) {
        uint64_t t[4];
        if (fseeko(r->f, (off_t)fsz - 32, SEEK_SET) || fread(t, 8, 4, r->f) != 4) die("seek");
        if (t[3] == IDXMAGIC && t[0] == h[0] && t[0] <= (fsz - 32) / 20
            && fsz - 32 - t[0] * 20 >= 8 + 128) { fsz -= 32 + t[0] * 20; r->hasidx = 1; }
    }
    r->fsz = fsz;
    ALLOCMAX = (size_t)1 << 38;
    if (fseeko(r->f, 8 + 128, SEEK_SET)) die("seek");
    uint64_t nml = h[1] & ~(TOKFLAG | INTFLAG | DROPPED | ALNFLAG | DBFLAG
                            | SIDEFLAG | LENFLAG | TAILFLAG | PERMFLAG
                            | PERMTOK | CARTFLAG);   /* h[1] also holds flags */

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
        unsigned long long lb = ZSTD_getFrameContentSize(lz, lzl);
        if (lb == ZSTD_CONTENTSIZE_UNKNOWN || lb == ZSTD_CONTENTSIZE_ERROR)
            die("label table has no usable frame size");
        char *lbuf = xmalloc((size_t)lb + 1);
        if (ZSTD_isError(ZSTD_decompress(lbuf, lb, lz, lzl)))
            die("label table decompress failed");
        lbuf[lb] = 0;

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
    if (names_tok && names_int) die("corrupt file: names flagged both tokenised and integer");
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
            die("this track holds featureIDs: decode it with --db features.tsv --out DIR, "
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
    TAILPRED_R = h[12]; TAILSKIP_R = h[13] & 0xff;
    if (TAILDEC && TAILSKIP_R != TAILSKIP) die("corrupt file: tail skip is not the format constant");
    r->taildec = TAILDEC; r->tailpred = h[12]; r->tailskip = TAILSKIP_R; r->rowhash = was_aln ? 0 : h[13] >> 8;
    if (was_aln && ALNPAF && (h[13] >> 8) != (ALNHASH & ((1ull << 56) - 1)))
        die("this PAF, reference or name map is not the one the track was encoded with");
    r->cart = (h[1] & CARTFLAG) != 0;
    if (r->cart && !was_db) die("corrupt file: vocabulary featureIDs without --db");
    h[1] &= ~(TOKFLAG | INTFLAG | DROPPED | ALNFLAG | DBFLAG | SIDEFLAG | LENFLAG
              | TAILFLAG | PERMFLAG | PERMTOK | CARTFLAG);
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

    r->nz_keep = sidebuf ? sidebuf : nz; r->nzl_keep = nzl_eff; r->names_ready = 0; r->names_tok = names_tok; r->names_int = names_int; r->dropped = dropped;
    r->name = NULL;                    /* 160 MB, allocated by ensurenames */

    /* Walk varints with a cursor, allocate nothing. */

    r->nichunk = h[10];
    r->ichunk = iz;                      /* the concatenated chunk frames */
    r->curchunk = UINT64_MAX;
    {   uint8_t *itz = xmalloc(h[11]);
        if (fread(itz, 1, h[11], r->f) != h[11]) die("short index table");
        unsigned long long tb = ZSTD_getFrameContentSize(itz, h[11]);
        if (tb == ZSTD_CONTENTSIZE_UNKNOWN || tb == ZSTD_CONTENTSIZE_ERROR)
            die("index table has no usable frame size");
        uint8_t *tbuf = xmalloc(tb);
        if (ZSTD_isError(ZSTD_decompress(tbuf, tb, itz, h[11])))
            die("index table decompress failed");
        free(itz);
        if (r->nichunk == SIZE_MAX) die("corrupt file: chunk count");
        r->coff = xmalloc(mul(r->nichunk + 1, 8));
        r->cpendv = xmalloc(mul(r->nichunk + 1, 8));
        if (!r->coff || !r->cpendv) die("oom");
        const uint8_t *q = tbuf, *qe = tbuf + tb;
        uint64_t acc = 0, mx = 0;
        for (uint64_t c = 0; c < r->nichunk; c++) {
            uint64_t cl = vgetb(&q, qe);
            r->cpendv[c] = vgetb(&q, qe);
            r->coff[c] = acc; acc += cl;
            if (cl > mx) mx = cl;
        }
        r->coff[r->nichunk] = acc;
        if (acc != h[2]) die("index chunk table does not match the index");
        if (r->nichunk != (r->nreads + ICHUNK - 1) / ICHUNK) die("corrupt file: index chunks do not match the sequence count");
        r->cpendv[r->nichunk] = r->rows;
        for (uint64_t c = 0; c < r->nichunk; c++)
            if (r->cpendv[c] > r->cpendv[c + 1] || (c == 0 && r->cpendv[0])) die("corrupt file: index chunk starts are not in order");
        free(tbuf);
        r->ibcap = 0; r->ibuf = NULL;
    }
    r->lencodec = (int)(uint32_t)h[5]; r->idcodec = (int)(uint32_t)h[9];
    if (nosz(r->lencodec) || nosz(r->idcodec)) die("corrupt file: a NOSZ codec cannot be size-checked");
    if ((uint32_t)h[5] >> 31 || (uint32_t)h[9] >> 31) die("corrupt file: codec id out of range");
    r->ci = 0; r->cpend = 0; r->cvalid = 0;

    unsigned long long sb = ZSTD_getFrameContentSize(sz, h[3]);
    if (sb == ZSTD_CONTENTSIZE_UNKNOWN || sb == ZSTD_CONTENTSIZE_ERROR)
        die("block sizes stream has no usable frame size");
    uint8_t *sbuf = xmalloc((size_t)sb);
    if (ZSTD_isError(ZSTD_decompress(sbuf, sb, sz, h[3])))
        die("block sizes decompress failed");
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
    if (!r->blkbuf) die("oom: block buffer");
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
    if (r->names_int) {
        size_t nb; names = intnamesd(r->nz_keep, r->nzl_keep, r->nreads, &nb);
    } else if (r->names_tok) {
        size_t nb; names = tokd(r->nz_keep, r->nzl_keep, &nb);
    } else {
        unsigned long long nb = ZSTD_getFrameContentSize(r->nz_keep, r->nzl_keep);
        if (nb == ZSTD_CONTENTSIZE_UNKNOWN || nb == ZSTD_CONTENTSIZE_ERROR)
            die("names stream has no usable frame size");
        names = xmalloc((size_t)nb + 1);
        if (ZSTD_isError(ZSTD_decompress(names, nb, r->nz_keep, r->nzl_keep)))
            die("names decompress failed");
        names[nb] = 0;
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
                die("--sort --no-names keeps the rank of every record and no "
                    "name text: give the original FASTA/FASTQ in "
                    "--names FILE, or --ordinal-names to see the ranks");
            }
            nat = readnamefile(r->nreads);
        }
        else if (r->perm_tok) { size_t nb; nat = tokd(r->pz, r->pzl, &nb); }
        else {
            unsigned long long nb = ZSTD_getFrameContentSize(r->pz, r->pzl);
            if (nb == ZSTD_CONTENTSIZE_UNKNOWN || nb == ZSTD_CONTENTSIZE_ERROR)
                die("input-order names have no usable frame size");
            nat = xmalloc((size_t)nb + 1);
            if (ZSTD_isError(ZSTD_decompress(nat, nb, r->pz, r->pzl)))
                die("input-order names decompress failed");
            nat[nb] = 0;
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

static void buildindex(const char *path) {
    R *r = ropen(path);
    ensurenames(r);
    if (r->hasidx) { fprintf(stderr, "karyozip: %s already carries an index\n", path); rclose(r); return; }
    if (r->nreads > 0xffffffffu)
        die("more than 2^32 sequences: the index stores the ordinal in 32 bits");

    uint64_t n = r->nreads ? r->nreads : 1;
    uint64_t *key = xmalloc(mul(n, 8)), *key2 = xmalloc(mul(n, 8));
    uint32_t *ord = xmalloc(mul(n, 4));
    uint64_t *pair = xmalloc(mul(n, 16));          /* sort key and ordinal together */
    for (uint64_t i = 0; i < r->nreads; i++) {
        pair[2 * i] = fnvstr(r->name[i]);
        pair[2 * i + 1] = i;
    }
    qsort(pair, r->nreads, 16, kcmp);
    uint64_t dup = 0;
    for (uint64_t i = 0; i < r->nreads; i++) {
        key[i] = pair[2 * i];
        ord[i] = (uint32_t)pair[2 * i + 1];
        key2[i] = h64b(r->name[ord[i]]) ^ (ord[i] * 0x9e3779b97f4a7c15ULL);   /* bound to its ordinal */
        if (i && key[i] == key[i - 1]) dup++;
    }
    free(pair);
    /* Appended to the track, magic last: a torn write leaves no magic and
       the reader simply does not see an index. One file, not two. */
    uint64_t off0 = (uint64_t)r->fsz;
    FILE *f = fopen(path, "r+b");
    if (!f || fseeko(f, 0, SEEK_END)) die("cannot append the index to the track");
    uint64_t tail[4] = { r->nreads, dup, r->mhash ^ namesrc(), IDXMAGIC };
    int bad = fwrite(key, 8, r->nreads, f) != r->nreads
           || fwrite(key2, 8, r->nreads, f) != r->nreads
           || fwrite(ord, 4, r->nreads, f) != r->nreads
           || fwrite(tail, 8, 4, f) != 4 || fflush(f);
    if (bad) { (void)!ftruncate(fileno(f), (off_t)off0); fclose(f); die("index write failed"); }
    if (fclose(f)) die("index write failed");
    free(key); free(key2); free(ord);
    fprintf(stderr, "index: %llu sequences, %llu B appended, %llu hash collisions\n",
            (unsigned long long)r->nreads,
            (unsigned long long)(32 + r->nreads * 20), (unsigned long long)dup);
    rclose(r);
}

/* Returns the ordinal, or -1 when absent. */
static int64_t idxlookup(R *r, const char *path, const char *nm) {
    int fd = open(path, O_RDONLY);
    if (fd < 0) return -1;
    struct stat st;
    if (fstat(fd, &st) || (size_t)st.st_size < 32) { close(fd); return -1; }
    void *m = mmap(NULL, (size_t)st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    if (m == MAP_FAILED) return -1;
    /* the trailer is the last four words, magic last. Nothing in the file is
       aligned, so every word is read with memcpy, never through a cast. */
    const char *tail = (const char *)m + st.st_size - 32;
    uint64_t n = ld64(tail);
    if (ld64(tail + 24) != IDXMAGIC || n != r->nreads || n > (uint64_t)st.st_size / 20
        || (uint64_t)st.st_size < 32 + n * 20
        || ld64(tail + 16) != (r->mhash ^ namesrc())) { munmap(m, (size_t)st.st_size); return -1; }
    const char *base = (const char *)m + st.st_size - 32 - n * 20;
    #define KEY(i) ld64(base + (i) * 8)
    #define KEY2(i) ld64(base + n * 8 + (i) * 8)
    #define ORD(i) ld32(base + n * 16 + (i) * 4)
    uint64_t want = fnvstr(nm), want2 = h64b(nm);
    int64_t out = -1;
    uint64_t lo = 0, hi = n;
    while (lo < hi) {
        uint64_t mid = lo + (hi - lo) / 2;
        if (KEY(mid) == want) {
            uint64_t a2 = mid;
            while (a2 > 0 && KEY(a2 - 1) == want) a2--;
            uint64_t b2 = mid + 1;
            while (b2 < n && KEY(b2) == want) b2++;

            if (b2 - a2 == 1 && KEY2(a2) == (want2 ^ (ORD(a2) * 0x9e3779b97f4a7c15ULL)) && ORD(a2) < r->nreads
                && (!r->names_ready || !strcmp(r->name[ORD(a2)], nm))) { out = (int64_t)ORD(a2); break; }
            if (b2 - a2 == 1) break;
            ensurenames(r);
            for (uint64_t j = a2; j < b2; j++)
                if (ORD(j) < r->nreads && KEY2(j) == (want2 ^ (ORD(j) * 0x9e3779b97f4a7c15ULL)) && !strcmp(r->name[ORD(j)], nm)) {
                    out = (int64_t)ORD(j); break;
                }
            break;
        }
        if (KEY(mid) < want) lo = mid + 1; else hi = mid;
    }
    #undef KEY
    #undef KEY2
    #undef ORD
    munmap(m, (size_t)st.st_size);
    return out;
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

static void loadblk(R *r, uint64_t b) {

    if (b >= r->nblk) die("corrupt index: a block number is past the end");
    if (r->cached == b) return;
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
    if (r->lencodec) { g = paydecomp(cb, la, r->lencodec, r->blkbuf, r->b[b].lraw);
                       if (!g && r->b[b].lraw) die("frame decode failed"); }
    else
    { g = ZSTD_decompressDCtx(r->d, r->blkbuf, r->blkcap, cb, la);
      if (ZSTD_isError(g)) die(ZSTD_getErrorName(g)); }
    if (r->idcodec) { h2 = paydecomp(cb + la, lb, r->idcodec, r->blkbuf + g, r->b[b].raw);
                      if (!h2 && r->b[b].raw) die("frame decode failed"); }
    else
    { h2 = ZSTD_decompressDCtx(r->d, r->blkbuf + g, r->blkcap - g, cb + la, lb);
      if (ZSTD_isError(h2)) die(ZSTD_getErrorName(h2)); }
    free(cb);
    if (g != r->b[b].lraw || h2 != r->b[b].raw)
        die("corrupt block: a stream decodes to a size the index does not claim");
    r->lp0 = r->blkbuf; r->ip0 = r->blkbuf + g;
    PAYEND = r->blkbuf + g + h2;
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

static inline char *putu64(char *o, uint64_t v) {
    char t[20]; int n = 0;
    do { t[n++] = (char)('0' + v % 10); v /= 10; } while (v);
    while (n) *o++ = t[--n];
    return o;
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
    char **lab = xmalloc(cap * sizeof *lab), tmp[4096];
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
    size_t nlen = strlen(seq), cap2 = 1 << 20; char *ob = xmalloc(cap2), *o = ob;
    for (uint64_t k = a; k < b; k++) {
        if ((size_t)(o - ob) + nlen + strlen(lab[k]) + 48 > cap2) {
            if (fwrite(ob, 1, (size_t)(o - ob), stdout) != (size_t)(o - ob)) die("write failed");
            o = ob;
        }
        o = putrow(o, seq, nlen, st[k], en[k], lab[k], strlen(lab[k]));
    }
    if (o > ob && fwrite(ob, 1, (size_t)(o - ob), stdout) != (size_t)(o - ob)) die("write failed");
    for (uint64_t k = 0; k < n; k++) free(lab[k]);
    free(st); free(en); free(lab); free(ob); rclose(r);
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
        if (len && line[len - 1] == '\n') line[--len] = 0;
        if (len && line[len - 1] == '\r') line[--len] = 0;
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
    free(line); fclose(f);
    if (!NNODE) die("feature set not found in the hierarchy");
    for (int a = 0; a < NNODE; a++)   /* every node but categorized has a parent */
        if (NPAR[a] < 0 && strcmp(NODE[a], "categorized")) die("hierarchy: a node without a parent that is not categorized");
    for (int a = 0; a < NNODE; a++) { int n = 0;
        for (int c = a; c >= 0; c = NPAR[c]) if (++n > NNODE) die("the hierarchy has a cycle"); }
    ROOTID = nodeid("categorized", 11);
}

static void buildanc(void) {
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
        die("the reference .kst must be a plain track: not --aln, --db, --sort or --names drop, names in the file");
    R *r = refopen(path);
    if (!r->rowhash) die("the reference .kst carries no row hash: encode it again with this version");
    ensurenames(r);
    REFLAB = xmalloc(mul(r->nlab ? r->nlab : 1, 4));
    for (uint64_t j = 0; j < r->nlab; j++) REFLAB[j] = (uint32_t)intern(r->lab[j], strlen(r->lab[j]));
    uint64_t tot = 0;
    for (uint64_t i = 0; i < r->nreads; i++) {
        if (rcon_get(r->name[i], 0)) die("reference track: a contig's rows are not together");
        rcon *rc = rcon_get(r->name[i], 1);
        rc->ri = (int64_t)i; rc->filled = 0;
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
        rc->st = xmalloc(mul(n, 8)); rc->lab = xmalloc(mul(n, 4)); rc->n = rc->cap = q.nruns;
        uint64_t s = 0;
        for (uint64_t k = 0; k < q.nruns; k++) {
            spannext(r, &sp, &lp, &ip); sp.left--;
            uint64_t e = s + rlen(&lp, rc->name, s, k + 1 == q.nruns && (uint64_t)rc->ri >= TAILSKIP_R), f = vget(&ip);
            if (f >= r->nlab) die("corrupt block: a label id past the table");
            if (e <= s) die("reference track: a run of no length");
            rc->st[k] = s; rc->lab[k] = REFLAB[f]; s = e;
        }
        rc->end = s; r->curblk = UINT64_MAX;
        LENDEC = ld; TAILDEC = td; TAILPRED_R = tp; TAILSKIP_R = ts; PAYEND = pe;
        rc->filled = 1;
    }
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
            const char *nm = (cf == ROOTID && cn) ? "novel" : NODE[cf];
            size_t ml = strlen(nm);
            if (v1 && cf != ROOTID && ml > 9 && !strcmp(nm + ml - 9, "_specific")) ml -= 9;
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
    free(r->blkbuf); free(r->b);
    free(r->nz_keep); free(r->pz); free(r->ibuf); free(r->ichunk);
    free(r->coff); free(r->cpendv); free(r->htab);
    free(r->sst); free(r->sen); free(r->sfeat); free(r->snov);
    free(r);
}

static void smooth_mt(const char *path, R *r0, int64_t first, int64_t last,
                      int *labmap, uint64_t lo, uint64_t hi, int v1, int nthr) {

    const int64_t CHUNK = 4096;
    int64_t nchunk = (last - first + CHUNK - 1) / CHUNK;
    _Atomic int failed = 0, wfail = 0;
#ifdef _OPENMP
#pragma omp parallel num_threads(nthr)
#endif
    {
        R *rr = ropen(path);

        if (rr) { rr->name = r0->name; rr->names_ready = 1; rr->name_borrowed = 1; }
        size_t cap = 1 << 20, len = 0;
        char *ob = malloc(cap);
        if (!ob || !rr) failed = 1;
#ifdef _OPENMP
#pragma omp for ordered schedule(dynamic, 1)
#endif
        for (int64_t c = 0; c < nchunk; c++) {
            if (failed) continue;
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
        if (rr) rclose(rr);
    }
    if (wfail) die("write failed");
    if (failed) die("threaded smoothing failed");
}

/* smooth TRACK HIER FS [SEQ [START END]] */
static void smoothcmd(const char *path, const char *hier, const char *fs,
                      const char *seq, uint64_t lo, uint64_t hi, int v1) {
    loadhier(hier, fs);
    buildanc();
    R *r = ropen_informed(path);
    ensurenames(r);   /* emit prints names */
    int *labmap = labnodes(r, fs);
    if (ANC) { free(ANC); free(LCA); buildanc(); }   /* unknown labels added above */

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
            char *e3 = NULL;
            (void)e3;
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
    /* --out DIR: one BED per selected set, staged then swapped in together.
       No --out: the selected columns, tab-joined, merged where they all
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
            if (pn < 0 || (size_t)pn >= sizeof pb - 16) die("--out: the output path is too long");
            if (!lstat(pb, &ts) && !S_ISREG(ts.st_mode)) die("--out: a set's BED exists and is not a regular file");
        }
        DBTMP = xcalloc(NSETK, sizeof *DBTMP); ndbtmp = NSETK; atexit(dbtmpclean);
        signal(SIGINT, dbtmpsig); signal(SIGTERM, dbtmpsig); signal(SIGHUP, dbtmpsig);
        for (size_t i = 0; i < NSETK; i++) {
            int pn = snprintf(pb, sizeof pb, "%s/%s.bed.tmp.%d", dir, fset[SETK[i]], (int)getpid());
            if (pn < 0 || (size_t)pn >= sizeof pb - 16) die("--out: the output path is too long");
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
int main(int argc, char **argv) {
    signal(SIGXFSZ, SIG_IGN);
    for (int fd = 0; fd < 3; fd++) if (fcntl(fd, F_GETFD) < 0) die("stdin, stdout or stderr is closed");
    atexit(profdump);
    atexit(emitflush);
    static const char *const USAGE = "usage:\n"
            "  karyozip encode [--block N|all] [--level N] [--slow: smaller, slower]\n"
            "                  [--names keep|drop|drop:FILE]\n"
            "                  [--sort | --sort=MEM]\n"
            "                  [--db features.tsv [--hier H]] [--aln P.paf --ref R.bed|R.kst [--aln-names F]]\n"
            "  karyozip decode F | --db features.tsv [--hier H] [--out DIR] [--set A,B] F | --aln P.paf --ref R F\n"
            "                  --names F reads the names, from drop:FILE or a FASTA, FASTQ or list\n"
            "                  --names ordinal accepts a track whose names were dropped\n"
            "                  --raw shows what is stored, reversing no transform\n"
            "  karyozip aln-pack P.paf REF.bed OUT   only what --aln reads\n"
            "  karyozip region F SEQ                 the whole sequence\n"
            "  karyozip region F SEQ POS             the run covering POS\n"
            "  karyozip region F SEQ START END [HALO]\n"
            "  karyozip smooth F HIER FS [SEQ [START END]]\n"
            "  karyozip index  F                     appends an index to the track itself\n"
            "  any command: --threads N, -t N  --profile  --tmp-dir DIR  --help, -h\n"
            "               --version prints the zstd, libbz3 and htscodecs versions\n"
            "  encode:      --block N|all  --level N  --slow  --sort[=MEM]\n"
            "               --names keep|drop|drop:FILE  --db F  --aln P --ref R\n"
            "               --tail 0|1  --tok-chunk N  --id-codec N  --len-codec N\n"
            "  decode:      --raw  --out DIR\n"
            "  readers:     --names FILE|ordinal  --db F  --aln P --ref R, as encode was given;\n"
            "               --db: --set A,B picks feature sets, default all; decode without --out and\n"
            "               region print the picked labels as columns, merged where they all agree\n"
            "  with --aln:  --aln-names FILE, on encode and decode, when the track\n"
            "               names rows by ordinal\n"
            "  smooth:      --max-gap N  --v1\n"
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
    static const struct { const char *f; int nv; } FLAGS[] = {
        {"--version",0},{"--threads",1},{"-t",1},{"--profile",0},{"--tmp-dir",1},{"--names",1},
        {"--tok-chunk",1},{"--tail",1},{"--id-codec",1},{"--len-codec",1},
        {"--max-gap",1},{"--v1",0},{"--aln-names",1},{"--sort",0},
        {"--slow",0},
        {"--block",1},{"--level",1},{"--db",1},{"--aln",1},
        {"--ref",1},{"--out",1},{"--set",1},{"--raw",0},{"--hier",1},
    };
    const char *cmd = "";   /* the first word that is not a flag or its value */
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (a[0] != '-' || !a[1]) { cmd = a; break; }
        int nv = -1;
        if (!strncmp(a, "--sort=", 7)) nv = 0;        /* the value is in the flag */
        else for (size_t k = 0; k < sizeof FLAGS / sizeof *FLAGS; k++)
            if (!strcmp(a, FLAGS[k].f)) { nv = FLAGS[k].nv; break; }
        if (nv < 0) { fprintf(stderr, "karyozip: unknown flag %s\n", a); exit(1); }
        i += nv;
    }
    int enc = !strcmp(cmd, "encode");
    int sm = !strcmp(cmd, "smooth");
    int dec = !strcmp(cmd, "decode");
    int rd = sm || dec || !strcmp(cmd, "region") || !strcmp(cmd, "index");
    int inf = enc || rd;                   /* --db and --aln: every reader but index */
    if (!strcmp(cmd, "index")) inf = 0;
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
            if (!strcmp(v, "keep")) { ok = enc; if (NAMEMODE == 1) die("--names keep and --names drop are exclusive"); NAMEMODE = 2; }
            else if (!strcmp(v, "drop")) { ok = enc; if (NAMEMODE == 2) die("--names keep and --names drop are exclusive"); NAMEMODE = 1; }
            else if (!strncmp(v, "drop:", 5)) { ok = enc; NAMES_OUT = v + 5;
                                                if (!*NAMES_OUT) die("--names drop:FILE needs a file"); }
            else if (!strcmp(v, "ordinal")) {
                if (enc) die("encode takes --names keep, drop or drop:FILE");
                ok = rd; ORDNAMES = 1; }
            else if (enc) die("encode takes --names keep, drop or drop:FILE; a name file is --aln-names");
            else { ok = rd; NAMESARG = v; }   /* the blob drop:FILE wrote, or a FASTA, FASTQ or list */
            take = 2; }
        else if (!strcmp(a, "--tok-chunk")) { ok = enc;
            TOKCHUNK = (uint32_t)int_arg(v, 1, 0x7fffffff, "--tok-chunk needs a number");
            if (TOKCHUNK > TOKMAX) TOKCHUNK = TOKMAX;
            TOKSET = 1; take = 2; }
        else if (!strcmp(a, "--tail")) { ok = enc; TAILFORCE = (int)int_arg(v, 0, 1, "--tail is 0 or 1"); take = 2; }
        else if (!strcmp(a, "--id-codec")) { ok = enc;
            IDCODEC = (int)codec_arg(v, "--id-codec needs a codec id"); take = 2; }
        else if (!strcmp(a, "--len-codec")) { ok = enc;
            LENCODEC = (int)codec_arg(v, "--len-codec needs a codec id"); take = 2; }
        else if (!strcmp(a, "--max-gap")) { ok = sm; MAXGAP = num_arg(v, UINT64_MAX, "--max-gap needs a number"); take = 2; }
        else if (!strcmp(a, "--v1")) { ok = sm; V1 = 1; take = 1; }
        else if (!strcmp(a, "--aln-names")) {   /* the alignment's own table: ordinal N is line N+1 */
            ok = inf; ALNNAMES = v; take = 2; }
        else if (!strcmp(a, "--sort")) { ok = enc; SORT = 1; take = 1; }
        else if (!strncmp(a, "--sort=", 7)) { ok = enc; SORT = 1; SORTMEM = mem_arg(a + 7); take = 1; }
        else if (!strcmp(a, "--slow")) { ok = enc; SLOW = 1; FASTMODE = 0; take = 1; }
        else if (!strcmp(a, "--block")) { ok = enc;
            BLOCK = !strcmp(v, "all") ? UINT64_MAX      /* one block, no random access */
                  : (uint64_t)int_arg(v, 1, 0x7fffffffffffffffL,
                                      "--block is a number of runs, 1 or more, or \"all\"");
            BLOCKSET = 1; take = 2; }
        else if (!strcmp(a, "--level")) { ok = enc;
            LEVEL = (int)int_arg(v, 1, 22, "--level is a zstd level, 1 to 22"); LEVELSET = 1; take = 2; }
        else if (!strcmp(a, "--db")) { ok = inf; DBFEATS = v; take = 2; }
        else if (!strcmp(a, "--hier")) { ok = inf; HIERTAB = v; take = 2; }
        else if (!strcmp(a, "--aln")) { ok = inf; ALNPAF = v; take = 2; }
        else if (!strcmp(a, "--ref")) { ok = inf; ALNREF = v; take = 2; }
        else if (!strcmp(a, "--out")) { ok = dec; OUTDIR = v; take = 2; }
        else if (!strcmp(a, "--set")) { ok = dec || !strcmp(cmd, "region"); SETNAME = v; take = 2; }
        else if (!strcmp(a, "--raw")) { ok = dec; RAW = 1; take = 1; }
        if (take && !ok) { fprintf(stderr, "karyozip: %s does not apply to %s\n", a, cmd); exit(1); }
        if (take) { for (int j = i; j + take <= argc; j++) argv[j] = argv[j + take]; argc -= take; i--; }
    }
    if (!*cmd) die(USAGE);   /* flags but no command: argv[1] is now past the end */
    for (int i = 2; i < argc; i++)   /* what the loop left is positional; "-" alone is stdin */
        if (argv[i][0] == '-' && argv[i][1]) { fprintf(stderr, "karyozip: unknown flag %s\n", argv[i]); exit(1); }
    if ((IDCODEC >= 0 && badcodec(IDCODEC)) || (LENCODEC >= 0 && badcodec(LENCODEC))) die("--id-codec and --len-codec take a codec id, not any number");
    if (nosz(IDCODEC) || nosz(LENCODEC)) die("a NOSZ codec (0x10) stores no size and is not accepted");
    if (noext(IDCODEC) || noext(LENCODEC)) die("arith X_EXT (0x04) needs libbz2 and is not accepted");
    if (NTHREADS > 1024) NTHREADS = 1024;   /* libgomp keeps team state on the stack */
#ifdef _OPENMP
    if (omp_get_max_threads() > 1024) omp_set_num_threads(1024);   /* OMP_NUM_THREADS too */
#endif


    WALL0 = nowsec();

    if (!strcmp(argv[1], "aln-pack")) {

        if (argc != 5 || argv[2][0] == '-' || argv[3][0] == '-' || argv[4][0] == '-')
            die("usage: karyozip aln-pack P.paf REF.bed OUT");
        loadref(argv[3]);
        alnpack(argv[2], argv[4]);
        return 0;
    }
    if (!strcmp(argv[1], "encode")) {

        char *freearg[64]; int nfree = 0;
        for (int i = 2; i < argc; i++) {
            if (nfree == 64) die("too many input files");
            freearg[nfree++] = argv[i];
        }
        if (!DBFEATS) {                 /* --db takes its BEDs positionally */
            if (nfree) die("without --db, encode reads one BED on stdin");
        }
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
            if (!nfree) die("--db: give one BED per feature set");
            encode_db(nfree, freearg);
        } else encode();
    }

    else if (!strcmp(argv[1], "decode") && argc > 2) {
        const char *f = NULL;
        for (int i = 2; i < argc; i++) {
            if (f) die("decode takes one track");
            f = argv[i];
        }
        if (!f) die("decode needs a track");
        loadinformed("decode ");
        if (DBFEATS) decode_db(f, OUTDIR); else dump(f);
    }

    /* one reader: SEQ alone is the whole sequence, SEQ POS the run covering
       POS, SEQ START END the overlap, plus an optional halo. */
    else if (!strcmp(argv[1], "region") && argc > 3) {
        uint64_t rlo = argc > 4 ? num_arg(argv[4], UINT64_MAX - 1, "region: start is not a number") : 0;
        uint64_t rhi = argc > 5 ? num_arg(argv[5], UINT64_MAX, "region: end is not a number") : 0;
        if (argc == 5) rhi = rlo + 1;                  /* one position */
        if (argc > 5 && rlo >= rhi) die("region: give START END with END past START");
        loadinformed("region ");
        region(argv[2], argv[3], rlo, rhi, argc > 6 ? num_arg(argv[6], UINT64_MAX, "region: halo is not a number") : 0);
    }
    else if (!strcmp(argv[1], "index") && argc > 2) { INDEXING = 1; buildindex(argv[2]); }
    else if (!strcmp(argv[1], "smooth") && argc > 4) {
        if (argc == 7 || argc > 8) die("smooth: give SEQ, or SEQ START END");
        uint64_t slo = argc > 6 ? num_arg(argv[6], UINT64_MAX, "smooth: start is not a number") : 0;
        uint64_t shi = argc > 7 ? num_arg(argv[7], UINT64_MAX, "smooth: end is not a number") : 0;
        if (argc > 7 && slo >= shi) die("smooth: end must be past start");
        if (DBFEATS && !HIERTAB) HIERTAB = argv[3];   /* the hierarchy it smooths with */
        loadinformed("smooth ");
        smoothcmd(argv[2], argv[3], argv[4], argc > 5 ? argv[5] : NULL, slo, shi, V1);
    }
    else die("bad command or missing arguments");
    return 0;
}
#endif   /* KARYOZIP_LIB: the library has no main */
