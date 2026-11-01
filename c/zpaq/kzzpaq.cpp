// karyozip's zpaq codec: libzpaq level 5 or the ilv model on the caller's pieces, side by side under OpenMP; frame: count, (size, block size) u32 LE each, the blocks
#include "libzpaq.h"
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>

namespace libzpaq { void error(const char* msg) { throw std::runtime_error(msg); } }

// -e high's model of an interleaved stream (k, then k label-length pairs, LEB128): HCOMP follows the fields, the read's earlier values and the reads before
// R0 next field, R1 pairs left, R2 k, R3 pair index, R4 previous label, R5 label, R6 length sum, R7 previous length, R8 previous k, R9-R10 partial value and shift, R11 and R20 labels before R4, R12 byte, R13 value, R14 state, R15 field + 4 on the last pair, R16 last label of the previous read, R17 k two reads back, R18-R19 first labels of this read and the previous, R21 hash of the read so far, R22-R25 where this read and the three before start in M
static const char ILV[] = R"zp(
comp 6 10 0 0 37
  0 icm 10
  1 isse 13 0
  2 isse 16 1
  3 isse 18 2
  4 isse 19 3
  5 isse 19 4
  6 icm 16
  7 icm 16
  8 icm 16
  9 cm 22 255
  10 cm 22 255
  11 match 22 24
  12 icm 16
  13 icm 16
  14 icm 16
  15 icm 16
  16 icm 18
  17 icm 16
  18 icm 16
  19 icm 16
  20 isse 18 5
  21 isse 18 20
  22 match 22 24
  23 match 22 24
  24 match 22 24
  25 match 22 24
  26 match 22 24
  27 match 22 24
  28 match 22 24
  29 match 22 24
  30 mix 12 0 30 12 255
  31 mix 18 0 30 12 255
  32 mix2 0 30 31 12 0
  33 sse 19 32 32 255
  34 sse 20 32 32 255
  35 mix2 0 33 34 12 0
  36 mix2 0 32 35 12 0
hcomp
  c++ *c=a r=a 12
  a> 127 ifl
    a=r 12 a&= 127 b=r 10 a<<=b b=r 9 a+=b r=a 9
    a=r 10 a+= 7 r=a 10
  elsel
    a=r 12 b=r 10 a<<=b b=r 9 a+=b r=a 13
    a=0 r=a 9 r=a 10
    a=r 0 a== 0 ifl
      a=r 8 r=a 17 a=r 2 r=a 8
      a=r 13 r=a 2 r=a 1
      a=0 r=a 3 r=a 6 r=a 21
      a=r 24 r=a 25 a=r 23 r=a 24 a=r 22 r=a 23 a=c r=a 22
      a= 255 r=a 4 r=a 7 r=a 11 r=a 20
      a=r 1 a> 0 if a= 1 else a=0 endif r=a 0
    elsel
      a== 1 if
        a=r 13 r=a 5 a=r 3 a== 0 if a=r 18 r=a 19 a=r 13 r=a 18 endif a= 2 r=a 0
      else
        a=r 13 r=a 7 b=r 6 a+=b r=a 6
        a=r 11 r=a 20 a=r 4 r=a 11 a=r 5 r=a 4
        a=r 3 a++ r=a 3
        a=r 1 a-- r=a 1
        a> 0 if a= 1 else a=r 5 r=a 16 a=0 endif r=a 0
      endif
    endif
  endif
  a=r 21 b=c hash r=a 21
  a=r 1 a== 1 if a= 4 else a=0 endif b=r 0 a+=b r=a 15
  b=r 10 a+=b a+=b r=a 14
  d= 0 *d=0 a=r 14 hashd a=r 3 a> 15 if a= 15 endif hashd a=r 2 a> 31 if a= 31 endif hashd
  d= 1 b=c a=r 14 hash *d=a
  d= 2 b-- hash *d=a
  d= 3 b-- hash *d=a
  d= 4 b-- hash *d=a
  d= 5 b-- hash b-- hash *d=a
  d= 6 *d=0 a=r 0 hashd
  a=r 0 a== 1 ifl
    a=r 4 hashd a=r 11 hashd a=r 15 hashd a=r 3 a> 7 if a= 7 endif hashd
  elsel
    a== 2 if
      a=r 5 hashd a=r 6 hashd a=r 15 hashd
    else
      a=r 8 hashd
    endif
  endif
  a=r 9 hashd a=r 10 hashd a=0 hashd
  d= 7 *d=0 a=r 0 hashd
  a=r 0 a== 1 ifl
    a=r 4 hashd a=r 7 hashd a=r 6 hashd
  elsel
    a== 2 if
      a=r 5 hashd a=r 4 hashd a=r 7 hashd a=r 15 hashd
    else
      a=r 8 hashd a=r 16 hashd
    endif
  endif
  a=r 9 hashd a=r 10 hashd a=0 hashd
  d= 8 *d=0 a=r 0 hashd
  a=r 0 a== 1 ifl
    a=r 4 hashd a=r 3 hashd a=r 1 hashd
  elsel
    a== 2 if
      a=r 5 hashd a=r 6 hashd a=r 7 hashd a=r 4 hashd
    else
      a=r 8 hashd a=r 16 hashd a=r 17 hashd
    endif
  endif
  a=r 9 hashd a=r 10 hashd a=0 hashd
  d= 6 a=*d d= 9 *d=a
  d= 7 a=*d d= 10 *d=a
  d= 11 a=r 14 b=c hash b-- hash b-- hash b-- hash b-- hash b-- hash b-- hash b-- hash *d=a
  d= 12 *d=0 a=r 0 hashd
  a=r 0 a== 1 ifl
    a=r 3 a== 0 if
      a=r 19 hashd a=r 16 hashd a=r 2 hashd
    else
      a=r 4 hashd a=r 11 hashd a=r 18 hashd
    endif
  elsel
    a== 2 if
      a=r 5 hashd a=r 7 hashd a=r 18 hashd a=r 6 a>>= 3 hashd
    else
      a=r 18 hashd a=r 16 hashd
    endif
  endif
  a=r 9 hashd a=r 10 hashd a=0 hashd
  d= 13 *d=0 a=r 14 hashd b=c b-- a=*b hashd b-- a=*b hashd
  d= 14 *d=0 a=r 14 hashd b=c b-- a=*b hashd b-- b-- a=*b hashd
  d= 15 *d=0 a=r 0 hashd
  a=r 0 a== 1 ifl
    a=r 4 hashd a=r 11 hashd a=r 20 hashd a=r 15 hashd
  elsel
    a== 2 if
      a=r 5 hashd a=r 4 hashd a=r 11 hashd a=r 7 hashd
    else
      a=r 16 hashd a=r 19 hashd a=r 8 hashd a=r 17 hashd
    endif
  endif
  a=r 9 hashd a=r 10 hashd a=0 hashd
  d= 16 *d=0 a=r 21 hashd a=r 14 hashd
  d= 17 *d=0 a=r 14 hashd a=c a++ b=r 22 a-=b b=r 23 a+=b b=a a=*b hashd a=r 3 hashd
  d= 18 *d=0 a=r 14 hashd a=c a++ b=r 22 a-=b b=r 23 a+=b b=a a=*b hashd a=*c hashd
  d= 19 *d=0 a=r 14 hashd a=c a++ b=r 22 a-=b b=r 23 a+=b b=a a=*b hashd a=c a++ b=r 22 a-=b b=r 24 a+=b b=a a=*b hashd
  d= 6 a=*d d= 20 *d=a
  d= 7 a=*d d= 21 *d=a
  d= 22 a=r 21 *d=a
  d= 23 *d=0 a=r 14 hashd a=c a++ b=r 22 a-=b b=r 23 a+=b b=a a=*b hashd a=r 3 hashd
  d= 24 *d=0 a=r 14 hashd a=c a++ b=r 22 a-=b b=r 23 a+=b b=a a=*b hashd a=*c hashd
  d= 25 a=r 14 b=c hash b-- hash b-- hash *d=a
  d= 26 *d=0 a=r 14 hashd a=r 4 hashd a=r 11 hashd a=r 20 hashd a=r 3 hashd
  d= 27 a=r 14 b=c hash *d=a
  d= 28 a=r 14 b=c hash b-- hash *d=a
  d= 7 a=*d d= 29 *d=a
  d= 30 a=r 3 a> 0 if a= 8 endif b=r 15 a+=b a<<= 8 *d=a
  d= 31 a=r 0 a<<= 8 a+=*c a<<= 8 *d=a
  d= 33 a=r 15 a<<= 8 a+=*c a<<= 8 *d=a
  d= 6 a=*d a>>= 20 a<<= 8 d= 34 *d=a
  halt
end
)zp";

namespace {
struct Out : libzpaq::Writer {   // a piece's zpaq block
  std::string s;
  void put(int c) { s.push_back((char)c); }
  void write(const char* b, int n) { s.append(b, n); }
};
struct In : libzpaq::Reader {    // a zpaq block, read back
  const unsigned char *p, *e;
  int get() { return p < e ? *p++ : -1; }
  int read(char* b, int n) { int k = e - p < n ? (int)(e - p) : n; memcpy(b, p, k); p += k; return k; }
};
struct Fill : libzpaq::Writer {  // its bytes into their place, never past it
  unsigned char *p, *e; bool over = false;
  void put(int c) { if (p < e) *p++ = (unsigned char)c; else over = true; }
  void write(const char* b, int n) { if (e - p < n) { over = true; n = (int)(e - p); } memcpy(p, b, n); p += n; }
};
void w32(unsigned char* q, uint32_t v) { for (int i = 0; i < 4; i++) q[i] = (unsigned char)(v >> 8 * i); }
uint32_t r32(const unsigned char* q) { return q[0] | q[1] << 8 | q[2] << 16 | (uint32_t)q[3] << 24; }
}

extern "C" {

// n bytes in ncut pieces of cut[i] bytes, by libzpaq's method ("5") or "ilv", into out of cap bytes: the frame's size, 0 on any failure
size_t kz_zpaq_encode(const unsigned char* in, size_t n, const uint32_t* cut, uint32_t ncut,
                      const char* method, unsigned char* out, size_t cap) {
  std::vector<size_t> off(ncut + 1, 0);
  for (uint32_t i = 0; i < ncut; i++) off[i + 1] = off[i] + cut[i];
  if (off[ncut] != n) return 0;
  std::vector<std::string> blk(ncut);
  int bad = 0;
#pragma omp parallel for schedule(dynamic, 1) if (ncut > 1 && !omp_in_parallel())
  for (long i = 0; i < (long)ncut; i++) {
    try {
      if (cut[i] > 0x7fffffffu) throw 1;
      Out o;
      if (!strcmp(method, "ilv")) {
        In r; r.p = in + off[i]; r.e = r.p + cut[i];
        libzpaq::Compressor c;
        c.setOutput(&o); c.setInput(&r);
        c.startBlock(ILV, 0);
        c.startSegment(0, 0);
        c.compress();
        c.endSegment(0);
        c.endBlock();
      } else {
        libzpaq::StringBuffer sb(cut[i]);
        sb.write((const char*)in + off[i], (int)cut[i]);
        libzpaq::compressBlock(&sb, &o, method, 0, 0, false);
      }
      blk[i].swap(o.s);
    } catch (...) {
#pragma omp atomic write
      bad = 1;
    }
  }
  size_t z = 4 + 8 * (size_t)ncut;
  for (const std::string& b : blk) { if (b.size() > 0xffffffffu) bad = 1; z += b.size(); }
  if (bad || z > cap) return 0;
  w32(out, ncut);
  unsigned char* q = out + 4 + 8 * (size_t)ncut;
  for (uint32_t i = 0; i < ncut; i++) {
    w32(out + 4 + 8 * (size_t)i, cut[i]); w32(out + 8 + 8 * (size_t)i, (uint32_t)blk[i].size());
    memcpy(q, blk[i].data(), blk[i].size()); q += blk[i].size();
  }
  return z;
}

// the frame back into raw bytes at dst: 0, or -1 when it does not come back to exactly raw bytes

int kz_zpaq_decode(const unsigned char* in, size_t n, unsigned char* dst, size_t raw) {
  if (n < 4) return -1;
  uint32_t np = r32(in);
  if ((n - 4) / 8 < np) return -1;
  std::vector<size_t> io(np + 1), oo(np + 1);
  io[0] = 4 + 8 * (size_t)np; oo[0] = 0;
  for (uint32_t i = 0; i < np; i++) { oo[i + 1] = oo[i] + r32(in + 4 + 8 * (size_t)i); io[i + 1] = io[i] + r32(in + 8 + 8 * (size_t)i); }
  if (io[np] != n || oo[np] != raw) return -1;
  int bad = 0;
#pragma omp parallel for schedule(dynamic, 1) if (np > 1 && !omp_in_parallel())
  for (long i = 0; i < (long)np; i++) {
    try {
      In r; r.p = in + io[i]; r.e = in + io[i + 1];
      Fill f; f.p = dst + oo[i]; f.e = dst + oo[i + 1];
      libzpaq::decompress(&r, &f);
      if (f.over || f.p != f.e) throw 1;
    } catch (...) {
#pragma omp atomic write
      bad = 1;
    }
  }
  return bad ? -1 : 0;
}

}
