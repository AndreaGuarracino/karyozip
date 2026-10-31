// karyozip's zpaq codec: libzpaq's level 5 on pieces the caller cuts,
// compressed and decompressed side by side under OpenMP. The frame: the
// piece count, then each piece's size and its zpaq block's size, all u32
// little endian, then the blocks in order: the same bytes on any number
// of threads.
#include "libzpaq.h"
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>

namespace libzpaq { void error(const char* msg) { throw std::runtime_error(msg); } }

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

// n bytes in ncut pieces of cut[i] bytes each, by libzpaq's method ("5" its
// slowest), into out of cap bytes: the frame's size, 0 on any failure
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
      libzpaq::StringBuffer sb(cut[i]);
      sb.write((const char*)in + off[i], (int)cut[i]);
      Out o;
      libzpaq::compressBlock(&sb, &o, method, 0, 0, false);
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

// the frame back into raw bytes at dst: 0, or -1 when it does not come back
// to exactly raw bytes
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
