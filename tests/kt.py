"""Forging helpers for the suite: parse a track, patch it, re-sign its hashes."""
import struct, subprocess, sys, random, os
FNV0 = 0xcbf29ce484222325; P = 0x100000001b3; M = (1 << 64) - 1
# the h[1] flag bits, as FLAGS in rust/src/kst.rs (run.sh checks the two agree); the rest: names length
FLAGS = sum(1 << b for b in (49, 50, 51, 52, 55, 56, 57, 58, 59, 60, 61, 62, 63))
def fnv(h, b):
    for x in b: h = ((h ^ x) * P) & M
    return h
def load(path):
    d = bytearray(open(path, "rb").read()); return d, list(struct.unpack_from("<16Q", d, 8))
def sign(d, h):   # hash header (h[15]=0), then lz..itz sections
    hz = list(h); hz[15] = 0; mh = fnv(FNV0, struct.pack("<16Q", *hz)); off = 136
    for n in (h[6], h[1] & ~FLAGS, h[2], h[3], h[11]): mh = fnv(mh, d[off:off + n]); off += n
    h[15] = mh; struct.pack_into("<16Q", d, 8, *h)
def vget(d, i):
    v = 0; s = 0
    while True:
        c = d[i]; i += 1; v |= (c & 0x7f) << s; s += 7
        if c < 0x80: return v, i
def vput(v):
    o = bytearray()
    while v >= 0x80: o.append(v & 0x7f | 0x80); v >>= 7
    o.append(v); return bytes(o)
def zst(W, body):   # via a file: frame keeps content size
    open(W + "/c", "wb").write(body)
    subprocess.run(["zstd", "-q", "-f", "-19", "--no-check", W + "/c", "-o", W + "/c.zst"], check=True)
    return open(W + "/c.zst", "rb").read()
def unzst(z):
    r = subprocess.run(["zstd", "-d", "-q", "-c"], input=bytes(z), capture_output=True)
    if r.returncode:
        raise RuntimeError("zstd -d failed: " + r.stderr.decode("utf-8", "replace").strip())
    return r.stdout
def gen(W, N):   # awkward valid BEDs for fuzz.sh, diff.sh
    random.seed(13)
    byt = [c for c in range(1, 256) if c not in (9, 10)]
    def nm(i):
        k = random.randrange(6)
        if k == 0: return str(random.choice([2**50 - 1, 2**50, 2**50 + 1, 2**53, 10**18 - 1])).encode()
        if k == 1: return ("0" * random.randint(1, 3) + str(random.randrange(10**6))).encode()
        if k == 2: return bytes(random.choice(byt) for _ in range(random.choice([1, 2, 63, 64, 511])))
        if k == 3: return ("r%d" % i).encode()
        if k == 4: return bytes(random.choice(byt) for _ in range(random.randint(1, 40)))
        return str(i).encode()
    for f in range(N):
        rows = []; seen = set()
        for i in range(random.randint(1, 60)):
            n = nm(i)
            while n in seen: n = nm(i + 1000)
            seen.add(n); p = 0
            for _ in range(random.randint(1, 50)):
                L = random.choice([1, 2, 7, 255, 256, 65535, 65536, random.randint(1, 3000000)])
                lab = bytes(random.choice([c for c in byt if c != 9] + [9]) for _ in range(random.randint(1, 12))) if random.random() < 0.3 else b"L%d" % random.randrange(5)
                rows.append(n + b"\t%d\t%d\t" % (p, p + L) + lab + b"\n"); p += L
        open(os.path.join(W, "ok%d.bed" % f), "wb").write(b"".join(rows))

def bgzf(src, dst):   # bgzip-style BGZF: 65280-byte blocks, empty EOF
    import zlib
    d = open(src, "rb").read(); o = bytearray()
    for i in list(range(0, len(d), 65280)) + [len(d)]:
        c = d[i:i + 65280]; z = zlib.compressobj(6, zlib.DEFLATED, -15); x = z.compress(c) + z.flush()
        o += b"\x1f\x8b\x08\x04" + bytes(6) + struct.pack("<HBBHH", 6, 66, 67, 2, len(x) + 25) + x
        o += struct.pack("<II", zlib.crc32(c), len(c))
    open(dst, "wb").write(o)

def blk(c, isz=None, crc=None, pre=b"", bsz=None):   # one BGZF block of c; pre: subfields before BC; the rest forged
    import zlib
    z = zlib.compressobj(6, zlib.DEFLATED, -15); x = z.compress(c) + z.flush(); xl = len(pre) + 6
    bs = 12 + xl + len(x) + 8 if bsz is None else bsz
    return (b"\x1f\x8b\x08\x04" + bytes(6) + struct.pack("<H", xl) + pre + struct.pack("<BBHH", 66, 67, 2, bs - 1) + x
            + struct.pack("<II", zlib.crc32(c) if crc is None else crc, len(c) if isz is None else isz))

def bgzfx(dst):   # forged BGZF: ISIZE bounds, subfields before BC, cuts inside the extra field
    def bed(n):   # a BED of exactly n bytes
        r = b"".join(b"r%d\t0\t1\tA\n" % i for i in range(n // 20)); pad = n - len(r) - 7
        return r + b"z\t0\t1\t" + b"B" * pad + b"\n"
    eof = blk(b""); w = lambda n, b: open(os.path.join(dst, n), "wb").write(b)
    for n in (65535, 65536, 65537):
        w("isz%d.bed" % n, bed(n)); w("isz%d.gz" % n, blk(bed(n)) + eof)
    w("huge.gz", b"".join(blk(b"r\t0\t1\tA\n", isz=0xffffffff) for _ in range(8)) + eof)     # each claims 4 GiB
    w("many.gz", blk(b"r\t0\t1\tA\n") * 1 + blk(b"", isz=65536, crc=0) * 300000 + eof)   # 20 GB claimed, each in bounds
    w("bsize.gz", blk(b"r\t0\t1\tA\n", pre=b"XY\x04\x00abcd", bsz=30) + eof)              # BSIZE ends inside its own header
    s = bed(1000); w("bc.bed", s)
    one = blk(s, pre=b"XY\x04\x00abcd"); w("bcafter.gz", one + blk(s[:0], pre=b"XY\x00\x00") + eof)
    for k in range(10, 12 + 14 + 3):
        w("cut1.%d.gz" % k, one[:k]); w("cut2.%d.gz" % k, blk(s) + one[:k])
    tail = b"\x1f\x8b\x08\x04" + bytes(6) + struct.pack("<H", 12) + b"XY\x02\x00ab" + b"BC\x02\x00"   # BC payload past the end
    w("bcend1.gz", tail); w("bcend2.gz", blk(s) + tail)

if __name__ == "__main__":   # kt.py IN OUT IDX VALUE ...   (IDX "lo9" keeps the high 32 bits)
    if sys.argv[1] == "gen": gen(sys.argv[2], int(sys.argv[3])); sys.exit(0)
    if sys.argv[1] == "bgzf": bgzf(sys.argv[2], sys.argv[3]); sys.exit(0)
    if sys.argv[1] == "bgzfx": bgzfx(sys.argv[2]); sys.exit(0)
    d, h = load(sys.argv[1])
    for i, v in zip(sys.argv[3::2], sys.argv[4::2]):
        v = int(v, 0)
        if i.startswith("lo"): j = int(i[2:]); h[j] = (h[j] & ~0xffffffff) | v
        else: h[int(i)] = v
    sign(d, h); open(sys.argv[2], "wb").write(d)
