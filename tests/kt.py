"""Forging helpers for the suite: parse a track, patch it, re-sign its hashes."""
import struct, subprocess, sys, random, os
FNV0 = 0xcbf29ce484222325; P = 0x100000001b3; M = (1 << 64) - 1
FLAGS = M ^ ((1 << 49) - 1)   # bits 49-63: flags; below: names length
def fnv(h, b):
    for x in b: h = ((h ^ x) * P) & M
    return h
def load(path):
    d = bytearray(open(path, "rb").read()); return d, list(struct.unpack_from("<16Q", d, 8))
def sign(d, h):   # hash header (h[15]=0), then lz..np sections
    hz = list(h); hz[15] = 0; mh = fnv(FNV0, struct.pack("<16Q", *hz)); off = 136
    for n in (h[6], h[1] & ~FLAGS, h[2], h[3], h[11], h[14]): mh = fnv(mh, d[off:off + n]); off += n
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

if __name__ == "__main__":   # kt.py IN OUT IDX VALUE ...   (IDX "lo9" keeps the high 32 bits)
    if sys.argv[1] == "gen": gen(sys.argv[2], int(sys.argv[3])); sys.exit(0)
    if sys.argv[1] == "bgzf": bgzf(sys.argv[2], sys.argv[3]); sys.exit(0)
    d, h = load(sys.argv[1])
    for i, v in zip(sys.argv[3::2], sys.argv[4::2]):
        v = int(v, 0)
        if i.startswith("lo"): j = int(i[2:]); h[j] = (h[j] & ~0xffffffff) | v
        else: h[int(i)] = v
    sign(d, h); open(sys.argv[2], "wb").write(d)
