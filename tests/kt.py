"""Forging helpers for the suite: parse a track, patch it, re-sign its hashes."""
import struct, subprocess, sys
FNV0 = 0xcbf29ce484222325; P = 0x100000001b3; M = (1 << 64) - 1
FLAGS = sum(1 << b for b in (63, 62, 60, 59, 58, 57, 56, 55, 54, 53))
def fnv(h, b):
    for x in b: h = ((h ^ x) * P) & M
    return h
def load(path):
    d = bytearray(open(path, "rb").read()); return d, list(struct.unpack_from("<16Q", d, 8))
def sign(d, h):   # header with h[15] = 0, then lz, nz, iz, sz, itz, np
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
def zst(W, body):   # from a file, so the frame carries its content size
    open(W + "/c", "wb").write(body)
    subprocess.run(["zstd", "-q", "-f", "-19", "--no-check", W + "/c", "-o", W + "/c.zst"], check=True)
    return open(W + "/c.zst", "rb").read()
def unzst(z):
    r = subprocess.run(["zstd", "-d", "-q", "-c"], input=bytes(z), capture_output=True)
    if r.returncode:
        raise RuntimeError("zstd -d failed: " + r.stderr.decode("utf-8", "replace").strip())
    return r.stdout
if __name__ == "__main__":   # kt.py IN OUT IDX VALUE ...   (IDX "lo9" keeps the high 32 bits)
    d, h = load(sys.argv[1])
    for i, v in zip(sys.argv[3::2], sys.argv[4::2]):
        v = int(v, 0)
        if i.startswith("lo"): j = int(i[2:]); h[j] = (h[j] & ~0xffffffff) | v
        else: h[int(i)] = v
    sign(d, h); open(sys.argv[2], "wb").write(d)
