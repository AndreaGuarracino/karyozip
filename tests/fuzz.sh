#!/bin/bash
# Sanitizer and leak checks, then generative tests
#   fuzz.sh BINARY [SANITIZER_BINARY] [N]
# N: cases per generator (default 60)
set -u
KT=$(realpath -m "$1"); SAN=$(realpath -m "${2:-$1}"); N=${3:-60}; J=$(nproc)   # forged files run J at a time
for t in python3 zstd; do command -v $t >/dev/null || { echo "fuzz: $t not on PATH"; exit 1; }; done
"$KT" --help >/dev/null 2>&1 || { echo "fuzz: $KT does not run"; exit 1; }   # a binary that cannot run refuses everything
[ -x "$SAN" ] || { echo "fuzz: FAILED, $SAN does not run"; exit 1; }
D=$(cd "$(dirname "$0")" && pwd)
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
fail=0

LK=$SAN
if [ "$SAN" = "$KT" ]; then
  command -v valgrind >/dev/null || { echo "leaks: FAILED, no valgrind for $KT"; exit 1; }
  # suppress glibc TLS of live threads
  printf '{\n  tls\n  Memcheck:Leak\n  match-leak-kinds: definite\n  fun:calloc\n  fun:allocate_dtv\n}\n' > "$W/supp"
  printf '#!/bin/bash\nulimit -n 4096 2>/dev/null\nexec valgrind -q --leak-check=full --errors-for-leak-kinds=definite --error-exitcode=99 --suppressions="%s" "%s" "$@"\n' "$W/supp" "$KT" > "$W/lk"; chmod +x "$W/lk"; LK=$W/lk
fi
lk() {   # fails on refusal unless "(refused)" in DESC
  local what=$1; shift
  printf '  %-46s ' "$what"
  "$@" >/dev/null 2>"$W/err"; local rc=$?
  if grep -qE "LeakSanitizer|definitely lost: [1-9]" "$W/err" || [ $rc = 99 ]; then echo "LEAKED"; fail=1
  elif grep -q "runtime error" "$W/err"; then echo "UNDEFINED BEHAVIOUR"; fail=1   # UBSan does not exit
  elif [ $rc = 0 ] && [[ "$what" != *"(refused)"* ]]; then echo clean
  elif [ $rc = 1 ] && [[ "$what" == *"(refused)"* ]]; then echo "clean (refused)"
  else echo "WRONG EXIT $rc"; fail=1; fi
}
L=$W/lk.d; mkdir -p "$L/out"
printf 'feature_set\tchild\tparent\nfs\tA\tcategorized\nfs\tB\tcategorized\nfs\tC\tcategorized\nfs\tD\tcategorized\nfs\tX\tcategorized\nfs\tZ\tcategorized\n' > "$L/h.tsv"
"$LK" compress < "$D/ls.bed" > "$L/p.kt" 2>/dev/null
"$LK" compress --names drop < "$D/ls.bed" > "$L/s.kt" 2>/dev/null; cut -f1 "$D/ls.bed" | uniq > "$L/s.names"
"$LK" compress --aln "$D/lp.paf" --ref "$D/lref.bed" < "$D/ls.bed" > "$L/a.kt" 2>/dev/null
DB="-i $D/db.acrocentric.bed -i $D/db.chromosome.bed -i $D/db.gene.bed -i $D/db.region.bed -i $D/db.repeat.bed -i $D/db.subtelomeric.bed"
"$LK" compress --db "$D/db.features.tsv" $DB > "$L/f.kt" 2>/dev/null
lk "compress"                       "$LK" compress < "$D/ls.bed"
lk "compress --effort low"          "$LK" --effort low compress < "$D/ls.bed"
lk "compress --block 262144"        "$LK" --block 262144 compress < "$D/ls.bed"
lk "compress --names drop"          "$LK" compress --names drop < "$D/ls.bed"
lk "compress --aln --ref"           "$LK" compress --aln "$D/lp.paf" --ref "$D/lref.bed" < "$D/ls.bed"
lk "compress --db"                  "$LK" compress --db "$D/db.features.tsv" $DB
lk "decompress"                       "$LK" decompress -i "$L/p.kt"
lk "decompress --raw"                 "$LK" decompress --raw -i "$L/p.kt"
lk "decompress --names F"             "$LK" --names "$L/s.names" decompress -i "$L/s.kt"
lk "decompress --aln --ref"           "$LK" decompress --aln "$D/lp.paf" --ref "$D/lref.bed" -i "$L/a.kt"
lk "decompress --db --out"            "$LK" decompress --db "$D/db.features.tsv" --out-dir "$L/out" -i "$L/f.kt"
lk "region, a name that is there" "$LK" query -i "$L/p.kt" -r r1
lk "region, a name that is not (refused)" "$LK" query -i "$L/p.kt" -r nosuchname
lk "region SEQ POS"               "$LK" query -i "$L/p.kt" -r r1:150-151
lk "index"                        "$LK" index -i "$L/p.kt"
lk "region, after index"          "$LK" query -i "$L/p.kt" -r r1:150-151
lk "decompress"                       "$LK" decompress -i "$L/p.kt" --hier "$L/h.tsv" --smooth --set fs
lk "query, one sequence"         "$LK" query -i "$L/p.kt" --hier "$L/h.tsv" --smooth --set fs -r r1
lk "aln-pack"                     "$LK" aln-pack -a "$D/lp.paf" --ref "$D/lref.bed" -o "$L/pk"
lk "--help"                       "$LK" --help
lk "--version"                    "$LK" --version
[ $fail = 0 ] && echo "leaks: clean" || echo "leaks: FOUND SOMETHING"
rm -rf "$L"

# 1. malformed input must be refused
python3 - "$W" "$N" <<'PY'
import random, sys, os
W, N = sys.argv[1], int(sys.argv[2])
random.seed(11)
bad = []
for i in range(N):
    k = random.randrange(7)
    if   k == 0: bad.append("a\t0\t10\tA\na\t20\t30\tB\n")        # gap
    elif k == 1: bad.append("a\t0\t10\tA\na\t9\t20\tB\n")         # overlap
    elif k == 2: bad.append("a\t0\t10\tA\nb\t0\t10\tB\na\t10\t20\tC\n")  # interleaved
    elif k == 3: bad.append("\t0\t10\tA\n")                        # empty name
    elif k == 4: bad.append("a\t0\t0\tA\n")                        # zero length
    elif k == 5: bad.append("a\t0\t%d\tA\n" % (2**64 + 7))         # overflow
    else:        bad.append("a\tx\t10\tA\n")                       # not a number
for i, b in enumerate(bad):
    open(os.path.join(W, "bad%d.bed" % i), "w").write(b)
PY
n=0
for f in "$W"/bad*.bed; do
  $KT compress --block 4 < "$f" > /dev/null 2>&1 && { echo "fuzz: malformed BED ACCEPTED ($f)"; fail=1; }
  n=$((n+1))
done
echo "fuzz: $n malformed encodes, all refused"

# 1b. awkward valid input must round-trip
python3 "$D/kt.py" gen "$W" "$N"
v=0
for f in "$W"/ok*.bed; do
  for m in "" "--block 4" "--block 262144" "--names drop" "--effort low --block 3"; do
    $KT compress $m < "$f" > "$W/v.kt" 2>"$W/err" || { echo "fuzz: valid BED REFUSED ($f $m): $(tail -1 "$W/err")"; fail=1; continue; }
    case "$m" in *drop*)  $KT decompress --raw -i "$W/v.kt" 2>/dev/null | cut -f2- | cmp -s - <(cut -f2- "$f") ;;
                 *)       $KT decompress -i "$W/v.kt" 2>/dev/null | cmp -s - "$f" ;; esac \
      || { echo "fuzz: valid BED did not come back ($f $m)"; fail=1; }
    v=$((v+1))
  done
done
echo "fuzz: $v valid encodes, all came back"

# 2. damaged files: no out-of-bounds reads
python3 - "$W" <<'PY'
import random, sys, os
W = sys.argv[1]; random.seed(12)
with open(os.path.join(W, "ok.bed"), "w") as f:
    for i in range(600):
        p = 0
        for c in sorted(random.sample(range(1, 200), 3)) + [200]:
            f.write("read%d\t%d\t%d\tL%d\n" % (i, p, c, random.randrange(6))); p = c
PY
# decoded rows must match; sanitizer reports fail
check() {
  local f=$1 what=$2 rc out
  out=$(ASAN_OPTIONS=detect_leaks=0 timeout 30 $SAN decompress --raw -i "$f" 2>"$W/err"); rc=$?
  if [ $rc -gt 1 ]; then echo "fuzz: $what crashed rc=$rc"; fail=1; return; fi
  grep -qE "AddressSanitizer|runtime error" "$W/err" && { echo "fuzz: $what sanitizer report"; fail=1; return; }
  if [ $rc -eq 0 ] && [ -s "$W/want" ]; then
    printf '%s\n' "$out" | cut -f2,3,4 | cmp -s - "$W/want" \
 || { echo "fuzz: $what decoded WRONG rows with exit 0"; fail=1; }
  fi
}
d=0
# every writer path
for m in "--block 32" "--block 262144" "--block 32 --effort low" "--block 32 --names drop"; do
  $KT compress $m < "$W/ok.bed" > "$W/t.kt" 2>/dev/null || { echo "fuzz: compress $m FAILED"; fail=1; continue; }
  sz=$(stat -c%s "$W/t.kt")
  $KT decompress --raw -i "$W/t.kt" 2>/dev/null | cut -f2,3,4 > "$W/want"
  for cut in 1 8 64 512; do
    [ "$cut" -ge "$sz" ] && continue
    head -c $((sz - cut)) "$W/t.kt" > "$W/m.kt"
    check "$W/m.kt" "truncation $m cut=$cut"
    d=$((d+1))
  done
  for i in $(seq 1 "$N"); do
    off=$(( (RANDOM * 32768 + RANDOM) % sz ))
    cp "$W/t.kt" "$W/m.kt"
    printf "\\x$(printf %02x $((RANDOM % 256)))" | dd of="$W/m.kt" bs=1 seek=$off conv=notrunc status=none 2>/dev/null
    check "$W/m.kt" "mutation $m off=$off"
    d=$((d+1))
  done
done
# damage --db and --aln tracks too
sets="acrocentric chromosome gene region repeat subtelomeric"; in=""; for f in $sets; do in="$in -i $D/db.$f.bed"; done
$KT compress --db "$D/db.features.tsv" $in > "$W/db.kt" 2>/dev/null
$KT compress --aln "$D/p.paf" --ref "$D/ref.bed" < "$D/s.bed" > "$W/al.kt" 2>/dev/null
for kind in db al; do
  sz=$(stat -c%s "$W/$kind.kt")
  for i in $(seq 1 "$N"); do
    off=$(( (RANDOM * 32768 + RANDOM) % sz ))
    cp "$W/$kind.kt" "$W/m.kt"
    printf "\\x$(printf %02x $((RANDOM % 256)))" | dd of="$W/m.kt" bs=1 seek=$off conv=notrunc status=none 2>/dev/null
    if [ $kind = db ]; then mkdir -p "$W/o"; rm -f "$W/o"/*
      out=$(ASAN_OPTIONS=detect_leaks=0 timeout 30 $SAN decompress --db "$D/db.features.tsv" --out-dir "$W/o" -i "$W/m.kt" 2>"$W/err"); rc=$?
    else
      out=$(ASAN_OPTIONS=detect_leaks=0 timeout 30 $SAN decompress --aln "$D/p.paf" --ref "$D/ref.bed" -i "$W/m.kt" 2>"$W/err"); rc=$?
    fi
    [ $rc -gt 1 ] && { echo "fuzz: $kind mutation crashed rc=$rc off=$off"; fail=1; }
    grep -qE "AddressSanitizer|runtime error" "$W/err" && { echo "fuzz: $kind mutation sanitizer report off=$off"; fail=1; }
    if [ $rc = 0 ]; then   # what decodes must be the input
      if [ $kind = db ]; then for f in $sets; do cmp -s "$W/o/$f.bed" "$D/db.$f.bed" || { echo "fuzz: db mutation decoded WRONG $f off=$off"; fail=1; }; done
      else printf '%s\n' "$out" | cmp -s - "$D/s.bed" || { echo "fuzz: al mutation decoded WRONG rows off=$off"; fail=1; }; fi
    fi
    d=$((d+1))
  done
done
echo "fuzz: $d damaged files, none crashed"

# 2b. damaged packed alignments must never crash
$KT compress --aln "$D/p.paf" --ref "$D/ref.bed" < "$D/s.bed" > "$W/pa.kt" 2>/dev/null
$KT aln-pack -a "$D/p.paf" --ref "$D/ref.bed" -o "$W/pa.pk" >/dev/null 2>&1
$KT decompress --aln "$D/p.paf" --ref "$D/ref.bed" -i "$W/pa.kt" 2>/dev/null > "$W/pa.want"
psz=$(stat -c%s "$W/pa.pk"); pd=0
for cut in $(seq 1 $psz) ; do head -c $((psz - cut)) "$W/pa.pk" > "$W/m.pk"
  out=$(ASAN_OPTIONS=detect_leaks=0 timeout 30 $SAN decompress --aln "$W/m.pk" --ref "$D/ref.bed" -i "$W/pa.kt" 2>"$W/err"); rc=$?
  [ $rc -gt 1 ] && { echo "fuzz: packed alignment cut=$cut crashed rc=$rc"; fail=1; }
  grep -qE "AddressSanitizer|runtime error" "$W/err" && { echo "fuzz: packed alignment cut=$cut sanitizer report"; fail=1; }
  [ $rc -eq 0 ] && { printf '%s\n' "$out" | cmp -s - "$W/pa.want" || { echo "fuzz: packed alignment cut=$cut decoded WRONG rows"; fail=1; }; }
  pd=$((pd+1)); done
for i in $(seq 1 $((N * 2))); do cp "$W/pa.pk" "$W/m.pk"; off=$(( (RANDOM * 32768 + RANDOM) % psz ))
  printf "\\x$(printf %02x $((RANDOM % 256)))" | dd of="$W/m.pk" bs=1 seek=$off conv=notrunc 2>/dev/null
  out=$(ASAN_OPTIONS=detect_leaks=0 timeout 30 $SAN decompress --aln "$W/m.pk" --ref "$D/ref.bed" -i "$W/pa.kt" 2>"$W/err"); rc=$?
  [ $rc -gt 1 ] && { echo "fuzz: packed alignment flip off=$off crashed rc=$rc"; fail=1; }
  grep -qE "AddressSanitizer|runtime error" "$W/err" && { echo "fuzz: packed alignment flip off=$off sanitizer report"; fail=1; }
  [ $rc -eq 0 ] && { printf '%s\n' "$out" | cmp -s - "$W/pa.want" || { echo "fuzz: packed alignment flip off=$off decoded WRONG rows"; fail=1; }; }
  pd=$((pd+1)); done
echo "fuzz: $pd damaged packed alignments, none crashed"

# 3. forged header counts, hash re-signed
$KT compress --block 32 < "$W/ok.bed" > "$W/t3.kt" 2>/dev/null
python3 - "$W" "$D" <<'PY' || { echo "fuzz: stage 3 generator FAILED"; exit 1; }
import struct, sys; sys.path.insert(0, sys.argv[2]); from kt import *
W = sys.argv[1]
src = bytearray(open(W + "/t3.kt", "rb").read()); h0 = list(struct.unpack_from("<16Q", src, 8))
k = 0
for i in range(15):
    for v in (0, 1, 2, 0x7fffffff, 0x80000000, 0xffffffff, 0x100000000, 0x7fffffffffffffff,
              0x8000000000000000, 0xffffffffffffffff, h0[i] ^ (1 << 40)):
        d = bytearray(src); h = list(h0); h[i] = v if i != 1 else (h0[1] & FLAGS) | (v & ~FLAGS)
        try: sign(d, h)
        except Exception: continue
        open("%s/forged%d.kt" % (W, k), "wb").write(d); k += 1
PY
printf 'feature_set\tchild\tparent\nfs\tL0\tcategorized\nfs\tL1\tcategorized\nfs\tL2\tcategorized\nfs\tL3\tcategorized\nfs\tL4\tcategorized\nfs\tL5\tcategorized\n' > "$W/hier.tsv"
one_forged() {   # three commands on one file
  local m=$1 e="$1.err"
  for c in "decompress --raw" "query --raw -r read1" "decompress --smooth --set fs"; do
    case $c in *smooth*) args="-i $m --hier $W/hier.tsv";; *) args="-i $m";; esac
    ASAN_OPTIONS=detect_leaks=0 timeout 30 $SAN $c $args > /dev/null 2>"$e"; rc=$?
    [ $rc -gt 1 ] && { echo "fuzz: forged $(basename $m) $c crashed rc=$rc"; return 1; }
    grep -qE "AddressSanitizer|runtime error" "$e" && { echo "fuzz: forged $(basename $m) $c sanitizer report"; return 1; }
  done
  return 0
}
export -f one_forged; export W KT SAN
f=$(ls "$W"/forged*.kt | wc -l); f=$((f * 3))
ls "$W"/forged*.kt | xargs -P "$J" -I{} bash -c 'one_forged "$@"' _ {} > "$W/forged.out" 2>&1
[ -s "$W/forged.out" ] && { cat "$W/forged.out"; fail=1; }
echo "fuzz: $f forged-header runs, none crashed"

# 4. forged sizes table, re-signed
awk 'BEGIN{for(i=0;i<300;i++){p=0;for(j=0;j<4;j++){c=p+29;printf "r%d\t%d\t%d\t%s\n",i,p,c,(i<50?"A":(j%2?"B":"C"));p=c}}}' > "$W/t4.bed"
$KT compress --block 200 < "$W/t4.bed" > "$W/t4.kt" 2>/dev/null
cut -f2,3,4 "$W/t4.bed" > "$W/want4"
python3 - "$W" "$D" <<'PY' || { echo "fuzz: stage 4 generator FAILED"; exit 1; }
import struct, sys; sys.path.insert(0, sys.argv[2]); from kt import *
W = sys.argv[1]
src = bytearray(open(W + "/t4.kt", "rb").read()); h = list(struct.unpack_from("<16Q", src, 8))
nzl = h[1] & ~FLAGS; o_lz = 136; o_nz = o_lz + h[6]; o_iz = o_nz + nzl; o_sz = o_iz + h[2]; o_end = o_sz + h[3]
tab = unzst(src[o_sz:o_end])
rows = []; i = 0
while i < len(tab):
    r = []
    for k in range(7): v, i = vget(tab, i); r.append(v)
    rows.append(r)
k = 0
for bi in sorted(set([0, len(rows) // 2, len(rows) - 1])):
    for f in range(7):
        for v in (0, 1, 2, 0xffffffff, 0x7fffffffffffffff, 0xffffffffffffffff, rows[bi][f] ^ 1, rows[bi][f] << 1, rows[bi][f] >> 1):
            R = [list(r) for r in rows]; R[bi][f] = v & M
            body = b"".join(vput(x) for r in R for x in r)
            open(W + "/tab", "wb").write(body)   # a file: frame keeps content size
            subprocess.run(["zstd", "-q", "-f", "-19", "--no-check", W + "/tab", "-o", W + "/tab.zst"], check=True)
            z = open(W + "/tab.zst", "rb").read()
            d = bytearray(src[:o_sz]) + z + src[o_end:]
            hh = list(h); hh[3] = len(z); hh[15] = 0
            mh = fnv(FNV0, struct.pack("<16Q", *hh)); off = 136
            for n in (hh[6], nzl, hh[2], hh[3], hh[11]):
                mh = fnv(mh, d[off:off + n]); off += n
            hh[15] = mh; struct.pack_into("<16Q", d, 8, *hh)
            open("%s/fsz%d.kt" % (W, k), "wb").write(d); k += 1
PY
one_fsz() {
  local m=$1 e="$1.err" out rc
  out=$(ASAN_OPTIONS=detect_leaks=0 timeout 30 $SAN decompress --raw -i "$m" 2>"$e"); rc=$?
  [ $rc -gt 1 ] && { echo "fuzz: forged sizes $(basename $m) crashed rc=$rc"; return 1; }
  grep -qE "AddressSanitizer|runtime error" "$e" && { echo "fuzz: forged sizes $(basename $m) sanitizer report"; return 1; }
  [ $rc -eq 0 ] && { printf '%s\n' "$out" | cut -f2,3,4 | cmp -s - "$W/want4" || { echo "fuzz: forged sizes $(basename $m) decoded WRONG rows with exit 0"; return 1; }; }
  return 0
}
export -f one_fsz
g=$(ls "$W"/fsz*.kt | wc -l)
ls "$W"/fsz*.kt | xargs -P "$J" -I{} bash -c 'one_fsz "$@"' _ {} > "$W/fsz.out" 2>&1
[ -s "$W/fsz.out" ] && { cat "$W/fsz.out"; fail=1; }
echo "fuzz: $g forged-sizes runs, none wrong"

# 5. forged index entries and chunk table
awk 'BEGIN{for(i=0;i<70000;i++){printf "s%d\t0\t7\tA\ns%d\t7\t9\tB\n",i,i}}' > "$W/t5.bed"
$KT compress --block 1000 < "$W/t5.bed" > "$W/t5.kt" 2>/dev/null
python3 - "$W" "$D" <<'PY' || { echo "fuzz: stage 5 generator FAILED"; exit 1; }
import struct, sys; sys.path.insert(0, sys.argv[2]); from kt import *
W = sys.argv[1]
src = bytearray(open(W + "/t5.kt", "rb").read()); h = list(struct.unpack_from("<16Q", src, 8))
nzl = h[1] & ~FLAGS; o_nz = 136 + h[6]; o_iz = o_nz + nzl; o_sz = o_iz + h[2]; o_it = o_sz + h[3]; o_np = o_it + h[11]
tab = unzst(src[o_it:o_np])
i = 0; chunks = []
while i < len(tab):
    cl, i = vget(tab, i); sp, i = vget(tab, i); chunks.append([cl, sp])
frames = []; off = 0
for cl, sp in chunks: frames.append(bytes(src[o_iz + off:o_iz + off + cl])); off += cl
def build(frames, chunks, hh):
    iz = b"".join(frames); table = b"".join(vput(c[0]) + vput(c[1]) for c in chunks); it = zst(W, table)
    hh = list(hh); hh[2] = len(iz); hh[11] = len(it); hh[15] = 0
    d = bytearray(src[:o_iz]) + iz + src[o_sz:o_it] + it + src[o_np:]
    mh = fnv(FNV0, struct.pack("<16Q", *hh)); off = 136
    for n in (hh[6], nzl, hh[2], hh[3], hh[11]): mh = fnv(mh, d[off:off + n]); off += n
    hh[15] = mh; struct.pack_into("<16Q", d, 8, *hh); return d
k = 0
ent = unzst(frames[0])
vals = []; i = 0
while i < len(ent): v, i = vget(ent, i); vals.append(v)
for pos in (0, 1, len(vals) // 2, len(vals) - 1):
    for v in (0, 1, 3, 5, 1 << 32, 1 << 40):
        nv = list(vals); nv[pos] = v
        fr = list(frames); fr[0] = zst(W, b"".join(vput(x) for x in nv)); ch = [list(c) for c in chunks]; ch[0][0] = len(fr[0])
        open("%s/fidx%d.kt" % (W, k), "wb").write(build(fr, ch, h)); k += 1
    nv = vals[:pos] + vals[pos + 1:]
    fr = list(frames); fr[0] = zst(W, b"".join(vput(x) for x in nv)); ch = [list(c) for c in chunks]; ch[0][0] = len(fr[0])
    open("%s/fidx%d.kt" % (W, k), "wb").write(build(fr, ch, h)); k += 1
for c in range(len(chunks)):
    for dv in (-1000, -1, 1, 1000, 1 << 40):
        ch = [list(x) for x in chunks]; ch[c][1] = (ch[c][1] + dv) & M
        open("%s/fidx%d.kt" % (W, k), "wb").write(build(frames, ch, h)); k += 1
hh = list(h); hh[0] = h[0] - 1; open("%s/fidx%d.kt" % (W, k), "wb").write(build(frames, chunks, hh)); k += 1
hh = list(h); hh[0] = h[0] + 1; open("%s/fidx%d.kt" % (W, k), "wb").write(build(frames, chunks, hh)); k += 1
hh = list(h); hh[4] = h[4] - 1; open("%s/fidx%d.kt" % (W, k), "wb").write(build(frames, chunks, hh)); k += 1
hh = list(h); hh[8] = 16; open("%s/fidx%d.kt" % (W, k), "wb").write(build(frames, chunks, hh)); k += 1
hh = list(h); hh[13] = 3; open("%s/fidx%d.kt" % (W, k), "wb").write(build(frames, chunks, hh)); k += 1
PY
$KT query -i "$W/t5.kt" -r s69999 2>/dev/null > "$W/want5"
printf 'feature_set\tchild\tparent\nfs\tA\tcategorized\nfs\tB\tcategorized\n' > "$W/hier.tsv"
$KT decompress -i "$W/t5.kt" --hier "$W/hier.tsv" --smooth --set fs 2>/dev/null > "$W/want5s"
one_fidx() {
  local m=$1 e="$1.err" out rc
  out=$(ASAN_OPTIONS=detect_leaks=0 timeout 30 $SAN query -i "$m" -r s69999 2>"$e"); rc=$?
  [ $rc -gt 1 ] && { echo "fuzz: forged index $(basename $m) crashed rc=$rc"; return 1; }
  grep -qE "AddressSanitizer|runtime error" "$e" && { echo "fuzz: forged index $(basename $m) sanitizer report"; return 1; }
  [ $rc -eq 0 ] && { printf '%s\n' "$out" | cmp -s - "$W/want5" || { echo "fuzz: forged index $(basename $m) WRONG rows with exit 0"; return 1; }; }
  out=$(ASAN_OPTIONS=detect_leaks=0 timeout 60 $SAN decompress -i "$m" --hier "$W/hier.tsv" --smooth --set fs 2>"$e"); rc=$?
  [ $rc -gt 1 ] && { echo "fuzz: forged index $(basename $m) decompress crashed rc=$rc"; return 1; }
  [ $rc -eq 0 ] && { printf '%s\n' "$out" | cmp -s - "$W/want5s" || { echo "fuzz: forged index $(basename $m) decompress WRONG rows with exit 0"; return 1; }; }
  return 0
}
export -f one_fidx
x=$(ls "$W"/fidx*.kt | wc -l)
ls "$W"/fidx*.kt | xargs -P "$J" -I{} bash -c 'one_fidx "$@"' _ {} > "$W/fidx.out" 2>&1
[ -s "$W/fidx.out" ] && { cat "$W/fidx.out"; fail=1; }
echo "fuzz: $x forged-index runs, none wrong"
exit $fail
