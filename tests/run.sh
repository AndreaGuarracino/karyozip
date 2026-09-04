#!/bin/bash
# Every case is computed by hand in aln.expected, not by agreeing with another
# implementation: the previous python and C both had the same reverse-strand
# off-by-one, so matching each other proved nothing.
#
#   q1  forward, exact
#   q2  MINUS strand: position 0 is the LAST target base, tp+n-1-pos
#   q3  50M20D50M, the deletion shifts the target
#   q4  50M20I50M, the insertion leaves query bases with no target
#   q5  no alignment at all
set -eu
KT=${KT:-$(dirname $0)/../build/karyozip}
D=$(dirname "$0")
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
fail=0

$KT encode --aln "$D/p.paf" --ref "$D/ref.bed" < "$D/s.bed" > "$W/d.kt" 2>/dev/null
$KT decode --raw "$W/d.kt" 2>/dev/null | cut -f2,3,4 > "$W/got"
cut -f2,3,4 "$D/aln.expected" > "$W/want"
cmp -s "$W/got" "$W/want" && echo "aln delta: correct" || { echo "aln delta: WRONG"; diff "$W/got" "$W/want" | head; fail=1; }

$KT decode --aln "$D/p.paf" --ref "$D/ref.bed" "$W/d.kt" 2>/dev/null > "$W/rt"
cmp -s "$W/rt" "$D/s.bed" && echo "aln restore: byte-identical" || { echo "aln restore: WRONG"; fail=1; }

$KT encode < "$D/s.bed" > "$W/p.kt" 2>/dev/null
$KT decode "$W/p.kt" 2>/dev/null | cmp -s - "$D/s.bed" && echo "plain: byte-identical" || { echo "plain: WRONG"; fail=1; }

for f in --sort "--sort --names drop" "--block all" "--block 4"; do
  ofl=""
  $KT encode $f < "$D/s.bed" > "$W/f.kt" 2>/dev/null || { echo "$f: ENCODE FAILED"; fail=1; continue; }
  case "$f" in
    # a dropped-names track refuses a plain decode: the stored ordinals are
    # not names, so asking for them is explicit
    *"--names drop"*|*--sort*)
      case "$f" in *"--names drop"*) ofl="--names ordinal" ;; esac
      $KT decode ${ofl:-} "$W/f.kt" 2>/dev/null | sort -k2 > "$W/a"
                           sort -k2 "$D/s.bed" | cut -f2,3,4 > "$W/b"
                           cut -f2,3,4 "$W/a" > "$W/a2"
                           cmp -s "$W/a2" "$W/b" && echo "$f: same rows" || { echo "$f: WRONG"; fail=1; } ;;
    *) $KT decode "$W/f.kt" 2>/dev/null | cmp -s - "$D/s.bed" && echo "$f: byte-identical" || { echo "$f: WRONG"; fail=1; } ;;
  esac
done
# --db: 6 BEDs in, one featureID track, 6 BEDs out, each byte compared. Two
# featureIDs can carry the same label in one set, so the merge on decode must
# compare the label TEXT; comparing pointers silently stopped merging.
sets="acrocentric chromosome gene region repeat subtelomeric"
in=""; for f in $sets; do in="$in $D/db.$f.bed"; done
$KT encode --db "$D/db.features.tsv" $in > "$W/f.kt" 2>/dev/null
mkdir -p "$W/o"
$KT decode --db "$D/db.features.tsv" --out "$W/o" "$W/f.kt" 2>/dev/null
n=0; for f in $sets; do cmp -s "$W/o/$f.bed" "$D/db.$f.bed" && n=$((n+1)); done
[ "$n" = 6 ] && echo "db: 6/6 identical" || { echo "db: only $n/6 identical"; fail=1; }

# give byte-identical output to the in-memory sort: that is
# the whole point of an external sort.
$KT encode --sort --names drop < "$D/s.bed" > "$W/mem.kt" 2>/dev/null
mkdir -p "$W/tmp"
$KT encode --tmp-dir "$W/tmp" --sort=1M --names drop < "$D/s.bed" > "$W/ext.kt" 2>/dev/null
cmp -s "$W/mem.kt" "$W/ext.kt" && echo "external sort: identical to in-memory" \
 || { echo "external sort: DIFFERS from in-memory"; fail=1; }
[ -z "$(ls -A "$W/tmp")" ] && echo "external sort: temp files cleaned up" \
 || { echo "external sort: temp files LEFT BEHIND"; fail=1; }

# --db with sequence names of DIFFERENT LENGTHS. encode_db cached the name
# length behind a pointer compare against a fixed buffer, so it froze on the
# first sequence and every later name was truncated or overran. The original
# db fixture used one name length and never caught it.
in2=""; for f in $sets; do in2="$in2 $D/db2.$f.bed"; done
$KT encode --db "$D/db2.features.tsv" $in2 > "$W/f2.kt" 2>/dev/null
mkdir -p "$W/o2"
$KT decode --db "$D/db2.features.tsv" --out "$W/o2" "$W/f2.kt" 2>/dev/null
n=0; for f in $sets; do cmp -s "$W/o2/$f.bed" "$D/db2.$f.bed" && n=$((n+1)); done
[ "$n" = 6 ] && echo "db mixed name lengths: 6/6 identical" \
 || { echo "db mixed name lengths: only $n/6"; fail=1; }

# A track written with a transform must not decode without it. Every one of
# these used to exit 0 and emit "=" rows, bare featureIDs, or fabricated
# ordinals: valid-looking output that was silently wrong.
must_fail() {  # description, then the command
  local what=$1; shift
  if "$@" >/dev/null 2>&1; then echo "$what: ACCEPTED, should have failed"; fail=1
  else echo "$what: refused, as it should"; fi
}
must_fail "decode --aln track without --aln" $KT decode "$W/d.kt"
must_fail "decode plain track with --aln" $KT decode --aln "$D/p.paf" --ref "$D/ref.bed" "$W/p.kt"
must_fail "decode --db track without --db" $KT decode "$W/f.kt"
must_fail "encode --ref without --aln" $KT encode --ref "$D/ref.bed"
must_fail "--aln-names without --aln" $KT encode --aln-names "$D/ord.names.txt" < "$D/s.bed"
must_fail "a name file given to encode as --names" $KT encode --names "$D/ord.names.txt" < "$D/s.bed"
must_fail "--sort= with no size" $KT encode --sort=
must_fail "the retired --sort-tmp-dir" $KT encode --sort-tmp-dir "$W"

# the sidecar --names drop:writes must be the one --names reads
$KT encode --names drop:"$W/side.bin" < "$D/s.bed" > "$W/sc.kt" 2>/dev/null
$KT decode --names "$W/side.bin" "$W/sc.kt" 2>/dev/null | cmp -s - "$D/s.bed" \
  && echo "sidecar: round-trips byte-identical" \
 || { echo "sidecar: WRONG"; fail=1; }
must_fail "decode a sidecar track without --names" $KT decode "$W/sc.kt"
# the delta is computed from the names, so a nameless aln track cannot be read
must_fail "encode --aln with --names drop and no name list" \
  $KT encode --aln "$D/p.paf" --ref "$D/ref.bed" --names drop
# sorting reorders records, so an ordinal name map names the wrong ones
must_fail "encode --aln --sort --names drop even with a name list" \
  $KT encode --aln "$D/p.paf" --ref "$D/ref.bed" --names "$D/s.bed" --sort --names drop

# --db with --aln. Two touching runs can both become "=" while the featureIDs
# behind them differ, and a guard that compared the LABEL rejected every such
# encode. The check belongs on the id, before the delta.
$KT encode --db "$D/db.features.tsv" --aln "$D/db.self.paf" --ref "$D/db.rall.bed" $in \
  > "$W/da.kt" 2>/dev/null
mkdir -p "$W/dao"
$KT decode --db "$D/db.features.tsv" --aln "$D/db.self.paf" --ref "$D/db.rall.bed" \
  --out "$W/dao" "$W/da.kt" 2>/dev/null
n=0; for f in $sets; do cmp -s "$W/dao/$f.bed" "$D/db.$f.bed" && n=$((n+1)); done
[ "$n" = 6 ] && echo "db+aln: 6/6 identical" || { echo "db+aln: only $n/6"; fail=1; }

# HKS tracks name sequences by ordinal, so --aln reaches the PAF only through
# --names FILE. Dropping the names is legal here and must still round-trip.
$KT encode --aln "$D/p.paf" --ref "$D/ref.bed" --aln-names "$D/ord.names.txt" \
  --names drop < "$D/ord.bed" > "$W/ord.kt" 2>/dev/null
$KT decode --aln "$D/p.paf" --ref "$D/ref.bed" --aln-names "$D/ord.names.txt" \
  "$W/ord.kt" 2>/dev/null > "$W/ord.out"
cmp -s "$D/ord.bed" "$W/ord.out" \
  && echo "aln --names drop with an ordinal map: round-trips" \
 || { echo "aln --names drop with an ordinal map: DIFFERS"; fail=1; }

# A read that hit no k-mer leaves no rows, so the ordinals in the track have
# gaps and the rank no longer equals the ordinal. Dropping the names then
# pointed --aln at the wrong read and decode failed on real ONT data.
$KT encode --aln "$D/p.paf" --ref "$D/ref.bed" --aln-names "$D/ord2.names.txt" \
  --names drop < "$D/ord2.bed" > "$W/ord2.kt" 2>/dev/null
$KT decode --aln "$D/p.paf" --ref "$D/ref.bed" --aln-names "$D/ord2.names.txt" \
  "$W/ord2.kt" 2>/dev/null > "$W/ord2.out"
cmp -s "$D/ord2.bed" "$W/ord2.out" \
  && echo "aln --names drop with gaps in the ordinals: round-trips" \
 || { echo "aln --names drop with gaps in the ordinals: DIFFERS"; fail=1; }

# A run that ends where the reference run ends stores no length at all. This
# fixture is built so 6 of its 8 rows do, which the s.bed fixture never did.
$KT encode --aln "$D/lp.paf" --ref "$D/lref.bed" --names keep < "$D/ls.bed" \
  > "$W/ls.kt" 2>/dev/null
$KT decode --aln "$D/lp.paf" --ref "$D/lref.bed" "$W/ls.kt" 2>/dev/null > "$W/ls.out"
cmp -s "$D/ls.bed" "$W/ls.out" \
  && echo "lengths the reference implies: round-trip" \
 || { echo "lengths the reference implies: DIFFERS"; fail=1; }
must_fail "decode a length-delta track without --aln" $KT decode "$W/ls.kt"
# Each block keeps the sentinel only where it shrinks the stream, so an 8-row
# file may store plain lengths. This one cannot: 3000 runs of random length
# that all end on a reference run end, so the sentinel column is all zeros.
awk 'BEGIN{srand(5);p=0;for(i=0;i<3000;i++){c=p+int(rand()*900)+1;printf "chrA\t%d\t%d\tL%d\n",p,c,i%4;p=c}print p > "'"$W"'/len"}' > "$W/big.ref"
L=$(cat "$W/len"); sed 's/^chrA/q/' "$W/big.ref" > "$W/big.q"
printf "q\t$L\t0\t$L\t+\tchrA\t$L\t0\t$L\t$L\t$L\t60\tcg:Z:${L}M\n" > "$W/big.paf"
$KT encode --aln "$W/big.paf" --ref "$W/big.ref" < "$W/big.q" > "$W/big.kt" 2>/dev/null
n=$($KT decode --raw "$W/big.kt" 2>/dev/null | awk -F'\t' '$2==$3' | wc -l)
[ "$n" = 3000 ] && echo "lengths left out: 3000 of 3000 rows" \
 || { echo "lengths left out: $n of 3000, expected 3000"; fail=1; }
$KT decode --aln "$W/big.paf" --ref "$W/big.ref" "$W/big.kt" 2>/dev/null | cmp -s - "$W/big.q" \
  && echo "sentinel block: round-trip" || { echo "sentinel block: DIFFERS"; fail=1; }
# and a block where plain lengths win keeps them, flagged in the sizes table
awk 'BEGIN{srand(6);p=0;for(i=0;i<3000;i++){c=p+int(rand()*900)+1;printf "q\t%d\t%d\tL%d\n",p,c,i%4;p=c}}' > "$W/off.q"
$KT encode --aln "$W/big.paf" --ref "$W/big.ref" < "$W/off.q" > "$W/off.kt" 2>/dev/null
$KT decode --aln "$W/big.paf" --ref "$W/big.ref" "$W/off.kt" 2>/dev/null | cmp -s - "$W/off.q" \
  && echo "plain-length block under --aln: round-trip" || { echo "plain-length block under --aln: DIFFERS"; fail=1; }

# Reads are all one length, so the last run of each is the total minus the
# rest. Storing that total as a delta costs nothing where a length would.
python3 - > "$W/tail.bed" <<'PY'
import random
random.seed(5)
for s in range(4000):
    cuts = sorted(random.sample(range(1, 116), 3)); p = 0
    for c in cuts + [116]:
        print(f"r{s}\t{p}\t{c}\tL{random.randint(0,5)}"); p = c
PY
a=$($KT --tail 0 encode --block 4096 < "$W/tail.bed" 2>/dev/null | wc -c)
b=$($KT --tail 1 encode --block 4096 < "$W/tail.bed" 2>/dev/null | wc -c)
$KT encode --block 4096 < "$W/tail.bed" > "$W/tail.kt" 2>/dev/null
$KT decode "$W/tail.kt" 2>/dev/null > "$W/tail.out"
cmp -s "$W/tail.bed" "$W/tail.out" \
  && echo "constant-length reads: round-trip" \
 || { echo "constant-length reads: DIFFERS"; fail=1; }
[ "$b" -lt "$a" ] && echo "constant-length reads: $a -> $b bytes" \
 || { echo "constant-length reads: no gain, $a -> $b"; fail=1; }
# random access must still land on the right rows
$KT region "$W/tail.kt" r2000 2>/dev/null > "$W/tq.out"
grep "^r2000	" "$W/tail.bed" > "$W/tq.exp"
cmp -s "$W/tq.out" "$W/tq.exp" && echo "constant-length reads: query exact" \
 || { echo "constant-length reads: query WRONG"; fail=1; }

# Trimmed reads vary in length, and then storing the total LOSES: measured
# -9.4% before the gate looked at the totals instead of the run count.
python3 - > "$W/vlen.bed" <<'PY'
import random
random.seed(5)
for s in range(4000):
    L = 116 - random.randint(0, 80)
    cuts = sorted(random.sample(range(1, L), 3)); p = 0
    for c in cuts + [L]:
        print(f"r{s}\t{p}\t{c}\tL{random.randint(0,5)}"); p = c
PY
a=$($KT --tail 0 encode --block 4096 < "$W/vlen.bed" 2>/dev/null | wc -c)
b=$($KT encode --block 4096 < "$W/vlen.bed" 2>/dev/null | wc -c)
[ "$a" = "$b" ] && echo "reads of mixed length: left alone, $a bytes" \
 || { echo "reads of mixed length: changed $a -> $b, gate is wrong"; fail=1; }
$KT encode --block 4096 < "$W/vlen.bed" > "$W/vlen.kt" 2>/dev/null
$KT decode "$W/vlen.kt" 2>/dev/null | cmp -s - "$W/vlen.bed" \
  && echo "reads of mixed length: round-trip" \
 || { echo "reads of mixed length: DIFFERS"; fail=1; }

# Every real short-read track is one length at 100.0% but still carries a few
# stray reads. A min/max gate threw the whole win away over one of them.
python3 - > "$W/out.bed" <<'PY'
import random
random.seed(11)
for s in range(4000):
    L = random.choice([8, 30, 69, 105]) if random.randint(0, 50) == 0 else 116
    cuts = sorted(random.sample(range(1, L), 3)); p = 0
    for c in cuts + [L]:
        print(f"r{s}\t{p}\t{c}\tL{random.randint(0,5)}"); p = c
PY
a=$($KT --tail 0 encode --block 4096 < "$W/out.bed" 2>/dev/null | wc -c)
b=$($KT encode --block 4096 < "$W/out.bed" 2>/dev/null | wc -c)
[ "$b" -lt "$a" ] && echo "uniform reads with stray lengths: $a -> $b bytes" \
 || { echo "uniform reads with stray lengths: gate gave up, $a -> $b"; fail=1; }
$KT encode --block 4096 < "$W/out.bed" > "$W/out.kt" 2>/dev/null
$KT decode "$W/out.kt" 2>/dev/null | cmp -s - "$W/out.bed" \
  && echo "uniform reads with stray lengths: round-trip" \
 || { echo "uniform reads with stray lengths: DIFFERS"; fail=1; }

# A row belongs to the record above it only if it carries the same name. The
# sort path keyed on start == 0 alone and filed this row under the wrong
# sequence, while the plain path rejected it.
printf 'a\t0\t10\tA\nb\t0\t10\tB\na\t10\t20\tC\n' > "$W/il.bed"
$KT encode --block 4 < "$W/il.bed" > /dev/null 2>&1 \
  && { echo "interleaved records, plain: ACCEPTED"; fail=1; } \
 || echo "interleaved records, plain: refused, as it should"
$KT encode --block 4 --tmp-dir "$W" --sort=1G < "$W/il.bed" > /dev/null 2>&1 \
  && { echo "interleaved records, --sort: ACCEPTED"; fail=1; } \
 || echo "interleaved records, --sort: refused, as it should"

# The rows tile each sequence gaplessly and that is what lets start go
# unstored. --sort took a gap and an overlap and moved the run instead.
printf 'a\t0\t10\tA\na\t20\t30\tB\n' > "$W/gap.bed"
printf 'a\t0\t10\tA\na\t9\t20\tB\n'  > "$W/ovl.bed"
for f in gap ovl; do
  for m in "" "--tmp-dir $W --sort=1G"; do
    t=plain; [ -n "$m" ] && t=--sort
    $KT encode --block 4 $m < "$W/$f.bed" > /dev/null 2>&1 \
      && { echo "$f, $t: ACCEPTED"; fail=1; } \
 || echo "$f, $t: refused, as it should"
  done
done

# --names keep on names that are already 0,1,2,... wrote a file it refused to read
printf '0\t0\t10\tA\n1\t0\t10\tB\n2\t0\t10\tC\n' > "$W/ord2.bed"
$KT encode --block 4 --names keep < "$W/ord2.bed" > "$W/ord2.kt" 2>/dev/null
$KT decode "$W/ord2.kt" 2>/dev/null | cmp -s - "$W/ord2.bed" \
  && echo "--names keep on ordinals: round-trips" \
 || { echo "--names keep on ordinals: DIFFERS"; fail=1; }

# a single line with no trailing newline was thrown away and still exited 0
printf 'a\t0\t10\tA' > "$W/nonl.bed"
$KT encode --block 4 < "$W/nonl.bed" > /dev/null 2>&1 \
  && { echo "input without a trailing newline: ACCEPTED"; fail=1; } \
 || echo "input without a trailing newline: refused, as it should"

timeout 20 $KT --tok-chunk 0 encode --block 4 < "$W/ord2.bed" > /dev/null 2>&1 \
  && { echo "--tok-chunk 0: ACCEPTED"; fail=1; } \
 || echo "--tok-chunk 0: refused, as it should"

# The --db parser used strtoull and took "x" as 0, an overflow as UINT64_MAX
# and a reverse interval as valid, while the plain parser rejected all three.
for bad in 'a\tx\t10\t1' '\t0\t10\t1' 'a\t0\t18446744073709551616\t1' 'a\t0\t20\t1\na\t20\t10\t2'; do
  printf "$bad\n" > "$W/dbbad.bed"
  $KT encode --db "$D/db.features.tsv" "$W/dbbad.bed" > /dev/null 2>&1 \
    && { echo "--db malformed input: ACCEPTED"; fail=1; } \
 || echo "--db malformed input: refused, as it should"
done

# smooth crashed with a double free on every real cell: the threaded path gives
# each worker the primary reader's name table and each one freed it. The
# harness reported it as a bare word no monitor counted.
# 10,000 sequences: the threaded path only splits work above 4,096, so 200
# never exercised it and a shared decode boundary between workers went unseen.
python3 - > "$W/sm.bed" <<'PY'
for i in range(10000):
    print("r%d\t0\t50\tacro" % i); print("r%d\t50\t100\tnovel" % i)
PY
printf 'feature_set\tchild\tparent\nacrocentric\tacro\tcategorized\nacrocentric\tnovel\tcategorized\n' > "$W/sm.tsv"
$KT encode --block 32 < "$W/sm.bed" > "$W/sm.kt" 2>/dev/null
$KT --threads 1 smooth "$W/sm.kt" "$W/sm.tsv" acrocentric > "$W/sm.1" 2>/dev/null || { echo "smooth, serial: FAILED"; fail=1; }
for t in 2 4 8; do
  $KT --threads $t smooth "$W/sm.kt" "$W/sm.tsv" acrocentric > "$W/sm.$t" 2>/dev/null \
 || { echo "smooth, $t threads: FAILED"; fail=1; continue; }
  cmp -s "$W/sm.1" "$W/sm.$t" && echo "smooth, $t threads: identical to serial" \
 || { echo "smooth, $t threads: DIFFERS from serial"; fail=1; }
done

# A CIGAR was trusted: "M" read as 0M, an overflow wrapped, "Q" was ignored,
# and consumed length was never checked against the record. Each one placed
# a delta somewhere else without a word.
T1=$(cut -f1 "$D/ref.bed" | head -1); printf 'a\t0\t10\tA\n' > "$W/cg.bed"
for cg in "cg:Z:M" "cg:Z:999999999999999999999999M" "cg:Z:5M5Q5M" "cg:Z:5M" "cg:Z:0M10M"; do
  printf "a\t10\t0\t10\t+\t$T1\t100\t0\t10\t10\t10\t60\t$cg\n" > "$W/cg.paf"
  $KT encode --aln "$W/cg.paf" --ref "$D/ref.bed" < "$W/cg.bed" > /dev/null 2>&1 \
    && { echo "bad CIGAR $cg: ACCEPTED"; fail=1; } || echo "bad CIGAR $cg: refused, as it should"
done
printf "a\t10\t0\t10\t+\t$T1\t100\t0\t10\t10\t10\t60\tcg:Z:10M\n" > "$W/cg.paf"
$KT encode --aln "$W/cg.paf" --ref "$D/ref.bed" < "$W/cg.bed" > /dev/null 2>&1 \
  && echo "good CIGAR: accepted" || { echo "good CIGAR: REFUSED"; fail=1; }

# An audit found these. A rows-per-byte header check refused a repetitive
# track; a features.tsv defining id 0 decoded as novel; a FASTQ read longer
# than 65535 bp desynchronised the name file; a stray <track>.names renamed an
# ordinal track; a label longer than 255 bytes overflowed the emit buffer.
awk 'BEGIN{for(i=0;i<20000;i++)printf "s\t%d\t%d\tA\n",i*10,i*10+10}' > "$W/rep.bed"
$KT encode < "$W/rep.bed" > "$W/rep.kt" 2>/dev/null
$KT decode "$W/rep.kt" 2>/dev/null | cmp -s - "$W/rep.bed" && echo "repetitive track: byte-identical" \
 || { echo "repetitive track: WRONG"; fail=1; }
$KT smooth "$W/rep.kt" "$W/sm.tsv" acrocentric > /dev/null 2>&1 && echo "repetitive smooth: ok" \
 || { echo "repetitive smooth: FAILED"; fail=1; }
{ head -1 "$D/db.features.tsv"; printf '0\tacro\tchr1\tintergenic\tarm\tLINE\tITS\n'; tail -n +2 "$D/db.features.tsv"; } > "$W/id0.tsv"
must_fail "features.tsv defining id 0" $KT encode --db "$W/id0.tsv" $in
LR=$(head -c 70000 /dev/zero | tr '\0' A)
printf 'long\t0\t10\tA\nnext\t0\t10\tB\n' > "$W/lr.bed"
printf '@long\n%s\n+\n%s\n@next\nACGT\n+\nIIII\n' "$LR" "$LR" > "$W/lr.fq"
$KT encode --sort --names drop < "$W/lr.bed" > "$W/lr.kt" 2>/dev/null
$KT --names "$W/lr.fq" decode "$W/lr.kt" 2>/dev/null | cmp -s - "$W/lr.bed" \
  && echo "long FASTQ read: names restored" || { echo "long FASTQ read: WRONG"; fail=1; }
printf '0\t0\t10\tA\n1\t0\t10\tB\n' > "$W/ord.bed"
$KT encode < "$W/ord.bed" > "$W/ord.kt" 2>/dev/null; printf 'x\ny\n' > "$W/ord.kt.names"
$KT decode "$W/ord.kt" 2>/dev/null | cmp -s - "$W/ord.bed" && echo "stray sidecar: ignored" \
 || { echo "stray sidecar: RENAMED"; fail=1; }
LL=$(head -c 300 /dev/zero | tr '\0' L)_specific
printf 's\t0\t10\t%s\n' "$LL" > "$W/ll.bed"
$KT encode < "$W/ll.bed" > "$W/ll.kt" 2>/dev/null
$KT decode "$W/ll.kt" > /dev/null 2>&1 && echo "long label v1: ok" || { echo "long label v1: CRASHED"; fail=1; }

# A second audit. --sort --aln looked up the PAF by rank, not name; a short
# features.tsv row left NULL labels that decode dereferenced; a wrong sidecar
# renamed every row; a plain name list lost everything after a space; and
# `encode FILE` truncated FILE as a sidecar before reading a byte.
printf 'chr1\t0\t100\tA\nchr1\t100\t200\tB\nchr1\t200\t300\tC\nchr1\t300\t400\tD\nchr1\t400\t500\tE\nchr1\t500\t600\tF\n' > "$W/sa.ref"
printf '1\t500\t0\t500\t+\tchr1\t600\t0\t550\t500\t550\t60\tcg:Z:100M50D400M\n0\t500\t0\t500\t+\tchr1\t600\t0\t500\t500\t500\t60\tcg:Z:500M\n' > "$W/sa.paf"
printf '1\t0\t100\tA\n1\t100\t200\tB\n1\t200\t300\tC\n1\t300\t400\tD\n1\t400\t500\tE\n0\t0\t100\tA\n0\t100\t200\tB\n0\t200\t300\tC\n0\t300\t400\tD\n0\t400\t500\tE\n' > "$W/sa.bed"
$KT encode --sort --aln "$W/sa.paf" --ref "$W/sa.ref" < "$W/sa.bed" > "$W/sa.kt" 2>/dev/null
$KT decode --aln "$W/sa.paf" --ref "$W/sa.ref" "$W/sa.kt" 2>/dev/null | sort | cmp -s - <(sort "$W/sa.bed") \
  && echo "sort with aln: same rows" || { echo "sort with aln: WRONG"; fail=1; }
printf 'featureID\tacrocentric\tchromosome\n1\tacro\tchr1\n2\tacro\n' > "$W/short.tsv"
printf 's\t0\t10\t1\ns\t10\t20\t2\n' > "$W/ids.bed"
must_fail "features.tsv row with a missing label" $KT encode --db "$W/short.tsv" "$W/ids.bed"
printf 'x1\t0\t10\tA\nx2\t0\t10\tB\n' > "$W/t1.bed"; printf 'y1\t0\t10\tA\ny2\t0\t10\tB\n' > "$W/t2.bed"
$KT encode --names drop:"$W/t1.names" < "$W/t1.bed" > "$W/t1.kt" 2>/dev/null
$KT encode --names drop:"$W/t2.names" < "$W/t2.bed" > "$W/t2.kt" 2>/dev/null
must_fail "sidecar from another track" $KT decode --names "$W/t2.names" "$W/t1.kt"
$KT decode --names "$W/t1.names" "$W/t1.kt" 2>/dev/null | cmp -s - "$W/t1.bed" \
  && echo "own sidecar: byte-identical" || { echo "own sidecar: WRONG"; fail=1; }
printf 'read one\t0\t10\tA\nread two\t0\t10\tB\n' > "$W/sp.bed"; printf 'read one\nread two\n' > "$W/sp.names.txt"
$KT encode --names drop < "$W/sp.bed" > "$W/sp.kt" 2>/dev/null
$KT --names "$W/sp.names.txt" decode "$W/sp.kt" 2>/dev/null | cmp -s - "$W/sp.bed" \
  && echo "name list with spaces: restored" || { echo "name list with spaces: WRONG"; fail=1; }
cp "$D/s.bed" "$W/keep.bed"; $KT encode "$W/keep.bed" < /dev/null > /dev/null 2>&1 && { echo "encode FILE: ACCEPTED"; fail=1; }
cmp -s "$W/keep.bed" "$D/s.bed" && echo "encode FILE: refused, input intact" || { echo "encode FILE: INPUT DESTROYED"; fail=1; }

# A third audit. region lost halo runs in an earlier block; a stale .idx
# answered for another sequence; an input label "=" under --aln decoded as
# the reference label; a sequence total past 2^63 broke the tail predictor;
# a cyclic hierarchy hung smooth; a full disk exited 0; a NUL in a row was
# cut; spill files outlived a failed encode; a different features.tsv at
# decode wrote wrong BEDs.
awk 'BEGIN{for(i=0;i<12;i++)printf "s\t%d\t%d\tL%d\n",i*10,i*10+10,i}' > "$W/halo.bed"
$KT encode --block 4 < "$W/halo.bed" > "$W/halo4.kt" 2>/dev/null; $KT encode < "$W/halo.bed" > "$W/halo1.kt" 2>/dev/null
$KT region "$W/halo4.kt" s 45 46 3 2>/dev/null | cmp -s - <($KT region "$W/halo1.kt" s 45 46 3 2>/dev/null) \
  && [ "$($KT region "$W/halo4.kt" s 45 46 3 2>/dev/null | wc -l)" = 7 ] \
  && echo "halo across blocks: 7 rows" || { echo "halo across blocks: WRONG"; fail=1; }
awk 'BEGIN{for(i=0;i<20;i++)printf "r%d\t0\t%d\tL%d\n",i,50+i,i%5}' > "$W/idx.bed"
$KT encode < "$W/idx.bed" > "$W/idx.kt" 2>/dev/null; $KT index "$W/idx.kt" >/dev/null 2>&1
$KT encode --sort --names keep < "$W/idx.bed" > "$W/idx.kt" 2>/dev/null
$KT region "$W/idx.kt" r7 2>/dev/null | cmp -s - <(grep -P '^r7\t' "$W/idx.bed") \
  && echo "stale idx: ignored" || { echo "stale idx: WRONG ROWS"; fail=1; }
$KT index "$W/idx.kt" >/dev/null 2>&1; $KT region "$W/idx.kt" r7 2>/dev/null | cmp -s - <(grep -P '^r7\t' "$W/idx.bed") \
  && echo "fresh idx: correct" || { echo "fresh idx: WRONG ROWS"; fail=1; }
T1=$(cut -f1 "$D/ref.bed" | head -1)
printf 'a\t0\t5\t=\na\t5\t10\tZ\n' > "$W/eq.bed"
printf "a\t10\t0\t10\t+\t$T1\t100\t0\t10\t10\t10\t60\tcg:Z:10M\n" > "$W/eq.paf"
must_fail "label = under --aln" $KT encode --aln "$W/eq.paf" --ref "$D/ref.bed" < "$W/eq.bed"
python3 -c "
for i in range(70): print(f's{i}\t0\t100\tA'); print(f's{i}\t100\t{2**63+100}\tB')" > "$W/t63.bed"
$KT encode < "$W/t63.bed" > "$W/t63.kt" 2>/dev/null
$KT decode "$W/t63.kt" 2>/dev/null | cmp -s - "$W/t63.bed" && echo "total past 2^63: byte-identical" \
 || { echo "total past 2^63: WRONG"; fail=1; }
printf 'feature_set\tchild\tparent\nacrocentric\tacro\tx\nacrocentric\tx\tacro\n' > "$W/cyc.tsv"
rc=0; timeout 10 $KT smooth "$W/p.kt" "$W/cyc.tsv" acrocentric >/dev/null 2>&1 || rc=$?
[ $rc = 1 ] && echo "cyclic hierarchy: refused" || { echo "cyclic hierarchy: rc=$rc"; fail=1; }
for c in "decode $W/p.kt" "region $W/p.kt q1" "smooth $W/sm.kt $W/sm.tsv acrocentric" "decode --raw $W/p.kt"; do
  $KT $c > /dev/full 2>/dev/null && { echo "$c to /dev/full: EXIT 0"; fail=1; } || echo "${c%% *} to /dev/full: exit 1"
done
printf 'a\t0\t10\tA\0B\nc\t0\t10\tB\n' > "$W/nul.bed"
must_fail "NUL inside a row" $KT encode < "$W/nul.bed"
mkdir -p "$W/sptmp"; { awk 'BEGIN{for(i=0;i<80000;i++)printf "s%d\t0\t10\tA\n",i}'; printf 'zz\t5\t10\tL0\n'; } > "$W/splate.bed"
$KT encode --tmp-dir "$W/sptmp" --sort=1M < "$W/splate.bed" > /dev/null 2> "$W/sp.err" && { echo "late overlap: ACCEPTED"; fail=1; }
grep -qi spill "$W/sp.err" || { echo "late overlap: NEVER SPILLED"; fail=1; }
[ "$(ls "$W/sptmp" | wc -l)" = 0 ] && echo "spill files after a failure: none" || { echo "spill files after a failure: $(ls "$W/sptmp" | wc -l) LEFT"; fail=1; }
sed '2s/acro/ACRO/' "$D/db.features.tsv" > "$W/other.tsv"
mkdir -p "$W/o2"; must_fail "decode --db with another features.tsv" $KT decode --db "$W/other.tsv" --out "$W/o2" "$W/f.kt"

# A block that keeps plain lengths compressed its ids from the wrong offset:
# every label in the block landed on another row, exit 0. These inputs make
# some blocks pick plain lengths and some the sentinel.
for seed in 1 3; do
python3 - "$W" "$seed" <<'PY'
import random, sys
W, seed = sys.argv[1], int(sys.argv[2]); random.seed(seed)
ref = []; p = 0
for i in range(3000):
    c = p + random.randint(1, 900); ref.append((p, c, "L%d" % (i % 4))); p = c
L = p
open(W + "/mx.ref", "w").write("".join("chrA\t%d\t%d\t%s\n" % r for r in ref))
paf = []; bed = []
for k in range(6):
    paf.append("s%d\t%d\t0\t%d\t%s\tchrA\t%d\t0\t%d\t%d\t%d\t60\tcg:Z:%dM\n" % (k, L, L, "-" if k % 2 else "+", L, L, L, L, L))
    p = 0
    for i in range(3000):
        c = p + (random.randint(1, 60) if i < 8 else random.randint(1, 900))
        bed.append("s%d\t%d\t%d\tL%d\n" % (k, p, c, random.randint(0, 3))); p = c
open(W + "/mx.paf", "w").write("".join(paf)); open(W + "/mx.bed", "w").write("".join(bed))
PY
$KT encode --aln "$W/mx.paf" --ref "$W/mx.ref" --block 3000 < "$W/mx.bed" > "$W/mx.kt" 2>/dev/null
$KT decode --aln "$W/mx.paf" --ref "$W/mx.ref" "$W/mx.kt" 2>/dev/null | cmp -s - "$W/mx.bed" \
  && echo "mixed sentinel blocks, seed $seed: byte-identical" || { echo "mixed sentinel blocks, seed $seed: WRONG"; fail=1; }
done
# On the minus strand a run end maps to a reference run start; the sentinel
# never fired there before.
printf 'chrB\t0\t100\tA\nchrB\t100\t250\tB\nchrB\t250\t400\tC\n' > "$W/mn.ref"
printf 'm\t400\t0\t400\t-\tchrB\t400\t0\t400\t400\t400\t60\tcg:Z:400M\n' > "$W/mn.paf"
printf 'm\t0\t150\tX\nm\t150\t300\tY\nm\t300\t400\tZ\n' > "$W/mn.bed"
$KT encode --aln "$W/mn.paf" --ref "$W/mn.ref" < "$W/mn.bed" > "$W/mn.kt" 2>/dev/null
n=$($KT decode --raw "$W/mn.kt" 2>/dev/null | awk -F'\t' '$2==$3' | wc -l)
$KT decode --aln "$W/mn.paf" --ref "$W/mn.ref" "$W/mn.kt" 2>/dev/null | cmp -s - "$W/mn.bed" && [ "$n" = 3 ] \
  && echo "minus-strand sentinel: 3 of 3 rows, round-trip" || { echo "minus-strand sentinel: $n of 3, WRONG"; fail=1; }

# A NOSZ codec carries no size, so a forged index could not be caught: refused
# from the environment and in a header.
must_fail "--id-codec 0x10, a NOSZ codec" $KT --id-codec 0x10 encode < "$D/s.bed"

# Round 32: flags with no command dereferenced argv[1] past the end of a
# compacted argv and segfaulted (exit 139). The argc<2 guard ran before the
# compaction, so it never saw these.
must_fail "--profile with no command" $KT --profile
must_fail "--threads with no command" $KT --threads 4
must_fail "--tmp-dir with no command" $KT --tmp-dir /tmp
must_fail "--slow with no command" $KT --slow
$KT --profile 2>&1 | grep -q '^karyozip: usage:' \
  && echo "flags with no command: usage, not a crash" \
  || { echo "flags with no command: WRONG"; fail=1; }
python3 "$D/kt.py" "$W/p.kt" "$W/nosz.kt" lo9 0x10
must_fail "a header naming a NOSZ codec" $KT decode "$W/nosz.kt"

#  ranked the plain length spelling in a scratch buffer sized for the
# sentinel one, which is up to 5 bytes a row shorter: with every run ending on
# a reference boundary and lengths near 2^31, the proxy ran out of room.
python3 - "$W" <<'PY'
import random, sys
W = sys.argv[1]; random.seed(9); p = 0; ref = []
for i in range(3000):
    c = p + random.randint(1 << 28, (1 << 31) - 1); ref.append((p, c)); p = c
open(W + "/huge.ref", "w").write("".join("chrH\t%d\t%d\tL%d\n" % (a, b, i % 3) for i, (a, b) in enumerate(ref)))
open(W + "/huge.q", "w").write("".join("h\t%d\t%d\tL%d\n" % (a, b, i % 3) for i, (a, b) in enumerate(ref)))
cg = "".join("%dM" % (b - a) for a, b in ref)   # one op per run: an op stays under 4 GB
open(W + "/huge.paf", "w").write("h\t%d\t0\t%d\t+\tchrH\t%d\t0\t%d\t%d\t%d\t60\tcg:Z:%s\n" % ((p,) * 6 + (cg,)))
PY
$KT encode  --aln "$W/huge.paf" --ref "$W/huge.ref" < "$W/huge.q" > "$W/huge.kt" 2>/dev/null \
  && $KT decode --aln "$W/huge.paf" --ref "$W/huge.ref" "$W/huge.kt" 2>/dev/null | cmp -s - "$W/huge.q" \
  && echo "fast with long plain lengths: byte-identical" || { echo "fast with long plain lengths: FAILED"; fail=1; }
# a codec id with bit 31 set slipped past the NOSZ gate as a negative int
python3 "$D/kt.py" "$W/p.kt" "$W/neg.kt" lo9 0x80010010
must_fail "a codec id with bit 31 set" $KT decode "$W/neg.kt"

# a packed alignment indexes the reference by row order: another order or
# another reference is refused rather than mislabelling every row
$KT aln-pack "$D/p.paf" "$D/ref.bed" "$W/p.kta" >/dev/null 2>&1
tac "$D/ref.bed" | awk -F'\t' '{print}' > "$W/ref.rev"
$KT encode --aln "$W/p.kta" --ref "$D/ref.bed" < "$D/s.bed" > "$W/pk.kt" 2>/dev/null \
  && $KT decode --aln "$W/p.kta" --ref "$D/ref.bed" "$W/pk.kt" 2>/dev/null | cmp -s - "$D/s.bed" \
  && echo "packed alignment: round-trip" || { echo "packed alignment: WRONG"; fail=1; }
must_fail "packed alignment with a reordered reference" $KT decode --aln "$W/p.kta" --ref "$W/ref.rev" "$W/pk.kt"

# a label count with no label table skipped the count check, and smooth then
# wrote a map sized by a wrapped multiply
$KT encode < /dev/null > "$W/mt.kt" 2>/dev/null
python3 "$D/kt.py" "$W/mt.kt" "$W/mt2.kt" 7 0x4000000000000000
must_fail "a label count with no label table" $KT smooth "$W/mt2.kt" "$W/sm.tsv" acrocentric

# Input-side findings: a featureID or set named twice, a name list whose names
# start with '>', a sort temp dir past 511 bytes, a closed stdout, a read error
# on stdin, a name that repeats later in the file, a 512-byte name under --sort.
printf 'featureID\tacrocentric\tchromosome\n1\tA\tX\n1\tB\tY\n' > "$W/dupid.tsv"
printf 's\t0\t10\tA\n' > "$W/d1.bed"; printf 's\t0\t10\tX\n' > "$W/d2.bed"
must_fail "features.tsv with a featureID twice" $KT encode --db "$W/dupid.tsv" "$W/d1.bed" "$W/d2.bed"
printf 'featureID\tx\tx\n1\tA\tX\n' > "$W/dupset.tsv"
must_fail "features.tsv with a set named twice" $KT encode --db "$W/dupset.tsv" "$W/d1.bed" "$W/d2.bed"
printf '>r1 extra\t0\t10\tA\n>r2\t0\t10\tB\n' > "$W/gt.bed"; printf '>r1 extra\n>r2\n' > "$W/gt.names"
$KT encode --names drop < "$W/gt.bed" > "$W/gt.kt" 2>/dev/null
$KT --names "$W/gt.names" decode "$W/gt.kt" 2>/dev/null | cmp -s - "$W/gt.bed" \
  && echo "name list starting with >: restored" || { echo "name list starting with >: WRONG"; fail=1; }
L="$W"; for i in 1 2 3 4 5; do L="$L/$(printf 'd%.0s' $(seq 1 120))"; done; mkdir -p "$L"
echo keep > "${L:0:511}" 2>/dev/null || true
awk 'BEGIN{for(i=0;i<80000;i++)printf "s%d\t0\t10\tA\n",i}' > "$W/long.bed"
$KT encode --tmp-dir "$L" --sort=1M --names drop < "$W/long.bed" > "$W/long.kt" 2>/dev/null \
  && $KT decode --names ordinal "$W/long.kt" 2>/dev/null | cut -f2,3,4 | cmp -s - <(cut -f2,3,4 "$W/long.bed") \
  && [ "$(cat "${L:0:511}" 2>/dev/null)" = keep ] \
  && echo "long sort temp dir: round-trip, bystander intact" || { echo "long sort temp dir: WRONG"; fail=1; }
$KT encode < "$D/s.bed" >&- 2>/dev/null && { echo "encode with stdout closed: EXIT 0"; fail=1; } || echo "encode with stdout closed: refused"
$KT encode < "$W" > "$W/dir.kt" 2>/dev/null && { echo "encode from a directory: EXIT 0"; fail=1; } || echo "encode from a directory: refused"
printf 'a\t0\t10\tA\nb\t0\t10\tB\na\t0\t10\tC\n' > "$W/dupname.bed"
must_fail "a sequence name that repeats later" $KT encode < "$W/dupname.bed"
must_fail "a sequence name that repeats later, sorted" $KT encode --sort < "$W/dupname.bed"
N512=$(head -c 512 /dev/zero | tr '\0' n); printf '%s\t0\t10\tA\n' "$N512" > "$W/n512.bed"
must_fail "a 512-byte name under --sort" $KT encode --sort --names keep < "$W/n512.bed"

# Round 14: two names with one FNV hash, a repeat of the second was missed;
# a spill left behind when the disk filled; CRLF and one-line name lists; a
# features.tsv whose first line is not the header.
printf 'FQ6SH1h576F\t0\t10\tL\nOyHUvXChfZK\t0\t10\tL\nzz\t0\t10\tL\nOyHUvXChfZK\t0\t10\tL\n' > "$W/coll.bed"
must_fail "a repeat behind a hash collision" $KT encode < "$W/coll.bed"
must_fail "a repeat behind a hash collision, sorted" $KT encode --sort < "$W/coll.bed"
printf 'FQ6SH1h576F\t0\t10\tL\nOyHUvXChfZK\t0\t10\tL\n' > "$W/coll2.bed"
$KT encode < "$W/coll2.bed" > "$W/coll2.kt" 2>/dev/null; $KT decode "$W/coll2.kt" 2>/dev/null | cmp -s - "$W/coll2.bed" \
  && echo "two names with one hash: accepted" || { echo "two names with one hash: WRONG"; fail=1; }
mkdir -p "$W/spfull"; (trap "" XFSZ; ulimit -f 32; $KT encode --tmp-dir "$W/spfull" --sort=1M < "$W/long.bed" > /dev/null 2>&1) || true
  [ "$(ls "$W/spfull" | wc -l)" = 0 ] && echo "spill on a full disk: removed" || { echo "spill on a full disk: LEFT"; fail=1; }
printf 'r1\r\nr2\r\n' > "$W/crlf.txt"; printf 'r1\t0\t10\tA\nr2\t0\t10\tB\n' > "$W/n2.bed"
$KT encode --names drop < "$W/n2.bed" > "$W/n2.kt" 2>/dev/null
$KT --names "$W/crlf.txt" decode "$W/n2.kt" 2>/dev/null | cmp -s - "$W/n2.bed" \
  && echo "CRLF name list: restored" || { echo "CRLF name list: WRONG"; fail=1; }
printf '>r1 extra\t0\t10\tA\n' > "$W/gt1.bed"; printf '>r1 extra\n' > "$W/gt1.txt"
$KT encode --names drop < "$W/gt1.bed" > "$W/gt1.kt" 2>/dev/null
$KT --names "$W/gt1.txt" decode "$W/gt1.kt" 2>/dev/null | cmp -s - "$W/gt1.bed" \
  && echo "one-line name list starting with >: restored" || { echo "one-line name list starting with >: WRONG"; fail=1; }
printf 'r1\n' > "$W/dupn.txt"; printf 'r1\nr1\n' > "$W/dupn2.txt"
must_fail "a name file with a repeat, on query" $KT --names "$W/dupn2.txt" region "$W/n2.kt" r1
printf 'featureID acrocentric chromosome\n1\tA\tX\n' > "$W/notab.tsv"; printf 's\t0\t10\t1\n' > "$W/ids1.bed"
must_fail "features.tsv whose first line has no tab" $KT encode --db "$W/notab.tsv" "$W/ids1.bed"

# Round 15: a .idx answered an absent name whose 64-bit hash matched a
# stored one, and one built under a name file answered under another.
printf 'aaa\t0\t10\tA\nYo5P2GphlrB\t0\t10\tB\nYo5P2GphlrB\t10\t30\tC\nzzz\t0\t10\tD\n' > "$W/hc.bed"
$KT encode < "$W/hc.bed" > "$W/hc.kt" 2>/dev/null; $KT index "$W/hc.kt" >/dev/null 2>&1
must_fail "indexed query of a colliding absent name" $KT region "$W/hc.kt" TQyOn5lCDT5
$KT region "$W/hc.kt" Yo5P2GphlrB 2>/dev/null | cmp -s - <(grep -P '^Yo5P2GphlrB\t' "$W/hc.bed") \
  && echo "indexed query of the stored name: correct" || { echo "indexed query of the stored name: WRONG"; fail=1; }
printf 'r0\t0\t10\tA\nr1\t0\t10\tB\nr2\t0\t10\tC\nr3\t0\t10\tD\n' > "$W/nx.bed"
$KT encode --names drop < "$W/nx.bed" > "$W/nx.kt" 2>/dev/null
printf 'r0\nr1\nr2\nr3\n' > "$W/A.txt"; printf 'r3\nr2\nr1\nr0\n' > "$W/B.txt"
$KT --names "$W/A.txt" index "$W/nx.kt" >/dev/null 2>&1
[ "$($KT --names "$W/B.txt" region "$W/nx.kt" r0 2>/dev/null | cut -f4)" = D ] \
  && echo "index under another name file: ignored" || { echo "index under another name file: WRONG"; fail=1; }

# Round 16: a name file rewritten in place, an index entry with its ordinal
# edited, a 5 MB label through smooth, query and decode, and --names ordinal
# on query.
cp "$W/B.txt" "$W/A.txt"
[ "$($KT --names "$W/A.txt" region "$W/nx.kt" r0 2>/dev/null | cut -f4)" = D ] \
  && echo "index after the name file changed: ignored" || { echo "index after the name file changed: WRONG"; fail=1; }
python3 - "$W" <<'PY'
import struct, sys
W = sys.argv[1]; d = bytearray(open(W + "/hc.kt.idx", "rb").read()); n = struct.unpack_from("<Q", d, 8)[0]
o = 32 + 16 * n; ords = list(struct.unpack_from("<%dI" % n, d, o)); ords[0], ords[-1] = ords[-1], ords[0]
struct.pack_into("<%dI" % n, d, o, *ords); open(W + "/hc.kt.idx", "wb").write(d)
PY
$KT region "$W/hc.kt" aaa 2>/dev/null | cmp -s - <(grep -P '^aaa\t' "$W/hc.bed") \
  && echo "index with edited ordinals: not trusted" || { echo "index with edited ordinals: WRONG"; fail=1; }
python3 -c "open('$W/five.bed', 'w').write('s\t0\t10\t%s\n' % ('z' * (5 << 20)))"
$KT encode < "$W/five.bed" > "$W/five.kt" 2>/dev/null
printf 'feature_set\tchild\tparent\nfs\t%s\tcategorized\n' "$(head -c $((5 << 20)) /dev/zero | tr '\0' z)" > "$W/five.tsv"
$KT smooth "$W/five.kt" "$W/five.tsv" fs 2>/dev/null | cmp -s - "$W/five.bed" && $KT decode "$W/five.kt" 2>/dev/null | cmp -s - "$W/five.bed" \
  && $KT region "$W/five.kt" s 2>/dev/null | cmp -s - "$W/five.bed" \
  && echo "5 MB label: smooth, decode, query" || { echo "5 MB label: WRONG"; fail=1; }
$KT encode --names drop < "$D/s.bed" > "$W/nn.kt" 2>/dev/null
$KT region --names ordinal "$W/nn.kt" 0 2>/dev/null | cut -f2,3,4 | cmp -s - <(grep -P '^q1\t' "$D/s.bed" | cut -f2,3,4) \
  && echo "query --names ordinal: correct" || { echo "query --names ordinal: WRONG"; fail=1; }

# Round 17: a 5 MB label with a long name overflowed decode's buffer, and a
# features.tsv label past 1 MB overflowed the per-set writer.
N300=$(head -c 300 /dev/zero | tr '\0' n); python3 -c "open('$W/fivename.bed', 'w').write('$N300\t0\t10\t%s\n' % ('z' * (5 << 20)))"
$KT encode < "$W/fivename.bed" > "$W/fivename.kt" 2>/dev/null
$KT decode "$W/fivename.kt" 2>/dev/null | cmp -s - "$W/fivename.bed" && echo "5 MB label, long name: decode" || { echo "5 MB label, long name: WRONG"; fail=1; }
python3 -c "open('$W/biglab.tsv', 'w').write('featureID\tacrocentric\tchromosome\n1\t%s\tX\n2\tB\tY\n' % ('L' * (2 << 20)))"
printf 's\t0\t10\t1\ns\t10\t20\t2\ns\t20\t30\t1\n' > "$W/bigids.bed"
$KT encode --db "$W/biglab.tsv" "$W/bigids.bed" > "$W/biglab.kt" 2>/dev/null; mkdir -p "$W/bl"
$KT decode --db "$W/biglab.tsv" --out "$W/bl" "$W/biglab.kt" 2>/dev/null && [ "$(wc -l < "$W/bl/acrocentric.bed")" = 3 ] \
  && echo "2 MB features label: decode --db" || { echo "2 MB features label: WRONG"; fail=1; }
must_fail "encode with a reader flag" $KT encode --names ordinal < "$D/s.bed"

# Round 18: decode sized its buffer without the name; a name file can bring
# a 5 MB name, which is refused, or a 1000-byte one, which is refused too.
python3 -c "open('$W/L.bed', 'w').write('0\t0\t10\t' + 'z' * (42 << 17) + '\n')"
$KT encode --names drop < "$W/L.bed" > "$W/L.kt" 2>/dev/null
printf '%s\n' "$(head -c 500 /dev/zero | tr '\0' n)" > "$W/n500.names"
$KT --names "$W/n500.names" decode "$W/L.kt" 2>/dev/null | cut -f1 | cmp -s - <(cut -f1 "$W/n500.names") \
  && echo "4.5 MB label with a 500-byte name: decode" || { echo "4.5 MB label with a 500-byte name: WRONG"; fail=1; }
printf '%s\n' "$(head -c 1000 /dev/zero | tr '\0' n)" > "$W/n1000.names"
must_fail "a name file name past 511 bytes" $KT --names "$W/n1000.names" decode "$W/L.kt"

# On-the-fly: the encoder must take rows as a producer writes them, from a
# pipe and from a slow writer, and give the same bytes as a file encode.
$KT encode < "$D/s.bed" > "$W/file.kt" 2>/dev/null
cat "$D/s.bed" | $KT encode > "$W/pipe.kt" 2>/dev/null
(while IFS= read -r l; do printf '%s\n' "$l"; sleep 0.02; done < "$D/s.bed") | $KT encode > "$W/slow.kt" 2>/dev/null
cmp -s "$W/file.kt" "$W/pipe.kt" && cmp -s "$W/file.kt" "$W/slow.kt" \
  && echo "encode from a pipe and a slow producer: byte-identical" || { echo "encode from a pipe: WRONG"; fail=1; }
awk 'BEGIN{srand(3);for(i=0;i<20000;i++){p=0;for(j=0;j<5;j++){c=p+int(rand()*300)+1;printf "read%d\t%d\t%d\tL%d\n",i,p,c,int(rand()*6);p=c}}}' > "$W/pbig.bed"
cat "$W/pbig.bed" | $KT encode --tmp-dir "$W" --sort=1M > "$W/pipesort.kt" 2>/dev/null
$KT encode --tmp-dir "$W" --sort=1M < "$W/pbig.bed" > "$W/filesort.kt" 2>/dev/null
cmp -s "$W/pipesort.kt" "$W/filesort.kt" && echo "sorted encode from a pipe: byte-identical" || { echo "sorted encode from a pipe: WRONG"; fail=1; }
$KT decode --raw "$W/f.kt" 2>/dev/null | $KT encode --db "$D/db.features.tsv" - > "$W/pipedb.kt" 2>/dev/null
$KT decode --raw "$W/f.kt" 2>/dev/null > "$W/merged.bed"; $KT encode --db "$D/db.features.tsv" "$W/merged.bed" > "$W/filedb.kt" 2>/dev/null
cmp -s "$W/pipedb.kt" "$W/filedb.kt" && echo "db encode from a pipe: byte-identical" || { echo "db encode from a pipe: WRONG"; fail=1; }

# Round 19: an all-gap PAF record desynchronised the packed alignment; an
# unsorted reference made encoder and decoder disagree; nothing bound a track
# to its PAF, reference and name map; a strand other than + or -.
printf 'chrA\t0\t100\tA\nchrA\t100\t200\tB\nchrB\t0\t100\tC\nchrB\t100\t200\tD\n' > "$W/r19.ref"
printf 'q0\t0\t0\t0\t+\tchrA\t200\t0\t100\t0\t100\t60\tcg:Z:100D\nq1\t200\t0\t200\t-\tchrB\t200\t0\t200\t200\t200\t60\tcg:Z:200M\n' > "$W/r19.paf"
printf 'q1\t0\t100\tD\nq1\t100\t200\tC\n' > "$W/r19.bed"
$KT encode --aln "$W/r19.paf" --ref "$W/r19.ref" < "$W/r19.bed" > "$W/r19p.kt" 2>/dev/null
$KT aln-pack "$W/r19.paf" "$W/r19.ref" "$W/r19.kta" >/dev/null 2>&1
$KT encode --aln "$W/r19.kta" --ref "$W/r19.ref" < "$W/r19.bed" > "$W/r19k.kt" 2>/dev/null
cmp -s "$W/r19p.kt" "$W/r19k.kt" && $KT decode --aln "$W/r19.kta" --ref "$W/r19.ref" "$W/r19p.kt" 2>/dev/null | cmp -s - "$W/r19.bed" \
  && echo "all-gap PAF record with a packed alignment: byte-identical" || { echo "all-gap PAF record: WRONG"; fail=1; }
printf 'chrA\t100\t200\tB\nchrA\t0\t100\tA\n' > "$W/unsorted.ref"
must_fail "an unsorted reference" $KT encode --aln "$W/r19.paf" --ref "$W/unsorted.ref" < "$W/r19.bed"
sed 's/\tD$/\tE/' "$W/r19.ref" > "$W/other.ref"
must_fail "decode with another reference" $KT decode --aln "$W/r19.paf" --ref "$W/other.ref" "$W/r19p.kt"
sed 's/200M/100M100M/' "$W/r19.paf" > "$W/other.paf"
must_fail "decode with another PAF" $KT decode --aln "$W/other.paf" --ref "$W/r19.ref" "$W/r19p.kt"
sed 's/\t-\t/\t*\t/' "$W/r19.paf" > "$W/star.paf"
must_fail "a PAF strand of *" $KT encode --aln "$W/star.paf" --ref "$W/r19.ref" < "$W/r19.bed"

# Round 20: an empty reference contig name crashed; a 300-byte contig name
# was refused; a PAF in another record order was refused at decode; aln-pack
# took a flag as its output name.
printf '\t0\t100\tA\n' > "$W/empty.ref"
must_fail "reference with an empty contig name" $KT encode --aln "$D/p.paf" --ref "$W/empty.ref" < "$D/s.bed"
N300=$(head -c 300 /dev/zero | tr '\0' N)
printf "$N300\t0\t100\tA\n$N300\t100\t200\tB\n" > "$W/n300.ref"
printf "x\t200\t0\t200\t+\t$N300\t200\t0\t200\t200\t200\t60\tcg:Z:200M\n" > "$W/n300.paf"
printf 'x\t0\t100\tA\nx\t100\t200\tB\n' > "$W/n300.bed"
$KT encode --aln "$W/n300.paf" --ref "$W/n300.ref" < "$W/n300.bed" > "$W/n300.kt" 2>/dev/null \
  && $KT decode --aln "$W/n300.paf" --ref "$W/n300.ref" "$W/n300.kt" 2>/dev/null | cmp -s - "$W/n300.bed" \
  && echo "300-byte reference contig name: round-trip" || { echo "300-byte reference contig name: WRONG"; fail=1; }
tac "$D/p.paf" > "$W/rev.paf"
$KT decode --aln "$W/rev.paf" --ref "$D/ref.bed" "$W/d.kt" 2>/dev/null | cmp -s - "$D/s.bed" \
  && echo "PAF in another record order: accepted" || { echo "PAF in another record order: WRONG"; fail=1; }
must_fail "aln-pack with a flag as output" $KT aln-pack "$D/p.paf" "$D/ref.bed" --names "$W/x.kta"

# Round 21: a CRLF hierarchy, a four-column row, a child listed twice, a
# root other than categorized, smooth region arguments, a halo that wraps,
# and the V1 suffix on a long label.
printf 'feature_set\tchild\tparent\r\nacrocentric\tacro\tcategorized\r\nacrocentric\tnovel\tcategorized\r\n' > "$W/crlf.tsv"
$KT smooth "$W/sm.kt" "$W/crlf.tsv" acrocentric 2>/dev/null | cmp -s - "$W/sm.1" \
  && echo "CRLF hierarchy: same rows" || { echo "CRLF hierarchy: WRONG"; fail=1; }
printf 'feature_set\tchild\tparent\nacrocentric\tacro\tcategorized\tx\n' > "$W/four.tsv"
must_fail "hierarchy row with four columns" $KT smooth "$W/sm.kt" "$W/four.tsv" acrocentric
printf 'feature_set\tchild\tparent\nacrocentric\tacro\tcategorized\nacrocentric\tacro\tnovel\n' > "$W/dupc.tsv"
must_fail "hierarchy child listed twice" $KT smooth "$W/sm.kt" "$W/dupc.tsv" acrocentric
printf 'feature_set\tchild\tparent\nacrocentric\tacro\troot\n' > "$W/root.tsv"
must_fail "hierarchy root that is not categorized" $KT smooth "$W/sm.kt" "$W/root.tsv" acrocentric
must_fail "smooth region with a junk start" $KT smooth "$W/sm.kt" "$W/sm.tsv" acrocentric r0 1e2 3e2
must_fail "smooth region with end before start" $KT smooth "$W/sm.kt" "$W/sm.tsv" acrocentric r0 300 100
$KT region "$W/halo1.kt" s 45 46 18446744073709551614 2>/dev/null | cmp -s - <($KT region "$W/halo1.kt" s 45 46 100 2>/dev/null) \
  && echo "halo of 2^64-2: whole sequence" || { echo "halo of 2^64-2: WRONG"; fail=1; }
L309=$(head -c 300 /dev/zero | tr '\0' L)_specific
printf 's\t0\t10\t%s\n' "$L309" > "$W/l309.bed"; $KT encode < "$W/l309.bed" > "$W/l309.kt" 2>/dev/null
printf 'feature_set\tchild\tparent\nfs\t%s\tcategorized\n' "$L309" > "$W/l309.tsv"
[ "$($KT --v1 smooth "$W/l309.kt" "$W/l309.tsv" fs 2>/dev/null | cut -f4 | wc -c)" = 301 ] \
  && echo "V1 strip on a 309-byte label: stripped" || { echo "V1 strip on a 309-byte label: WRONG"; fail=1; }

# Round 22: a refused decode --db truncated the old per-set BEDs; touching
# same-label runs in one set were merged on decode; "+1" was a featureID;
# 33 sets with long labels and a 300-byte label were refused; a last line
# without a newline was taken.
mkdir -p "$W/keepdir"; printf 'PRECIOUS\n' > "$W/keepdir/acrocentric.bed"
$KT decode --db "$D/db.features.tsv" --out "$W/keepdir" "$W/p.kt" >/dev/null 2>&1 && { echo "decode --db of a plain track: EXIT 0"; fail=1; }
[ "$(cat "$W/keepdir/acrocentric.bed")" = PRECIOUS ] && echo "refused decode --db: old BEDs intact" || { echo "refused decode --db: OLD BEDS DESTROYED"; fail=1; }
printf 'featureID\ts1\ts2\n1\tx\tp\n2\ty\tp\n' > "$W/two.tsv"
printf 's\t0\t10\tx\ns\t10\t20\ty\n' > "$W/n1.s1.bed"; printf 's\t0\t10\tp\ns\t10\t20\tp\n' > "$W/n1.s2.bed"
must_fail "touching same-label runs in one set" $KT encode --db "$W/two.tsv" "$W/n1.s1.bed" "$W/n1.s2.bed"
printf 's\t0\t10\t+1\ns\t10\t20\t2\n' > "$W/plus.bed"
must_fail "featureID +1 in a merged track" $KT encode --db "$W/two.tsv" "$W/plus.bed"
L201=$(head -c 201 /dev/zero | tr '\0' L); { printf 'featureID'; for i in $(seq 33); do printf '\ts%d' $i; done; printf '\n1'; for i in $(seq 33); do printf '\t%s%d' "$L201" $i; done; printf '\n'; } > "$W/t33.tsv"
f33=""; for i in $(seq 33); do printf 's\t0\t10\t%s%d\n' "$L201" $i > "$W/t33.$i.bed"; f33="$f33 $W/t33.$i.bed"; done
$KT encode --db "$W/t33.tsv" $f33 > "$W/t33.kt" 2>/dev/null && echo "33 sets with 201-byte labels: encoded" || { echo "33 sets with 201-byte labels: REFUSED"; fail=1; }
L300=$(head -c 300 /dev/zero | tr '\0' L); printf "featureID\ts1\ts2\n1\t$L300\tp\n" > "$W/l300.tsv"
printf "s\t0\t10\t$L300\n" > "$W/l300.s1.bed"; printf 's\t0\t10\tp\n' > "$W/l300.s2.bed"; mkdir -p "$W/l300o"
$KT encode --db "$W/l300.tsv" "$W/l300.s1.bed" "$W/l300.s2.bed" > "$W/l300.kt" 2>/dev/null \
  && $KT decode --db "$W/l300.tsv" --out "$W/l300o" "$W/l300.kt" 2>/dev/null && cmp -s "$W/l300o/s1.bed" "$W/l300.s1.bed" \
  && echo "300-byte label in a set: round-trip" || { echo "300-byte label in a set: WRONG"; fail=1; }
head -c -1 "$D/db.gene.bed" > "$W/nonl.gene.bed"
must_fail "a set BED without a final newline" $KT encode --db "$D/db.features.tsv" "$D/db.acrocentric.bed" "$D/db.chromosome.bed" "$W/nonl.gene.bed" "$D/db.region.bed" "$D/db.repeat.bed" "$D/db.subtelomeric.bed"

# Round 23: a first sequence named \x01 matched the "no sequence yet"
# sentinel of --db; a long --out path was truncated silently; MAXGAP near 2^64
# wrapped; region with END before START printed nothing.
printf '\001\t0\t10\t1\n\001\t10\t20\t2\nb\t0\t10\t1\n' > "$W/c3.bed"
$KT encode --db "$W/two.tsv" "$W/c3.bed" > "$W/c3.kt" 2>/dev/null && $KT decode --raw "$W/c3.kt" 2>/dev/null | cmp -s - "$W/c3.bed" \
  && echo "first sequence named 0x01 under --db: round-trip" || { echo "first sequence named 0x01 under --db: WRONG"; fail=1; }
LD="$W"; while [ ${#LD} -lt 4080 ]; do LD="$LD/$(printf 'd%.0s' $(seq 1 100))"; done; mkdir -p "$LD" 2>/dev/null
must_fail "decode --db --out with a 4080-byte path" $KT decode --db "$D/db.features.tsv" --out "$LD" "$W/f.kt"
$KT --max-gap 18446744073709551615 smooth "$W/sm.kt" "$W/sm.tsv" acrocentric 2>/dev/null | cmp -s - "$W/sm.1" \
  && echo "MAXGAP of 2^64-1: same as default" || { echo "MAXGAP of 2^64-1: WRONG"; fail=1; }
must_fail "region with END before START" $KT region "$W/p.kt" q1 150 20

# Round 24: 100,000 threads crashed in libgomp; a file-size limit killed
# encode by signal; --slow used zstd-22 over an explicit --level; the block
# temp file ignored TMPDIR.
must_fail "--threads 100000" $KT -t 100000 encode < "$D/s.bed"
rc=0; (ulimit -f 1; $KT encode < "$W/pbig.bed" > "$W/lim.kt" 2> "$W/lim.err") || rc=$?
[ $rc = 1 ] && echo "file-size limit: exit 1, not a signal" || { echo "file-size limit: rc=$rc"; fail=1; }
$KT --len-codec 0 --id-codec 0 encode --level 1 --block all < "$W/pbig.bed" > "$W/l5.kt" 2>/dev/null
$KT --len-codec 0 --id-codec 0 encode --level 1 --block all --slow < "$W/pbig.bed" > "$W/l5s.kt" 2>/dev/null
[ "$(stat -c%s "$W/l5.kt")" = "$(stat -c%s "$W/l5s.kt")" ] && echo "--slow keeps --level: same size" || { echo "--slow keeps --level: WRONG"; fail=1; }
must_fail "--tmp-dir that does not exist" $KT --tmp-dir /nonexistent/dir encode < "$D/s.bed"
$KT --tmp-dir "$W" encode < "$D/s.bed" 2>/dev/null | cmp -s - "$W/p.kt" && echo "--tmp-dir honoured: same bytes" || { echo "--tmp-dir honoured: WRONG"; fail=1; }

# Round 25: OMP_NUM_THREADS bypassed the thread cap; a stream of exactly
# 2^31-1 bytes under rANS wrote an unreadable file; arith X_EXT double-freed;
# numeric options took junk suffixes.
OMP_NUM_THREADS=100000 $KT encode < "$D/s.bed" 2>/dev/null | cmp -s - "$W/p.kt" && echo "OMP_NUM_THREADS 100000: capped" || { echo "OMP_NUM_THREADS 100000: WRONG"; fail=1; }
must_fail "--id-codec with arith X_EXT" $KT --id-codec 0x10004 encode < "$D/s.bed"
must_fail "--level abc" $KT encode --level abc < "$D/s.bed"
must_fail "--block 4abc" $KT encode --block 4abc < "$D/s.bed"
must_fail "-t 4x" $KT -t 4x encode < "$D/s.bed"
must_fail "--sort=nan" $KT encode --tmp-dir "$W" --sort=nan < "$D/s.bed"

# Round 26: an index built while its name file was replaced answered from the
# new file; a failed decode --db destroyed the previous output; a signal left
# spill files; junk after the unit.
cp "$W/A.txt" "$W/A2.txt"; $KT --names "$W/A2.txt" index "$W/nx.kt" >/dev/null 2>&1; cp "$W/B.txt" "$W/A2.tmp"; mv -f "$W/A2.tmp" "$W/A2.txt"
[ "$($KT --names "$W/A2.txt" region "$W/nx.kt" r0 2>/dev/null | cut -f4)" = D ] \
  && echo "index built before a name file swap: ignored" || { echo "index built before a name file swap: WRONG"; fail=1; }
mkdir -p "$W/good"; $KT decode --db "$D/db.features.tsv" --out "$W/good" "$W/f.kt" >/dev/null 2>&1; m1=$(cat "$W/good"/*.bed | md5sum)
cp "$W/f.kt" "$W/fbad.kt"; sz=$(stat -c%s "$W/fbad.kt"); printf '\xff' | dd of="$W/fbad.kt" bs=1 seek=$((sz-100)) conv=notrunc status=none
$KT decode --db "$D/db.features.tsv" --out "$W/good" "$W/fbad.kt" >/dev/null 2>&1 && { echo "decode of a damaged db track: EXIT 0"; fail=1; }
[ "$(cat "$W/good"/*.bed | md5sum)" = "$m1" ] && [ -z "$(ls "$W/good" | grep tmp)" ] \
  && echo "failed decode --db: previous output intact, no temporaries" || { echo "failed decode --db: OUTPUT DAMAGED"; fail=1; }
mkdir -p "$W/sigspill"; ( (awk '{print; if (NR % 2000 == 0) system("sleep 0.05")}' "$W/long.bed"; sleep 3) | $KT encode --tmp-dir "$W/sigspill" --sort=1M > /dev/null 2>&1 & p=$!; sleep 1.5; kill -INT $p 2>/dev/null; wait $p ) 2>/dev/null || true
[ "$(ls "$W/sigspill" | wc -l)" = 0 ] && echo "SIGINT during a spilling sort: no files left" || { echo "SIGINT during a spilling sort: FILES LEFT"; fail=1; }
must_fail "--sort=1M2" $KT encode --tmp-dir "$W" --sort=1M2 < "$D/s.bed"

# Round 27: a name file on a pipe lost its first lines to a failed rewind; a
# failed rename left the db output half-updated; a signal left db temporaries;
# its cast.
printf 'r0\nr1\nr2\nr3\n' > "$W/pn.txt"
$KT decode --names ordinal "$W/nx.kt" 2>/dev/null | cut -f2,3,4 > "$W/pn.want"
$KT --names <(cat "$W/pn.txt") decode "$W/nx.kt" 2>/dev/null | cmp -s - "$W/nx.bed" \
  && echo "name file from a pipe: restored" || { echo "name file from a pipe: WRONG"; fail=1; }
printf '@r0\nACGT\n+\nIIII\n@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIIII\n@r3\nACGT\n+\nIIII\n' > "$W/pn.fq"
$KT --names <(cat "$W/pn.fq") decode "$W/nx.kt" 2>/dev/null | cmp -s - "$W/nx.bed" \
  && echo "FASTQ from a pipe: restored" || { echo "FASTQ from a pipe: WRONG"; fail=1; }
mkdir -p "$W/mix"; for f in acrocentric chromosome gene region repeat subtelomeric; do echo "OLD_$f" > "$W/mix/$f.bed"; done
mkdir -p "$W/mix/gene.bed.blocked" && mv "$W/mix/gene.bed" "$W/mix/gene.keep" && mkdir "$W/mix/gene.bed"
must_fail "decode --db over a directory in --out" $KT decode --db "$D/db.features.tsv" --out "$W/mix" "$W/f.kt"
[ "$(cat "$W/mix/acrocentric.bed")" = OLD_acrocentric ] && [ -z "$(ls "$W/mix" | grep -E 'tmp|old')" ] \
  && echo "refused db decode: all old files kept, no leftovers" || { echo "refused db decode: HALF UPDATED"; fail=1; }
mkdir -p "$W/dbsig"; ( $KT decode --db "$D/db.features.tsv" --out "$W/dbsig" "$W/f.kt" >/dev/null 2>&1 & p=$!; sleep 0.02; kill -INT $p 2>/dev/null; wait $p ) 2>/dev/null || true
[ -z "$(ls "$W/dbsig" | grep tmp)" ] && echo "SIGINT during decode --db: no temporaries left" || { echo "SIGINT during decode --db: TEMPORARIES LEFT"; fail=1; }
must_fail "--sort=1e18T" $KT encode --tmp-dir "$W" --sort=1e18T < "$D/s.bed"

# The README says these pairs are exclusive and these ranges are the limits.
must_fail "--names drop with --names keep" $KT encode --names drop --names keep < "$D/s.bed"
must_fail "--names keep with --names drop" $KT encode --names keep --names drop < "$D/s.bed"
must_fail "the retired --one-block" $KT encode --one-block < "$D/s.bed"
must_fail "--block junk" $KT encode --block xyz < "$D/s.bed"
must_fail "--threads 1025" $KT -t 1025 encode < "$D/s.bed"
must_fail "--level 23" $KT encode --level 23 < "$D/s.bed"
must_fail "--level -1" $KT encode --level -1 < "$D/s.bed"
$KT -t 1024 encode --level 22 < "$D/s.bed" > /dev/null 2>&1 && echo "the documented limits: accepted" || { echo "the documented limits: REFUSED"; fail=1; }

# Every flag is parsed once, before the command, so it works on either side of
# it and a flag on the wrong command is refused instead of read as an argument.
$KT --slow encode < "$D/s.bed" > "$W/side1.kt" 2>/dev/null
$KT encode --slow < "$D/s.bed" > "$W/side2.kt" 2>/dev/null
cmp -s "$W/side1.kt" "$W/side2.kt" && echo "a flag before the command: same bytes" \
  || { echo "a flag before the command: DIFFERENT"; fail=1; }
# A names file that is not the one drop:FILE wrote is refused by name, not as
# a corrupt track.
$KT encode --names "drop:$W/side.names" < "$D/s.bed" > "$W/side.kt" 2>/dev/null
awk '{print $1}' "$D/s.bed" | uniq | sed 's/^/>/' > "$W/side.fa"
$KT decode --names "$W/side.fa" "$W/side.kt" 2>&1 | grep -q "does not belong to this track" \
  && echo "a foreign names file: named as such" || { echo "a foreign names file: WRONG MESSAGE"; fail=1; }
$KT decode --names "$W/side.names" "$W/side.kt" 2>/dev/null | cmp -s - "$D/s.bed" \
  && echo "the sidecar it wrote: round trips" || { echo "the sidecar it wrote: BROKEN"; fail=1; }

must_fail "--len-codec 999999" $KT encode --len-codec 999999 < "$D/s.bed"
must_fail "--id-codec 999999" $KT encode --id-codec 999999 < "$D/s.bed"
# A sorted track decodes in sorted order, never the input order.
$KT encode --sort < "$D/s.bed" > "$W/so.kt" 2>/dev/null
$KT decode "$W/so.kt" 2>/dev/null | cmp -s - "$D/s.bed" \
  && { echo "sorted decode: CLAIMS INPUT ORDER"; fail=1; } || echo "sorted decode: sorted order, as documented"
$KT decode "$W/so.kt" 2>/dev/null | sort | cmp -s - <(sort "$D/s.bed") \
  && echo "sorted decode: same rows" || { echo "sorted decode: ROWS DIFFER"; fail=1; }

must_fail "--preset, a flag that was removed" $KT encode --preset illumina
# A file that is not KST is named as such, and a future version is told apart
# from a foreign file.
printf 'KTRACK02rest' > "$W/old.kst"
$KT decode "$W/old.kst" 2>&1 | grep -q 'not "KST00001"' \
  && echo "a foreign magic: named" || { echo "a foreign magic: WRONG MESSAGE"; fail=1; }
printf 'KST00099rest' > "$W/v99.kst"
$KT decode "$W/v99.kst" 2>&1 | grep -q "KST version 00099" \
  && echo "a future KST version: named" || { echo "a future KST version: WRONG MESSAGE"; fail=1; }

must_fail "unknown flag" $KT --bogus encode
$KT --help > "$W/help.txt" 2>/dev/null && [ -s "$W/help.txt" ] \
  && echo "--help: on stdout, exit 0" || { echo "--help: WRONG"; fail=1; }
must_fail "query --db" $KT region "$W/side1.kt" --db x
must_fail "encode --out" $KT encode --out "$W" < "$D/s.bed"
must_fail "encode --names ordinal" $KT encode --names ordinal < "$D/s.bed"
must_fail "decode --tail" $KT decode --tail 1 "$W/side1.kt"
must_fail "decode two tracks" $KT decode "$W/side1.kt" "$W/side2.kt"
must_fail "the retired --sort-mem" $KT encode --sort-mem 1M

# htscodecs wrote "compiled without libbz2" for 32 candidates that were never
# meant to reach arith: nothing from htscodecs may reach stderr.
$KT encode < "$D/s.bed" 2>&1 >/dev/null | grep -qi "htscodecs" \
  && { echo "encode stderr: HTSCODECS COMPLAINED"; fail=1; } || echo "encode stderr: clean"

# Generative stage. The checks above are one fixed case per bug already found;
# this one searches input space instead of remembering it.
if bash "$D/fuzz.sh" "$KT"; then echo "fuzz: clean"; else echo "fuzz: FOUND SOMETHING"; fail=1; fi

exit $fail
