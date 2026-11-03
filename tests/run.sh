#!/bin/bash
# aln.expected cases, all computed by hand:
#   q1  forward, exact
#   q2  minus strand, maps to tp+n-1-pos
#   q3  50M20D50M, the deletion shifts the target
#   q4  50M20I50M, insertion bases have no target
#   q5  no alignment at all
#   run.sh BINARY LIBTEST [SANITIZER_BINARY]
# LIBTEST: C library test, compared to BINARY
set -eu
KT=$(realpath -m "$1"); LT=$(realpath -m "$2"); SAN=$(realpath -m "${3:-$1}")
D=$(dirname "$0")
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
src=$D; D=$W/fx; mkdir -p "$D"; cp "$src"/fixtures.txt "$src"/*.py "$src"/*.sh "$D"/   # the fixtures, unpacked
awk '/^@@ /{f=d"/"$2; next} {print > f}' d="$D" "$D/fixtures.txt"
fail=0
ok() { echo "$1"; }
no() { echo "$1"; fail=1; }
must_fail() {  # DESC CMD...: expect exit 1, message
  local what=$1 rc=0; shift
  "$@" >/dev/null 2>"$W/mf.err" || rc=$?
  if [ $rc = 1 ] && grep -q '^karyozip: ' "$W/mf.err"; then echo "$what: refused ($(grep -m1 '^karyozip: ' "$W/mf.err" | cut -c11-90))"
  else echo "$what: WRONG, exit $rc: $(head -c 200 "$W/mf.err")"; fail=1; fi
}

$KT compress --aln "$D/p.paf" --ref "$D/ref.bed" < "$D/s.bed" > "$W/d.kt" 2>/dev/null
$KT decompress --raw -i "$W/d.kt" 2>/dev/null | cut -f2,3,4 > "$W/got"
cut -f2,3,4 "$D/aln.expected" > "$W/want"
cmp -s "$W/got" "$W/want" && echo "aln delta: correct" || { echo "aln delta: WRONG"; diff "$W/got" "$W/want" | head; fail=1; }

$KT decompress --aln "$D/p.paf" --ref "$D/ref.bed" -i "$W/d.kt" 2>/dev/null > "$W/rt"
cmp -s "$W/rt" "$D/s.bed" && ok "aln restore: byte-identical" || no "aln restore: WRONG"
# a .kst reference: same bytes, same rows
$KT compress < "$D/ref.bed" > "$W/ref.kst" 2>/dev/null
$KT compress --aln "$D/p.paf" --ref "$W/ref.kst" < "$D/s.bed" 2>/dev/null | cmp -s - "$W/d.kt" && ok "aln with a .kst reference: same bytes" || no "aln with a .kst reference: WRONG"
$KT decompress --aln "$D/p.paf" --ref "$W/ref.kst" -i "$W/d.kt" 2>/dev/null | cmp -s - "$D/s.bed" && ok "aln restore from a .kst reference: byte-identical" || no "aln restore from a .kst reference: WRONG"
# .kst reference across many contigs and blocks
$KT compress --block 4 < "$D/db.rall.bed" > "$W/rall.kst" 2>/dev/null
$KT compress --block 4 --aln "$D/db.self.paf" --ref "$D/db.rall.bed" < "$D/db.rall.bed" > "$W/rall.aln.kst" 2>/dev/null
$KT decompress --aln "$D/db.self.paf" --ref "$D/db.rall.bed" -i "$W/rall.aln.kst" > "$W/rall.want" 2>/dev/null
$KT decompress --aln "$D/db.self.paf" --ref "$W/rall.kst" -i "$W/rall.aln.kst" 2>/dev/null | cmp -s - "$W/rall.want" \
  && ok "a .kst reference read across contigs: same rows" || no "a .kst reference read across contigs: WRONG"
must_fail "an --aln track as reference" $KT decompress --aln "$D/p.paf" --ref "$W/d.kt" -i "$W/d.kt"
# eight readers at once share no state
$KT decompress -i "$W/rall.kst" > "$W/conc.want" 2>/dev/null
seq=$(head -1 "$D/db.rall.bed" | cut -f1)
$KT query -i "$W/rall.kst" -r "$seq" > "$W/conc.rwant" 2>/dev/null
for i in 1 2 3 4 5 6 7 8; do
  ( $KT decompress -i "$W/rall.kst" > "$W/c.d.$i" 2>/dev/null
    $KT query -i "$W/rall.kst" -r "$seq" > "$W/c.r.$i" 2>/dev/null
    $KT decompress --aln "$D/db.self.paf" --ref "$W/rall.kst" -i "$W/rall.aln.kst" > "$W/c.a.$i" 2>/dev/null ) &
done
wait
bad=0
for i in 1 2 3 4 5 6 7 8; do
  cmp -s "$W/c.d.$i" "$W/conc.want" || bad=1
  cmp -s "$W/c.r.$i" "$W/conc.rwant" || bad=1
  cmp -s "$W/c.a.$i" "$W/rall.want" || bad=1
done
[ $bad = 0 ] && ok "eight readers at once: same bytes" || no "eight readers at once: WRONG"
# review cases of 2026-09-07, one each
python3 - "$W" <<'PY'
import sys; W = sys.argv[1]; L = "L" * 2100000
open(W + "/rv.tsv", "w").write("featureID\tA\tB\n1\tx\ty\n2\t%s\tz\n" % L)
open(W + "/rv1.bed", "w").write("s1\t0\t10\tx\ns1\t10\t20\t%s\n" % L); open(W + "/rv2.bed", "w").write("s1\t0\t10\ty\ns1\t10\t20\tz\n")
open(W + "/rv.big", "w").write("1125899906842625\t0\t10\tx\n"); open(W + "/rv.names", "w").write("q1\n\nq2\n")
PY
$KT --db "$W/rv.tsv" compress -i "$W/rv1.bed" -i "$W/rv2.bed" > "$W/rv.kt" 2>/dev/null
[ "$($KT --db "$W/rv.tsv" query -i "$W/rv.kt" -r s1:12-13 2>/dev/null | wc -c)" = 2100012 ] && ok "decompress --db, a 2100000-byte label: printed" || no "decompress --db, a 2100000-byte label: WRONG"
$KT compress < "$W/rv.big" > "$W/rv.big.kt" 2>/dev/null; $KT decompress -i "$W/rv.big.kt" 2>/dev/null | cmp -s - "$W/rv.big" && ok "a name past 2^50: round trip" || no "a name past 2^50: WRONG"
must_fail "a blank --aln-names line" $KT compress --aln "$D/p.paf" --ref "$D/ref.bed" --aln-names "$W/rv.names" < "$D/s.bed"
gzip -c "$D/ref.bed" | head -c 40 > "$W/rv.gz"
must_fail "a cut .gz reference" $KT compress --aln "$D/p.paf" --ref "$W/rv.gz" < "$D/s.bed"
$KT aln-pack -a "$D/p.paf" --ref "$D/ref.bed" -o "$W/rv.pk" >/dev/null 2>&1; printf '\0' >> "$W/rv.pk"
must_fail "a packed alignment with a byte over" $KT decompress --aln "$W/rv.pk" --ref "$D/ref.bed" -i "$W/d.kt"

$KT compress < "$D/s.bed" > "$W/p.kt" 2>/dev/null
$KT compress < "$D/db.rall.bed" > "$W/omp.ref" 2>/dev/null
$KT decompress -i "$W/p.kt" 2>/dev/null | cmp -s - "$D/s.bed" && ok "plain: byte-identical" || no "plain: WRONG"

# dropped names need --names FILE
cut -f1 "$D/s.bed" | uniq > "$W/s.names"
for f in "--names drop" "--block 262144" "--block 4"; do
  ofl=""; [ "$f" = "--names drop" ] && ofl="--names $W/s.names"
  $KT compress $f < "$D/s.bed" > "$W/f.kt" 2>/dev/null || { echo "$f: ENCODE FAILED"; fail=1; continue; }
  $KT decompress $ofl -i "$W/f.kt" 2>/dev/null | cmp -s - "$D/s.bed" && ok "$f: byte-identical" || no "$f: WRONG"
done
# --db: 6 BEDs in, 6 out, byte-compared
sets="acrocentric chromosome gene region repeat subtelomeric"
in=""; for f in $sets; do in="$in -i $D/db.$f.bed"; done
$KT compress --db "$D/db.features.tsv" $in > "$W/f.kt" 2>/dev/null
mkdir -p "$W/o"
$KT decompress --db "$D/db.features.tsv" --out-dir "$W/o" -i "$W/f.kt" 2>/dev/null
n=0; for f in $sets; do cmp -s "$W/o/$f.bed" "$D/db.$f.bed" && n=$((n+1)); done
[ "$n" = 6 ] && ok "db: 6/6 identical" || no "db: only $n/6 identical"

# --db with mixed sequence name lengths
in2=""; for f in $sets; do in2="$in2 -i $D/db2.$f.bed"; done
$KT compress --db "$D/db2.features.tsv" $in2 > "$W/f2.kt" 2>/dev/null
mkdir -p "$W/o2"
$KT decompress --db "$D/db2.features.tsv" --out-dir "$W/o2" -i "$W/f2.kt" 2>/dev/null
n=0; for f in $sets; do cmp -s "$W/o2/$f.bed" "$D/db2.$f.bed" && n=$((n+1)); done
[ "$n" = 6 ] && ok "db mixed name lengths: 6/6 identical" || no "db mixed name lengths: only $n/6"

# a transform track refuses decode without it
must_fail "decompress --aln track without --aln" $KT decompress -i "$W/d.kt"
must_fail "decompress plain track with --aln" $KT decompress --aln "$D/p.paf" --ref "$D/ref.bed" -i "$W/p.kt"
must_fail "decompress --db track without --db" $KT decompress -i "$W/f.kt"
must_fail "compress --ref without --aln" $KT compress --ref "$D/ref.bed"
must_fail "--aln-names without --aln" $KT compress --aln-names "$D/ord.names.txt" < "$D/s.bed"
must_fail "a name file given to compress as --names" $KT compress --names "$D/ord.names.txt" < "$D/s.bed"

# --aln needs names to compute deltas
must_fail "compress --aln with --names drop and no name list" \
  $KT compress --aln "$D/p.paf" --ref "$D/ref.bed" --names drop

# --db with --aln: check ids, not labels
$KT compress --db "$D/db.features.tsv" --aln "$D/db.self.paf" --ref "$D/db.rall.bed" $in \
  > "$W/da.kt" 2>/dev/null
mkdir -p "$W/dao"
$KT decompress --db "$D/db.features.tsv" --aln "$D/db.self.paf" --ref "$D/db.rall.bed" --out-dir "$W/dao" -i "$W/da.kt" 2>/dev/null
n=0; for f in $sets; do cmp -s "$W/dao/$f.bed" "$D/db.$f.bed" && n=$((n+1)); done
[ "$n" = 6 ] && ok "db+aln: 6/6 identical" || no "db+aln: only $n/6"
# smooth and region on every informed form
{ printf 'feature_set\tchild\tparent\n'; for l in $(cut -f4 "$D/db.acrocentric.bed" "$D/s.bed" "$D/ref.bed" | sort -u); do printf 'acrocentric\t%s\tcategorized\n' "$l"; done; } > "$W/inf.tsv"
$KT compress < "$D/db.acrocentric.bed" > "$W/pa.kt" 2>/dev/null
$KT decompress -i "$W/pa.kt" --hier "$W/inf.tsv" --smooth --set acrocentric > "$W/inf.want" 2>/dev/null
$KT decompress -i "$W/p.kt" --hier "$W/inf.tsv" --smooth --set acrocentric > "$W/inf.wants" 2>/dev/null
$KT --db "$D/db.features.tsv" decompress -i "$W/f.kt" --hier "$W/inf.tsv" --smooth --set acrocentric 2>/dev/null | cmp -s - "$W/inf.want" && ok "decompress --db: same as plain" || no "decompress --db: DIFFERS"
$KT --db "$D/db.features.tsv" --aln "$D/db.self.paf" --ref "$D/db.rall.bed" -t 4 decompress -i "$W/da.kt" --hier "$W/inf.tsv" --smooth --set acrocentric 2>/dev/null | cmp -s - "$W/inf.want" && ok "decompress --db --aln, 4 threads: same as plain" || no "decompress --db --aln: DIFFERS"
$KT --aln "$D/p.paf" --ref "$D/ref.bed" decompress -i "$W/d.kt" --hier "$W/inf.tsv" --smooth --set acrocentric 2>/dev/null | cmp -s - "$W/inf.wants" && ok "decompress --aln: same as plain" || no "decompress --aln: DIFFERS"
sq=$(head -1 "$D/db.acrocentric.bed" | cut -f1)
$KT query -i "$W/pa.kt" -r "$sq":0-100000 > "$W/inf.want" 2>/dev/null
$KT --db "$D/db.features.tsv" query -i "$W/f.kt" -r "$sq":0-100000 2>/dev/null > "$W/inf.got"
[ "$(awk -F'\t' '{print NF}' "$W/inf.got" | sort -u)" = 9 ] && ok "decompress --db: 6 label columns" || no "decompress --db: wrong column count"
awk -F'\t' 'NR>1 && $1==p && $2==e && $4==l {e=$3; next} NR>1 {print p"\t"s"\t"e"\t"l} {p=$1; s=$2; e=$3; l=$4} END{print p"\t"s"\t"e"\t"l}' <(cut -f1-4 "$W/inf.got") | cmp -s - "$W/inf.want" && ok "decompress --db: acrocentric column, merged, same as plain" || no "decompress --db: DIFFERS"
$KT --db "$D/db.features.tsv" --aln "$D/db.self.paf" --ref "$D/db.rall.bed" query -i "$W/da.kt" -r "$sq":0-100000 2>/dev/null | cmp -s - "$W/inf.got" && ok "decompress --db --aln: same as --db" || no "decompress --db --aln: DIFFERS"
$KT --aln "$D/p.paf" --ref "$D/ref.bed" query -i "$W/d.kt" -r q1 2>/dev/null | cmp -s <($KT query -i "$W/p.kt" -r q1 2>/dev/null) && ok "decompress --aln: same as plain" || no "decompress --aln: DIFFERS"
# --set: one set equals its own track
n=0; for f in $sets; do $KT --db "$D/db.features.tsv" --set $f decompress -i "$W/f.kt" 2>/dev/null | cmp -s - "$D/db.$f.bed" && n=$((n+1)); done
[ "$n" = 6 ] && ok "decompress --set: 6/6 equal to the set's track" || no "decompress --set: only $n/6"
$KT --db "$D/db.features.tsv" --aln "$D/db.self.paf" --ref "$D/db.rall.bed" --set gene decompress -i "$W/da.kt" 2>/dev/null | cmp -s - "$D/db.gene.bed" && ok "decompress --set on db+aln: same" || no "decompress --set on db+aln: DIFFERS"
n=0; for r in "-r $sq" "-r $sq:100-101" "-r $sq:0-100000"; do $KT --db "$D/db.features.tsv" --set acrocentric query $r -i "$W/f.kt" 2>/dev/null | cmp -s - <($KT query $r -i "$W/pa.kt" 2>/dev/null) && [ -n "$($KT query $r -i "$W/pa.kt" 2>/dev/null)" ] && n=$((n+1)); done
[ "$n" = 3 ] && ok "query --set: 3/3 forms equal the set's track" || no "query --set: only $n/3"
$KT --db "$D/db.features.tsv" --set gene,acrocentric decompress -i "$W/f.kt" 2>/dev/null > "$W/ga.bed"
[ "$(awk -F'\t' '{print NF}' "$W/ga.bed" | sort -u)" = 5 ] && awk -F'\t' 'NR>1 && $1==p && $2==e && $5==l {e=$3; next} NR>1 {print p"\t"s"\t"e"\t"l} {p=$1; s=$2; e=$3; l=$5} END{print p"\t"s"\t"e"\t"l}' "$W/ga.bed" | cmp -s - "$D/db.acrocentric.bed" && ok "decompress --set gene,acrocentric: two columns, merged jointly" || no "decompress --set gene,acrocentric: WRONG"
$KT --db "$D/db.features.tsv" decompress -i "$W/f.kt" 2>/dev/null | cmp -s - <($KT --db "$D/db.features.tsv" --set all decompress -i "$W/f.kt" 2>/dev/null) && ok "decompress --db on stdout: all sets, same as --set all" || no "decompress --db on stdout: DIFFERS from --set all"
# --hier: outside tuples get 2^26 + rank
{ printf 'feature_set\tchild\tparent\n'; for f in $sets; do for l in $(cut -f4 "$D/db.$f.bed" | sort -u | grep -v '^novel$'); do printf '%s\t%s\tcategorized\n' "$f" "$l"; done; done; printf 'repeat\tLINE\trepeat\nrepeat\trepeat\tcategorized\n'; } | awk -F'\t' '!s[$1 FS $2]++' > "$W/hv.tsv"
mkdir -p "$W/hv"; for f in $sets; do cp "$D/db.$f.bed" "$W/hv/$f.bed"; done
awk -F'\t' 'NR==5{$4="repeat"} {print}' OFS='\t' "$D/db.repeat.bed" > "$W/hv/repeat.bed"
hin=$(for f in $sets; do echo "-i $W/hv/$f.bed"; done)
must_fail "a tuple outside features.tsv without --hier" $KT compress --db "$D/db.features.tsv" $hin
$KT compress --db "$D/db.features.tsv" --hier "$W/hv.tsv" $hin > "$W/hv.kt" 2>/dev/null
$KT --raw decompress -i "$W/hv.kt" 2>/dev/null | awk -F'\t' '$4>=67108864' | grep -q . && ok "--hier: a vocabulary featureID stored" || no "--hier: NO vocabulary featureID"
# featureID 2^26-1 round trips, 2^26 refused
printf 'featureID\tacrocentric\tchromosome\n67108863\tacro\tchr1\n' > "$W/big.tsv"
printf 'a\t0\t10\t67108863\n' > "$W/big.bed"
$KT compress --db "$W/big.tsv" -i "$W/big.bed" > "$W/big.kt" 2>/dev/null && $KT decompress --db "$W/big.tsv" --set acrocentric -i "$W/big.kt" 2>/dev/null | grep -q acro \
  && ok "featureID 2^26-1: round trip" || no "featureID 2^26-1: WRONG"
printf 'featureID\tacrocentric\tchromosome\n67108864\tacro\tchr1\n' > "$W/over.tsv"
must_fail "featureID 2^26" $KT compress --db "$W/over.tsv" -i "$W/big.bed"
mkdir -p "$W/hvo"; $KT decompress --db "$D/db.features.tsv" --hier "$W/hv.tsv" --out-dir "$W/hvo" -i "$W/hv.kt" 2>/dev/null
n=0; for f in $sets; do cmp -s "$W/hvo/$f.bed" "$W/hv/$f.bed" && n=$((n+1)); done
[ "$n" = 6 ] && ok "--hier: 6/6 restored" || no "--hier: only $n/6 restored"
$KT compress --db "$D/db.features.tsv" --hier "$W/hv.tsv" $in 2>/dev/null | cmp -s - "$W/f.kt" && ok "--hier: a presmoothed track is byte-identical" || no "--hier: presmoothed track CHANGED"
must_fail "a vocabulary track decoded without --hier" $KT decompress --db "$D/db.features.tsv" --out-dir "$W/hvo" -i "$W/hv.kt"
sed '$d' "$W/hv.tsv" > "$W/hv2.tsv"
must_fail "a vocabulary track with another hierarchy" $KT decompress --db "$D/db.features.tsv" --hier "$W/hv2.tsv" --out-dir "$W/hvo" -i "$W/hv.kt"
must_fail "--hier without --db" $KT --hier "$W/hv.tsv" decompress -i "$W/hv.kt"
$KT compress < "$W/hv/repeat.bed" > "$W/hvp.kt" 2>/dev/null
$KT --db "$D/db.features.tsv" decompress -i "$W/hv.kt" --hier "$W/hv.tsv" --smooth --set repeat 2>/dev/null | cmp -s - <($KT decompress -i "$W/hvp.kt" --hier "$W/hv.tsv" --smooth --set repeat 2>/dev/null) && ok "--hier: decompress --db same as the plain track" || no "--hier: decompress --db DIFFERS"
$KT --db "$D/db.features.tsv" --hier "$W/hv.tsv" --set repeat decompress -i "$W/hv.kt" 2>/dev/null | cmp -s - "$W/hv/repeat.bed" && ok "--hier: --set repeat gives the input" || no "--hier: --set repeat DIFFERS"
mkdir -p "$W/so"; $KT --db "$D/db.features.tsv" --set repeat,gene --out-dir "$W/so/new" decompress -i "$W/f.kt" 2>/dev/null; [ "$(ls "$W/so/new" | tr '\n' ' ')" = "gene.bed repeat.bed " ] && cmp -s "$W/so/new/gene.bed" "$D/db.gene.bed" && ok "decompress --set --out: only the picked sets, directory created" || no "decompress --set --out: WRONG"
# index needs no --aln or --db companions
cp "$W/d.kt" "$W/di.kt"; cp "$W/f.kt" "$W/fi.kt"
$KT index -i "$W/di.kt" >/dev/null 2>&1 && $KT index -i "$W/fi.kt" >/dev/null 2>&1 && $KT --aln "$D/p.paf" --ref "$D/ref.bed" query -i "$W/di.kt" -r q1 2>/dev/null | cmp -s - <($KT --aln "$D/p.paf" --ref "$D/ref.bed" query -i "$W/d.kt" -r q1 2>/dev/null) && ok "index on --aln and --db tracks: built, same answers" || no "index on --aln and --db tracks: FAILED"
must_fail "--set without --db" $KT --set gene decompress -i "$W/pa.kt"
must_fail "--set not in features.tsv" $KT --db "$D/db.features.tsv" --set gene,nosuch decompress -i "$W/f.kt"
must_fail "--set names a set twice" $KT --db "$D/db.features.tsv" --set gene,gene decompress -i "$W/f.kt"
must_fail "--set all with a name" $KT --db "$D/db.features.tsv" --set all,gene decompress -i "$W/f.kt"
{ printf 'feature_set\tchild\tparent\n'; for s in $(head -1 "$D/db.features.tsv" | cut -f2-); do cut -f4 "$D/db.$s.bed" | sort -u | grep -vx novel | sed "s/^/$s\t/; s/\$/\tcategorized/"; done; } > "$W/full.tsv"
rm -rf "$W/sod"; $KT --db "$D/db.features.tsv" decompress -i "$W/f.kt" --hier "$W/full.tsv" --smooth --out-dir "$W/sod" 2>/dev/null; sok=1
for s in $(head -1 "$D/db.features.tsv" | cut -f2-); do [ -s "$W/sod/$s.bed" ] || sok=0
  $KT --db "$D/db.features.tsv" decompress -i "$W/f.kt" --hier "$W/full.tsv" --smooth --set $s 2>/dev/null | cmp -s - "$W/sod/$s.bed" || sok=0; done
[ $sok = 1 ] && ok "--smooth --out-dir: every set, each as --set gives it" || no "--smooth --out-dir: WRONG"
must_fail "--smooth of a -d track, no --set nor --out-dir" $KT --db "$D/db.features.tsv" decompress -i "$W/hv.kt" --hier "$W/hv.tsv" --smooth
must_fail "--smooth --out-dir without --db" $KT decompress -i "$W/p.kt" --hier "$W/inf.tsv" --smooth --out-dir "$W/sod2"
must_fail "query --smooth without -r" $KT query -i "$W/p.kt" --hier "$W/inf.tsv" --smooth --set acrocentric
must_fail "query --smooth --out-dir" $KT --db "$D/db.features.tsv" query -i "$W/hv.kt" -r x --hier "$W/hv.tsv" --smooth --out-dir "$W/sod2"
must_fail "--smooth with two sets" $KT --db "$D/db.features.tsv" decompress -i "$W/f.kt" --hier "$W/inf.tsv" --smooth --set gene,acrocentric
must_fail "decompress an --aln track without --aln" $KT decompress -i "$W/d.kt" --hier "$W/inf.tsv" --smooth --set acrocentric
must_fail "decompress a --db track without --db" $KT decompress -i "$W/f.kt" --hier "$W/inf.tsv" --smooth --set acrocentric
must_fail "query --db on a plain track" $KT --db "$D/db.features.tsv" query -i "$W/p.kt" -r q1
must_fail "decompress --aln without --ref" $KT --aln "$D/p.paf" decompress -i "$W/d.kt" --hier "$W/inf.tsv" --smooth --set acrocentric
must_fail "index --aln" $KT --aln "$D/p.paf" --ref "$D/ref.bed" index -i "$W/d.kt"

# ordinal names reach the PAF via --aln-names
$KT compress --aln "$D/p.paf" --ref "$D/ref.bed" --aln-names "$D/ord.names.txt" --names drop < "$D/ord.bed" > "$W/ord.kt" 2>/dev/null
$KT decompress --aln "$D/p.paf" --ref "$D/ref.bed" --aln-names "$D/ord.names.txt" -i "$W/ord.kt" 2>/dev/null > "$W/ord.out"
cmp -s "$D/ord.bed" "$W/ord.out" && ok "aln --names drop with an ordinal map: round-trips" || no "aln --names drop with an ordinal map: DIFFERS"

# gaps in the ordinals, names dropped
$KT compress --aln "$D/p.paf" --ref "$D/ref.bed" --aln-names "$D/ord2.names.txt" --names drop < "$D/ord2.bed" > "$W/ord2.kt" 2>/dev/null
$KT decompress --aln "$D/p.paf" --ref "$D/ref.bed" --aln-names "$D/ord2.names.txt" -i "$W/ord2.kt" 2>/dev/null > "$W/ord2.out"
cmp -s "$D/ord2.bed" "$W/ord2.out" && ok "aln --names drop with gaps in the ordinals: round-trips" || no "aln --names drop with gaps in the ordinals: DIFFERS"

# runs ending on reference run ends
$KT compress --aln "$D/lp.paf" --ref "$D/lref.bed" < "$D/ls.bed" \
  > "$W/ls.kt" 2>/dev/null
$KT decompress --aln "$D/lp.paf" --ref "$D/lref.bed" -i "$W/ls.kt" 2>/dev/null > "$W/ls.out"
cmp -s "$D/ls.bed" "$W/ls.out" && ok "lengths the reference implies: round-trip" || no "lengths the reference implies: DIFFERS"
must_fail "decompress a length-delta track without --aln" $KT decompress -i "$W/ls.kt"
# 3000 runs, all on reference run ends
awk 'BEGIN{srand(5);p=0;for(i=0;i<3000;i++){c=p+int(rand()*900)+1;printf "chrA\t%d\t%d\tL%d\n",p,c,i%4;p=c}print p > "'"$W"'/len"}' > "$W/big.ref"
L=$(cat "$W/len"); sed 's/^chrA/q/' "$W/big.ref" > "$W/big.q"
printf "q\t$L\t0\t$L\t+\tchrA\t$L\t0\t$L\t$L\t$L\t60\tcg:Z:${L}M\n" > "$W/big.paf"
$KT compress --aln "$W/big.paf" --ref "$W/big.ref" < "$W/big.q" > "$W/big.kt" 2>/dev/null
n=$($KT decompress --raw -i "$W/big.kt" 2>/dev/null | awk -F'\t' '$2==$3' | wc -l)
[ "$n" = 3000 ] && ok "lengths left out: 3000 of 3000 rows" || no "lengths left out: $n of 3000, expected 3000"
$KT decompress --aln "$W/big.paf" --ref "$W/big.ref" -i "$W/big.kt" 2>/dev/null | cmp -s - "$W/big.q" && ok "sentinel block: round-trip" || no "sentinel block: DIFFERS"
# a block where plain lengths win
awk 'BEGIN{srand(6);p=0;for(i=0;i<3000;i++){c=p+int(rand()*900)+1;printf "q\t%d\t%d\tL%d\n",p,c,i%4;p=c}}' > "$W/off.q"
$KT compress --aln "$W/big.paf" --ref "$W/big.ref" < "$W/off.q" > "$W/off.kt" 2>/dev/null
$KT decompress --aln "$W/big.paf" --ref "$W/big.ref" -i "$W/off.kt" 2>/dev/null | cmp -s - "$W/off.q" && ok "plain-length block under --aln: round-trip" || no "plain-length block under --aln: DIFFERS"

# one-length reads: last run from the total
python3 - > "$W/tail.bed" <<'PY'
import random
random.seed(5)
for s in range(4000):
    cuts = sorted(random.sample(range(1, 116), 3)); p = 0
    for c in cuts + [116]:
        print(f"r{s}\t{p}\t{c}\tL{random.randint(0,5)}"); p = c
PY
$KT compress --block 4096 < "$W/tail.bed" > "$W/tail.kt" 2>/dev/null
$KT decompress -i "$W/tail.kt" 2>/dev/null > "$W/tail.out"
cmp -s "$W/tail.bed" "$W/tail.out" && ok "constant-length reads: round-trip" || no "constant-length reads: DIFFERS"
# random access lands on the right rows
$KT query -i "$W/tail.kt" -r r2000 2>/dev/null > "$W/tq.out"
grep "^r2000	" "$W/tail.bed" > "$W/tq.exp"
cmp -s "$W/tq.out" "$W/tq.exp" && ok "constant-length reads: query exact" || no "constant-length reads: query WRONG"

# mixed-length reads: storing the total loses
python3 - > "$W/vlen.bed" <<'PY'
import random
random.seed(5)
for s in range(4000):
    L = 116 - random.randint(0, 80)
    cuts = sorted(random.sample(range(1, L), 3)); p = 0
    for c in cuts + [L]:
        print(f"r{s}\t{p}\t{c}\tL{random.randint(0,5)}"); p = c
PY
$KT compress --block 4096 < "$W/vlen.bed" > "$W/vlen.kt" 2>/dev/null
$KT decompress -i "$W/vlen.kt" 2>/dev/null | cmp -s - "$W/vlen.bed" && ok "reads of mixed length: round-trip" || no "reads of mixed length: DIFFERS"

# one-length reads with 1 in 51 strays
python3 - > "$W/out.bed" <<'PY'
import random
random.seed(11)
for s in range(4000):
    L = random.choice([8, 30, 69, 105]) if random.randint(0, 50) == 0 else 116
    cuts = sorted(random.sample(range(1, L), 3)); p = 0
    for c in cuts + [L]:
        print(f"r{s}\t{p}\t{c}\tL{random.randint(0,5)}"); p = c
PY
$KT compress --block 4096 < "$W/out.bed" > "$W/out.kt" 2>/dev/null
$KT decompress -i "$W/out.kt" 2>/dev/null | cmp -s - "$W/out.bed" && ok "uniform reads with stray lengths: round-trip" || no "uniform reads with stray lengths: DIFFERS"

# interleaved records refused
printf 'a\t0\t10\tA\nb\t0\t10\tB\na\t10\t20\tC\n' > "$W/il.bed"
must_fail "interleaved records" $KT compress --block 4 < "$W/il.bed"

# gaps and overlaps refused
printf 'a\t0\t10\tA\na\t20\t30\tB\n' > "$W/gap.bed"
printf 'a\t0\t10\tA\na\t9\t20\tB\n'  > "$W/ovl.bed"
for f in gap ovl; do must_fail "$f" $KT compress --block 4 < "$W/$f.bed"; done

# names that are already 0,1,2,...
printf '0\t0\t10\tA\n1\t0\t10\tB\n2\t0\t10\tC\n' > "$W/ord2.bed"
$KT compress --block 4 < "$W/ord2.bed" > "$W/ord2.kt" 2>/dev/null
$KT decompress -i "$W/ord2.kt" 2>/dev/null | cmp -s - "$W/ord2.bed" && ok "ordinal names: round-trip" || no "ordinal names: DIFFERS"

# a last line with no newline: refused
printf 'a\t0\t10\tA' > "$W/nonl.bed"
must_fail "input without a trailing newline" $KT compress --block 4 < "$W/nonl.bed"


# --db parser refuses what plain parser refuses
for bad in 'a\tx\t10\t1' '\t0\t10\t1' 'a\t0\t18446744073709551616\t1' 'a\t0\t20\t1\na\t20\t10\t2'; do
  printf "$bad\n" > "$W/dbbad.bed"
  must_fail "--db malformed input" $KT compress --db "$D/db.features.tsv" -i "$W/dbbad.bed"
done

# threaded smooth equals serial on 10,000 sequences
python3 - > "$W/sm.bed" <<'PY'
for i in range(10000):
    print("r%d\t0\t50\tacro" % i); print("r%d\t50\t100\tnovel" % i)
PY
printf 'feature_set\tchild\tparent\nacrocentric\tacro\tcategorized\nacrocentric\tnovel\tcategorized\n' > "$W/sm.tsv"
must_fail "decompress --db a set not in features.tsv" $KT --db "$D/db.features.tsv" decompress -i "$W/f.kt" --hier "$W/sm.tsv" --smooth --set nosuch
$KT compress --block 32 < "$W/sm.bed" > "$W/sm.kt" 2>/dev/null
$KT --threads 1 decompress -i "$W/sm.kt" --hier "$W/sm.tsv" --smooth --set acrocentric > "$W/sm.1" 2>/dev/null || no "decompress, serial: FAILED"
for t in 2 4 8; do
  $KT --threads $t decompress -i "$W/sm.kt" --hier "$W/sm.tsv" --smooth --set acrocentric > "$W/sm.$t" 2>/dev/null || { echo "decompress, $t threads: FAILED"; fail=1; continue; }
  cmp -s "$W/sm.1" "$W/sm.$t" && ok "decompress, $t threads: identical to serial" || no "decompress, $t threads: DIFFERS from serial"
done

# malformed CIGARs refused, a good one accepted
T1=$(cut -f1 "$D/ref.bed" | head -1); printf 'a\t0\t10\tA\n' > "$W/cg.bed"
for cg in "cg:Z:M" "cg:Z:999999999999999999999999M" "cg:Z:5M5Q5M" "cg:Z:5M" "cg:Z:0M10M"; do
  printf "a\t10\t0\t10\t+\t$T1\t100\t0\t10\t10\t10\t60\t$cg\n" > "$W/cg.paf"
  must_fail "bad CIGAR $cg" $KT compress --aln "$W/cg.paf" --ref "$D/ref.bed" < "$W/cg.bed"
done
printf "a\t10\t0\t10\t+\t$T1\t100\t0\t10\t10\t10\t60\tcg:Z:10M\n" > "$W/cg.paf"
$KT compress --aln "$W/cg.paf" --ref "$D/ref.bed" < "$W/cg.bed" > /dev/null 2>&1 && ok "good CIGAR: accepted" || no "good CIGAR: REFUSED"

# first audit: one case per finding
awk 'BEGIN{for(i=0;i<20000;i++)printf "s\t%d\t%d\tA\n",i*10,i*10+10}' > "$W/rep.bed"
$KT compress < "$W/rep.bed" > "$W/rep.kt" 2>/dev/null
$KT decompress -i "$W/rep.kt" 2>/dev/null | cmp -s - "$W/rep.bed" && ok "repetitive track: byte-identical" || no "repetitive track: WRONG"
printf 'feature_set\tchild\tparent\nacrocentric\tA\tcategorized\n' > "$W/rep.tsv"
$KT decompress -i "$W/rep.kt" --hier "$W/rep.tsv" --smooth --set acrocentric > /dev/null 2>&1 && ok "repetitive decompress: ok" || no "repetitive decompress: FAILED"
{ head -1 "$D/db.features.tsv"; printf '0\tacro\tchr1\tintergenic\tarm\tLINE\tITS\n'; tail -n +2 "$D/db.features.tsv"; } > "$W/id0.tsv"
must_fail "features.tsv defining id 0" $KT compress --db "$W/id0.tsv" $in
LR=$(head -c 70000 /dev/zero | tr '\0' A)
printf 'long\t0\t10\tA\nnext\t0\t10\tB\n' > "$W/lr.bed"
printf '@long\n%s\n+\n%s\n@next\nACGT\n+\nIIII\n' "$LR" "$LR" > "$W/lr.fq"
$KT compress --names drop < "$W/lr.bed" > "$W/lr.kt" 2>/dev/null
$KT --names "$W/lr.fq" decompress -i "$W/lr.kt" 2>/dev/null | cmp -s - "$W/lr.bed" && ok "long FASTQ read: names restored" || no "long FASTQ read: WRONG"
printf '0\t0\t10\tA\n1\t0\t10\tB\n' > "$W/ord.bed"
LL=$(head -c 300 /dev/zero | tr '\0' L)_specific
printf 's\t0\t10\t%s\n' "$LL" > "$W/ll.bed"
$KT compress < "$W/ll.bed" > "$W/ll.kt" 2>/dev/null
$KT decompress -i "$W/ll.kt" > /dev/null 2>&1 && ok "long label v1: ok" || no "long label v1: CRASHED"

# second audit: one case per finding
printf 'chr1\t0\t100\tA\nchr1\t100\t200\tB\nchr1\t200\t300\tC\nchr1\t300\t400\tD\nchr1\t400\t500\tE\nchr1\t500\t600\tF\n' > "$W/sa.ref"
printf '1\t500\t0\t500\t+\tchr1\t600\t0\t550\t500\t550\t60\tcg:Z:100M50D400M\n0\t500\t0\t500\t+\tchr1\t600\t0\t500\t500\t500\t60\tcg:Z:500M\n' > "$W/sa.paf"
printf '1\t0\t100\tA\n1\t100\t200\tB\n1\t200\t300\tC\n1\t300\t400\tD\n1\t400\t500\tE\n0\t0\t100\tA\n0\t100\t200\tB\n0\t200\t300\tC\n0\t300\t400\tD\n0\t400\t500\tE\n' > "$W/sa.bed"
printf 'featureID\tacrocentric\tchromosome\n1\tacro\tchr1\n2\tacro\n' > "$W/short.tsv"
printf 's\t0\t10\t1\ns\t10\t20\t2\n' > "$W/ids.bed"
must_fail "features.tsv row with a missing label" $KT compress --db "$W/short.tsv" -i "$W/ids.bed"
printf 'read one\t0\t10\tA\nread two\t0\t10\tB\n' > "$W/sp.bed"; printf 'read one\nread two\n' > "$W/sp.names.txt"
$KT compress --names drop < "$W/sp.bed" > "$W/sp.kt" 2>/dev/null
$KT --names "$W/sp.names.txt" decompress -i "$W/sp.kt" 2>/dev/null | cmp -s - "$W/sp.bed" && ok "name list with spaces: restored" || no "name list with spaces: WRONG"
cp "$D/s.bed" "$W/keep.bed"; $KT compress "$W/keep.bed" < /dev/null > /dev/null 2>&1 && { echo "compress FILE: ACCEPTED"; fail=1; }
cmp -s "$W/keep.bed" "$D/s.bed" && ok "compress FILE: refused, input intact" || no "compress FILE: INPUT DESTROYED"

# third audit: one case per finding
awk 'BEGIN{for(i=0;i<20;i++)printf "r%d\t0\t%d\tL%d\n",i,50+i,i%5}' > "$W/idx.bed"
$KT compress < "$W/idx.bed" > "$W/idx.kt" 2>/dev/null; $KT index -i "$W/idx.kt" >/dev/null 2>&1
tac "$W/idx.bed" | $KT compress > "$W/idx.kt" 2>/dev/null
$KT query -i "$W/idx.kt" -r r7 2>/dev/null | cmp -s - <(grep -P '^r7\t' "$W/idx.bed") && ok "stale idx: ignored" || no "stale idx: WRONG ROWS"
$KT index -i "$W/idx.kt" >/dev/null 2>&1; $KT query -i "$W/idx.kt" -r r7 2>/dev/null | cmp -s - <(grep -P '^r7\t' "$W/idx.bed") && ok "fresh idx: correct" || no "fresh idx: WRONG ROWS"
T1=$(cut -f1 "$D/ref.bed" | head -1)
printf 'a\t0\t5\t=\na\t5\t10\tZ\n' > "$W/eq.bed"
printf "a\t10\t0\t10\t+\t$T1\t100\t0\t10\t10\t10\t60\tcg:Z:10M\n" > "$W/eq.paf"
must_fail "label = under --aln" $KT compress --aln "$W/eq.paf" --ref "$D/ref.bed" < "$W/eq.bed"
python3 -c "
for i in range(70): print(f's{i}\t0\t100\tA'); print(f's{i}\t100\t{2**63+100}\tB')" > "$W/t63.bed"
$KT compress < "$W/t63.bed" > "$W/t63.kt" 2>/dev/null
$KT decompress -i "$W/t63.kt" 2>/dev/null | cmp -s - "$W/t63.bed" && ok "total past 2^63: byte-identical" || no "total past 2^63: WRONG"
printf 'feature_set\tchild\tparent\nacrocentric\tacro\tx\nacrocentric\tx\tacro\n' > "$W/cyc.tsv"
rc=0; timeout 10 $KT decompress -i "$W/p.kt" --hier "$W/cyc.tsv" --smooth --set acrocentric >/dev/null 2>&1 || rc=$?
[ $rc = 1 ] && ok "cyclic hierarchy: refused" || no "cyclic hierarchy: rc=$rc"
for c in "decompress -i $W/p.kt" "query -i $W/p.kt -r q1" "decompress -i $W/sm.kt --hier $W/sm.tsv --smooth --set acrocentric" "decompress --raw -i $W/p.kt"; do
  $KT $c > /dev/full 2>/dev/null && { echo "$c to /dev/full: EXIT 0"; fail=1; } || echo "${c%% *} to /dev/full: exit 1"
done
printf 'a\t0\t10\tA\0B\nc\t0\t10\tB\n' > "$W/nul.bed"
must_fail "NUL inside a row" $KT compress < "$W/nul.bed"
sed '2s/acro/ACRO/' "$D/db.features.tsv" > "$W/other.tsv"
mkdir -p "$W/o2"; must_fail "decompress --db with another features.tsv" $KT decompress --db "$W/other.tsv" --out-dir "$W/o2" -i "$W/f.kt"

# blocks mixing plain lengths and the sentinel
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
$KT compress --aln "$W/mx.paf" --ref "$W/mx.ref" --block 3000 < "$W/mx.bed" > "$W/mx.kt" 2>/dev/null
$KT decompress --aln "$W/mx.paf" --ref "$W/mx.ref" -i "$W/mx.kt" 2>/dev/null | cmp -s - "$W/mx.bed" \
  && ok "mixed sentinel blocks, seed $seed: byte-identical" || no "mixed sentinel blocks, seed $seed: WRONG"
done
# minus-strand sentinel on reference run starts
printf 'chrB\t0\t100\tA\nchrB\t100\t250\tB\nchrB\t250\t400\tC\n' > "$W/mn.ref"
printf 'm\t400\t0\t400\t-\tchrB\t400\t0\t400\t400\t400\t60\tcg:Z:400M\n' > "$W/mn.paf"
printf 'm\t0\t150\tX\nm\t150\t300\tY\nm\t300\t400\tZ\n' > "$W/mn.bed"
$KT compress --aln "$W/mn.paf" --ref "$W/mn.ref" < "$W/mn.bed" > "$W/mn.kt" 2>/dev/null
n=$($KT decompress --raw -i "$W/mn.kt" 2>/dev/null | awk -F'\t' '$2==$3' | wc -l)
$KT decompress --aln "$W/mn.paf" --ref "$W/mn.ref" -i "$W/mn.kt" 2>/dev/null | cmp -s - "$W/mn.bed" && [ "$n" = 3 ] \
  && ok "minus-strand sentinel: 3 of 3 rows, round-trip" || no "minus-strand sentinel: $n of 3, WRONG"

# NOSZ codecs carry no size: refused

# flags with no command: usage, no crash
must_fail "--threads with no command" $KT --threads 4
must_fail "--tmp-dir with no command" $KT --tmp-dir /tmp
must_fail "--effort low with no command" $KT --effort low
$KT --profile 2>&1 | grep -q '^karyozip: usage:' && ok "flags with no command: usage, not a crash" || no "flags with no command: WRONG"
python3 "$D/kt.py" "$W/p.kt" "$W/nosz.kt" lo9 0x10
must_fail "a header naming a NOSZ codec" $KT decompress -i "$W/nosz.kt"

# plain lengths near 2^31 under --aln
python3 - "$W" <<'PY'
import random, sys
W = sys.argv[1]; random.seed(9); p = 0; ref = []
for i in range(3000):
    c = p + random.randint(1 << 28, (1 << 31) - 1); ref.append((p, c)); p = c
open(W + "/huge.ref", "w").write("".join("chrH\t%d\t%d\tL%d\n" % (a, b, i % 3) for i, (a, b) in enumerate(ref)))
open(W + "/huge.q", "w").write("".join("h\t%d\t%d\tL%d\n" % (a, b, i % 3) for i, (a, b) in enumerate(ref)))
cg = "".join("%dM" % (b - a) for a, b in ref)   # one op per run, each under 2^32
open(W + "/huge.paf", "w").write("h\t%d\t0\t%d\t+\tchrH\t%d\t0\t%d\t%d\t%d\t60\tcg:Z:%s\n" % ((p,) * 6 + (cg,)))
PY
$KT compress --aln "$W/huge.paf" --ref "$W/huge.ref" < "$W/huge.q" > "$W/huge.kt" 2>/dev/null && $KT decompress --aln "$W/huge.paf" --ref "$W/huge.ref" -i "$W/huge.kt" 2>/dev/null | cmp -s - "$W/huge.q" \
  && ok "fast with long plain lengths: byte-identical" || no "fast with long plain lengths: FAILED"
# codec id with bit 31 set: refused
python3 "$D/kt.py" "$W/p.kt" "$W/neg.kt" lo9 0x80010010
must_fail "a codec id with bit 31 set" $KT decompress -i "$W/neg.kt"

# a packed alignment needs the same reference
{ cat "$D/ref.bed"; printf 'chr2\t0\t50\tB\n'; } > "$W/ref2.bed"; { printf 'chr2\t0\t50\tB\n'; cat "$D/ref.bed"; } > "$W/ref.rev"
$KT aln-pack -a "$D/p.paf" --ref "$W/ref2.bed" -o "$W/p.kta" >/dev/null 2>&1
$KT compress --aln "$W/p.kta" --ref "$W/ref2.bed" < "$D/s.bed" > "$W/pk.kt" 2>/dev/null && $KT decompress --aln "$W/p.kta" --ref "$W/ref2.bed" -i "$W/pk.kt" 2>/dev/null | cmp -s - "$D/s.bed" \
  && ok "packed alignment: round-trip" || no "packed alignment: WRONG"
must_fail "packed alignment with a reordered reference" $KT decompress --aln "$W/p.kta" --ref "$W/ref.rev" -i "$W/pk.kt"
{ printf 'KSTAL2'; tail -c +7 "$W/p.kta"; } > "$W/old.kta"
must_fail "a packed alignment with another magic" $KT decompress --aln "$W/old.kta" --ref "$W/ref2.bed" -i "$W/pk.kt"

# label count with no label table: refused
$KT compress < /dev/null > "$W/mt.kt" 2>/dev/null
python3 "$D/kt.py" "$W/mt.kt" "$W/mt2.kt" 7 0x4000000000000000
must_fail "a label count with no label table" $KT decompress -i "$W/mt2.kt" --hier "$W/sm.tsv" --smooth --set acrocentric

# input-side findings, one case each
printf 'featureID\tacrocentric\tchromosome\n1\tA\tX\n1\tB\tY\n' > "$W/dupid.tsv"
printf 's\t0\t10\tA\n' > "$W/d1.bed"; printf 's\t0\t10\tX\n' > "$W/d2.bed"
must_fail "features.tsv with a featureID twice" $KT compress --db "$W/dupid.tsv" -i "$W/d1.bed" -i "$W/d2.bed"
printf 'featureID\tx\tx\n1\tA\tX\n' > "$W/dupset.tsv"
must_fail "features.tsv with a set named twice" $KT compress --db "$W/dupset.tsv" -i "$W/d1.bed" -i "$W/d2.bed"
printf '>r1 extra\t0\t10\tA\n>r2\t0\t10\tB\n' > "$W/gt.bed"; printf '>r1 extra\n>r2\n' > "$W/gt.names"
$KT compress --names drop < "$W/gt.bed" > "$W/gt.kt" 2>/dev/null
$KT --names "$W/gt.names" decompress -i "$W/gt.kt" 2>/dev/null | cmp -s - "$W/gt.bed" && ok "name list starting with >: restored" || no "name list starting with >: WRONG"
L="$W"; for i in 1 2 3 4 5; do L="$L/$(printf 'd%.0s' $(seq 1 120))"; done; mkdir -p "$L"
echo keep > "${L:0:511}" 2>/dev/null || true
awk 'BEGIN{for(i=0;i<80000;i++)printf "s%d\t0\t10\tA\n",i}' > "$W/long.bed"
$KT compress --tmp-dir "$L" --names drop < "$W/long.bed" > "$W/long.kt" 2>/dev/null \
  && $KT decompress --names <(cut -f1 "$W/long.bed") -i "$W/long.kt" 2>/dev/null | cmp -s - "$W/long.bed" && [ "$(cat "${L:0:511}" 2>/dev/null)" = keep ] \
  && ok "long temp dir: round-trip, bystander intact" || no "long temp dir: WRONG"
$KT compress < "$D/s.bed" >&- 2>/dev/null && { echo "compress with stdout closed: EXIT 0"; fail=1; } || echo "compress with stdout closed: refused"
$KT compress < "$W" > "$W/dir.kt" 2>/dev/null && { echo "compress from a directory: EXIT 0"; fail=1; } || echo "compress from a directory: refused"
printf 'a\t0\t10\tA\nb\t0\t10\tB\na\t0\t10\tC\n' > "$W/dupname.bed"
must_fail "a sequence name that repeats later" $KT compress < "$W/dupname.bed"
N512=$(head -c 512 /dev/zero | tr '\0' n); printf '%s\t0\t10\tA\n' "$N512" > "$W/n512.bed"
must_fail "a 512-byte name" $KT compress < "$W/n512.bed"

# round 14 findings, one case each
printf 'FQ6SH1h576F\t0\t10\tL\nOyHUvXChfZK\t0\t10\tL\nzz\t0\t10\tL\nOyHUvXChfZK\t0\t10\tL\n' > "$W/coll.bed"
must_fail "a repeat behind a hash collision" $KT compress < "$W/coll.bed"
printf 'FQ6SH1h576F\t0\t10\tL\nOyHUvXChfZK\t0\t10\tL\n' > "$W/coll2.bed"
$KT compress < "$W/coll2.bed" > "$W/coll2.kt" 2>/dev/null; $KT decompress -i "$W/coll2.kt" 2>/dev/null | cmp -s - "$W/coll2.bed" \
  && ok "two names with one hash: accepted" || no "two names with one hash: WRONG"
printf 'r1\r\nr2\r\n' > "$W/crlf.txt"; printf 'r1\t0\t10\tA\nr2\t0\t10\tB\n' > "$W/n2.bed"
$KT compress --names drop < "$W/n2.bed" > "$W/n2.kt" 2>/dev/null
$KT --names "$W/crlf.txt" decompress -i "$W/n2.kt" 2>/dev/null | cmp -s - "$W/n2.bed" && ok "CRLF name list: restored" || no "CRLF name list: WRONG"
printf '>r1 extra\t0\t10\tA\n' > "$W/gt1.bed"; printf '>r1 extra\n' > "$W/gt1.txt"
$KT compress --names drop < "$W/gt1.bed" > "$W/gt1.kt" 2>/dev/null
$KT --names "$W/gt1.txt" decompress -i "$W/gt1.kt" 2>/dev/null | cmp -s - "$W/gt1.bed" && ok "one-line name list starting with >: restored" || no "one-line name list starting with >: WRONG"
printf 'r1\n' > "$W/dupn.txt"; printf 'r1\nr1\n' > "$W/dupn2.txt"
must_fail "a name file with a repeat, on query" $KT --names "$W/dupn2.txt" query -i "$W/n2.kt" -r r1
printf 'featureID acrocentric chromosome\n1\tA\tX\n' > "$W/notab.tsv"; printf 's\t0\t10\t1\n' > "$W/ids1.bed"
must_fail "features.tsv whose first line has no tab" $KT compress --db "$W/notab.tsv" -i "$W/ids1.bed"

# round 15 findings, one case each
printf 'aaa\t0\t10\tA\nYo5P2GphlrB\t0\t10\tB\nYo5P2GphlrB\t10\t30\tC\nzzz\t0\t10\tD\n' > "$W/hc.bed"
$KT compress < "$W/hc.bed" > "$W/hc.kt" 2>/dev/null; $KT index -i "$W/hc.kt" >/dev/null 2>&1
must_fail "indexed query of a colliding absent name" $KT query -i "$W/hc.kt" -r TQyOn5lCDT5
$KT query -i "$W/hc.kt" -r Yo5P2GphlrB 2>/dev/null | cmp -s - <(grep -P '^Yo5P2GphlrB\t' "$W/hc.bed") && ok "indexed query of the stored name: correct" || no "indexed query of the stored name: WRONG"
printf 'r0\t0\t10\tA\nr1\t0\t10\tB\nr2\t0\t10\tC\nr3\t0\t10\tD\n' > "$W/nx.bed"
$KT compress --names drop < "$W/nx.bed" > "$W/nx.kt" 2>/dev/null
printf 'r0\nr1\nr2\nr3\n' > "$W/A.txt"; printf 'r3\nr2\nr1\nr0\n' > "$W/B.txt"
$KT --names "$W/A.txt" index -i "$W/nx.kt" >/dev/null 2>&1
[ "$($KT --names "$W/B.txt" query -i "$W/nx.kt" -r r0 2>/dev/null | cut -f4)" = D ] && ok "index under another name file: ignored" || no "index under another name file: WRONG"

# round 16 findings, one case each
cp "$W/B.txt" "$W/A.txt"
[ "$($KT --names "$W/A.txt" query -i "$W/nx.kt" -r r0 2>/dev/null | cut -f4)" = D ] && ok "index after the name file changed: ignored" || no "index after the name file changed: WRONG"
python3 - "$W" <<'PY'
import struct, sys
# swap the index's first and last ordinal
W = sys.argv[1]; d = bytearray(open(W + "/hc.kt.kzi", "rb").read())
n, k, _, pay = struct.unpack_from("<4Q", d, len(d) - 40)
nb = 1 << k; base = len(d) - 40 - pay + nb + (nb + 1023) // 1024 * 8
w = max(1, (n - 1).bit_length()); ent = int.from_bytes(d[base:len(d) - 48], "little")
get = lambda i: (ent >> (i * (32 + w) + 32)) & ((1 << w) - 1)
put = lambda e, i, v: (e & ~(((1 << w) - 1) << (i * (32 + w) + 32))) | (v << (i * (32 + w) + 32))
a, z = get(0), get(n - 1); ent = put(put(ent, 0, z), n - 1, a)
d[base:len(d) - 48] = ent.to_bytes(len(d) - 48 - base, "little"); open(W + "/hc.kt.kzi", "wb").write(d)
PY
$KT query -i "$W/hc.kt" -r aaa 2>/dev/null | cmp -s - <(grep -P '^aaa\t' "$W/hc.bed") && ok "index with edited ordinals: not trusted" || no "index with edited ordinals: WRONG"
python3 -c "open('$W/five.bed', 'w').write('s\t0\t10\t%s\n' % ('z' * (5 << 20)))"
$KT compress < "$W/five.bed" > "$W/five.kt" 2>/dev/null
printf 'feature_set\tchild\tparent\nfs\t%s\tcategorized\n' "$(head -c $((5 << 20)) /dev/zero | tr '\0' z)" > "$W/five.tsv"
$KT decompress -i "$W/five.kt" --hier "$W/five.tsv" --smooth --set fs 2>/dev/null | cmp -s - "$W/five.bed" && $KT decompress -i "$W/five.kt" 2>/dev/null | cmp -s - "$W/five.bed" \
  && $KT query -i "$W/five.kt" -r s 2>/dev/null | cmp -s - "$W/five.bed" && ok "5 MB label: decompress, decompress, query" || no "5 MB label: WRONG"
$KT compress --names drop < "$D/s.bed" > "$W/nn.kt" 2>/dev/null
$KT query --names <(cut -f1 "$D/s.bed" | uniq) -i "$W/nn.kt" -r q1 2>/dev/null | cmp -s - <(grep -P '^q1\t' "$D/s.bed") \
  && ok "query --names FILE on dropped names: correct" || no "query --names FILE on dropped names: WRONG"
must_fail "decompress of dropped names without --names" $KT decompress -i "$W/nn.kt"
must_fail "query of dropped names without --names" $KT query -i "$W/nn.kt" -r 0
must_fail "decompress --names ordinal (no such file)" $KT decompress --names ordinal -i "$W/nn.kt"

# round 17 findings, one case each
N300=$(head -c 300 /dev/zero | tr '\0' n); python3 -c "open('$W/fivename.bed', 'w').write('$N300\t0\t10\t%s\n' % ('z' * (5 << 20)))"
$KT compress < "$W/fivename.bed" > "$W/fivename.kt" 2>/dev/null
$KT decompress -i "$W/fivename.kt" 2>/dev/null | cmp -s - "$W/fivename.bed" && ok "5 MB label, long name: decompress" || no "5 MB label, long name: WRONG"
python3 -c "open('$W/biglab.tsv', 'w').write('featureID\tacrocentric\tchromosome\n1\t%s\tX\n2\tB\tY\n' % ('L' * (2 << 20)))"
printf 's\t0\t10\t1\ns\t10\t20\t2\ns\t20\t30\t1\n' > "$W/bigids.bed"
$KT compress --db "$W/biglab.tsv" -i "$W/bigids.bed" > "$W/biglab.kt" 2>/dev/null; mkdir -p "$W/bl"
$KT decompress --db "$W/biglab.tsv" --out-dir "$W/bl" -i "$W/biglab.kt" 2>/dev/null && [ "$(wc -l < "$W/bl/acrocentric.bed")" = 3 ] \
  && ok "2 MB features label: decompress --db" || no "2 MB features label: WRONG"

# round 18 findings, one case each
python3 -c "open('$W/L.bed', 'w').write('0\t0\t10\t' + 'z' * (42 << 17) + '\n')"
$KT compress --names drop < "$W/L.bed" > "$W/L.kt" 2>/dev/null
printf '%s\n' "$(head -c 500 /dev/zero | tr '\0' n)" > "$W/n500.names"
$KT --names "$W/n500.names" decompress -i "$W/L.kt" 2>/dev/null | cut -f1 | cmp -s - <(cut -f1 "$W/n500.names") \
  && ok "4.5 MB label with a 500-byte name: decompress" || no "4.5 MB label with a 500-byte name: WRONG"
printf '%s\n' "$(head -c 1000 /dev/zero | tr '\0' n)" > "$W/n1000.names"
must_fail "a name file name past 511 bytes" $KT --names "$W/n1000.names" decompress -i "$W/L.kt"

# pipe or slow producer equals file encode
cat "$D/s.bed" | $KT compress > "$W/pipe.kt" 2>/dev/null
(while IFS= read -r l; do printf '%s\n' "$l"; sleep 0.02; done < "$D/s.bed") | $KT compress > "$W/slow.kt" 2>/dev/null
cmp -s "$W/p.kt" "$W/pipe.kt" && cmp -s "$W/p.kt" "$W/slow.kt" && ok "compress from a pipe and a slow producer: byte-identical" || no "compress from a pipe: WRONG"
awk 'BEGIN{srand(3);for(i=0;i<20000;i++){p=0;for(j=0;j<5;j++){c=p+int(rand()*300)+1;printf "read%d\t%d\t%d\tL%d\n",i,p,c,int(rand()*6);p=c}}}' > "$W/pbig.bed"
$KT decompress --raw -i "$W/f.kt" 2>/dev/null | $KT compress --db "$D/db.features.tsv" -i - > "$W/pipedb.kt" 2>/dev/null
$KT decompress --raw -i "$W/f.kt" 2>/dev/null > "$W/merged.bed"; $KT compress --db "$D/db.features.tsv" -i "$W/merged.bed" > "$W/filedb.kt" 2>/dev/null
cmp -s "$W/pipedb.kt" "$W/filedb.kt" && ok "db compress from a pipe: byte-identical" || no "db compress from a pipe: WRONG"

# round 19 findings, one case each
printf 'chrA\t0\t100\tA\nchrA\t100\t200\tB\nchrB\t0\t100\tC\nchrB\t100\t200\tD\n' > "$W/r19.ref"
printf 'q0\t0\t0\t0\t+\tchrA\t200\t0\t100\t0\t100\t60\tcg:Z:100D\nq1\t200\t0\t200\t-\tchrB\t200\t0\t200\t200\t200\t60\tcg:Z:200M\n' > "$W/r19.paf"
printf 'q1\t0\t100\tD\nq1\t100\t200\tC\n' > "$W/r19.bed"
$KT compress --aln "$W/r19.paf" --ref "$W/r19.ref" < "$W/r19.bed" > "$W/r19p.kt" 2>/dev/null
$KT aln-pack -a "$W/r19.paf" --ref "$W/r19.ref" -o "$W/r19.kta" >/dev/null 2>&1
$KT compress --aln "$W/r19.kta" --ref "$W/r19.ref" < "$W/r19.bed" > "$W/r19k.kt" 2>/dev/null
cmp -s "$W/r19p.kt" "$W/r19k.kt" && $KT decompress --aln "$W/r19.kta" --ref "$W/r19.ref" -i "$W/r19p.kt" 2>/dev/null | cmp -s - "$W/r19.bed" \
  && ok "all-gap PAF record with a packed alignment: byte-identical" || no "all-gap PAF record: WRONG"
printf 'chrA\t100\t200\tB\nchrA\t0\t100\tA\n' > "$W/unsorted.ref"
must_fail "an unsorted reference" $KT compress --aln "$W/r19.paf" --ref "$W/unsorted.ref" < "$W/r19.bed"
sed 's/\tD$/\tE/' "$W/r19.ref" > "$W/other.ref"
must_fail "decompress with another reference" $KT decompress --aln "$W/r19.paf" --ref "$W/other.ref" -i "$W/r19p.kt"
sed 's/200M/100M100M/' "$W/r19.paf" > "$W/other.paf"
must_fail "decompress with another PAF" $KT decompress --aln "$W/other.paf" --ref "$W/r19.ref" -i "$W/r19p.kt"
sed 's/\t-\t/\t*\t/' "$W/r19.paf" > "$W/star.paf"
must_fail "a PAF strand of *" $KT compress --aln "$W/star.paf" --ref "$W/r19.ref" < "$W/r19.bed"

# round 20 findings, one case each
printf '\t0\t100\tA\n' > "$W/empty.ref"
must_fail "reference with an empty contig name" $KT compress --aln "$D/p.paf" --ref "$W/empty.ref" < "$D/s.bed"
N300=$(head -c 300 /dev/zero | tr '\0' N)
printf "$N300\t0\t100\tA\n$N300\t100\t200\tB\n" > "$W/n300.ref"
printf "x\t200\t0\t200\t+\t$N300\t200\t0\t200\t200\t200\t60\tcg:Z:200M\n" > "$W/n300.paf"
printf 'x\t0\t100\tA\nx\t100\t200\tB\n' > "$W/n300.bed"
$KT compress --aln "$W/n300.paf" --ref "$W/n300.ref" < "$W/n300.bed" > "$W/n300.kt" 2>/dev/null \
  && $KT decompress --aln "$W/n300.paf" --ref "$W/n300.ref" -i "$W/n300.kt" 2>/dev/null | cmp -s - "$W/n300.bed" \
  && ok "300-byte reference contig name: round-trip" || no "300-byte reference contig name: WRONG"
tac "$D/p.paf" > "$W/rev.paf"
$KT decompress --aln "$W/rev.paf" --ref "$D/ref.bed" -i "$W/d.kt" 2>/dev/null | cmp -s - "$D/s.bed" && ok "PAF in another record order: accepted" || no "PAF in another record order: WRONG"
must_fail "--names on aln-pack" $KT aln-pack --names "$W/x.kta" -a "$D/p.paf" --ref "$D/ref.bed"

# round 21 findings, one case each
printf 'feature_set\tchild\tparent\r\nacrocentric\tacro\tcategorized\r\nacrocentric\tnovel\tcategorized\r\n' > "$W/crlf.tsv"
$KT decompress -i "$W/sm.kt" --hier "$W/crlf.tsv" --smooth --set acrocentric 2>/dev/null | cmp -s - "$W/sm.1" && ok "CRLF hierarchy: same rows" || no "CRLF hierarchy: WRONG"
printf 'feature_set\tchild\tparent\nacrocentric\tacro\tcategorized\tx\n' > "$W/four.tsv"
must_fail "hierarchy row with four columns" $KT decompress -i "$W/sm.kt" --hier "$W/four.tsv" --smooth --set acrocentric
printf 'feature_set\tchild\tparent\nacrocentric\tacro\tcategorized\nacrocentric\tacro\tnovel\n' > "$W/dupc.tsv"
must_fail "hierarchy child listed twice" $KT decompress -i "$W/sm.kt" --hier "$W/dupc.tsv" --smooth --set acrocentric
printf 'feature_set\tchild\tparent\nacrocentric\tacro\troot\n' > "$W/root.tsv"
must_fail "hierarchy root that is not categorized" $KT decompress -i "$W/sm.kt" --hier "$W/root.tsv" --smooth --set acrocentric
must_fail "-r with a non-numeric range is a name" $KT query -i "$W/sm.kt" --hier "$W/sm.tsv" --smooth --set acrocentric -r r0:1e2-3e2
must_fail "query region with end before start" $KT query -i "$W/sm.kt" --hier "$W/sm.tsv" --smooth --set acrocentric -r r0:300-100
L309=$(head -c 300 /dev/zero | tr '\0' L)_specific
printf 's\t0\t10\t%s\n' "$L309" > "$W/l309.bed"; $KT compress < "$W/l309.bed" > "$W/l309.kt" 2>/dev/null
printf 'feature_set\tchild\tparent\nfs\t%s\tcategorized\n' "$L309" > "$W/l309.tsv"
[ "$($KT --v1 decompress -i "$W/l309.kt" --hier "$W/l309.tsv" --smooth --set fs 2>/dev/null | cut -f4 | wc -c)" = 301 ] && ok "V1 strip on a 309-byte label: stripped" || no "V1 strip on a 309-byte label: WRONG"

printf 'feature_set\tchild\tparent\nfs\ta\tcategorized\n' > "$W/miss.tsv"
printf 's\t0\t10\ta\ns\t10\t20\tb\ns\t20\t30\tnovel\n' | $KT compress > "$W/miss.kt" 2>/dev/null
must_fail "decompress: a track label missing from the hierarchy" $KT decompress -i "$W/miss.kt" --hier "$W/miss.tsv" --smooth --set fs
printf 'fs\tb\tcategorized\n' >> "$W/miss.tsv"
$KT decompress -i "$W/miss.kt" --hier "$W/miss.tsv" --smooth --set fs > /dev/null 2>&1 && ok "decompress: novel need not be in the hierarchy" || no "decompress: novel REFUSED"
# --chunk: KaryoScope 6186079 chunked smoothing
python3 -c "
import random; random.seed(5); p=0
with open('$W/ch.bed','w') as f:
    for i in range(6000):
        l=random.randint(1,40); f.write('s1\t%d\t%d\t%s\n'%(p,p+l,random.choice('a b c d G1 G2 a a c novel'.split()))); p+=l
"
printf 'feature_set\tchild\tparent\nfs\tG1\tcategorized\nfs\tG2\tcategorized\nfs\ta\tG1\nfs\tb\tG1\nfs\tc\tG2\nfs\td\tG2\n' > "$W/ch.tsv"
$KT compress < "$W/ch.bed" > "$W/ch.kst" 2>/dev/null
$KT decompress -i "$W/ch.kst" --hier "$W/ch.tsv" --smooth --set fs > "$W/ch.whole" 2>/dev/null
$KT decompress -i "$W/ch.kst" --hier "$W/ch.tsv" --smooth --set fs --chunk 64 2>/dev/null | cmp -s - "$W/ch.whole" && no "--chunk 64: SAME as whole" || ok "--chunk 64: chunks change the smoothing"
$KT decompress -i "$W/ch.kst" --hier "$W/ch.tsv" --smooth --set fs --chunk 5000 2>/dev/null | cmp -s - "$W/ch.whole" && ok "--chunk 5000 on 6000 runs: under 6250, smoothed whole" || no "--chunk 5000: WRONG"
$KT decompress -i "$W/ch.kst" --hier "$W/ch.tsv" --smooth 2>/dev/null | cmp -s - "$W/ch.whole" && ok "--smooth without --set: the one set that fits" || no "--smooth without --set: WRONG"
printf 'feature_set\tchild\tparent\nfs\ta\tcategorized\nzz\ta\tcategorized\n' > "$W/amb.tsv"
printf 's\t0\t10\ta\ns\t10\t20\tnovel\n' | $KT compress > "$W/amb.kt" 2>/dev/null
must_fail "--smooth without --set, two sets fit" $KT decompress -i "$W/amb.kt" --hier "$W/amb.tsv" --smooth
$KT compress < "$D/ls.bed" > "$W/nofit.kt" 2>/dev/null
must_fail "--smooth without --set, no set fits" $KT decompress -i "$W/nofit.kt" --hier "$W/ch.tsv" --smooth
must_fail "--chunk 0" $KT decompress -i "$W/ch.kst" --hier "$W/ch.tsv" --smooth --set fs --chunk 0
must_fail "--chunk on decompress" $KT decompress -i "$W/ch.kst" --chunk 64

# round 22 findings, one case each
mkdir -p "$W/keepdir"; printf 'PRECIOUS\n' > "$W/keepdir/acrocentric.bed"
$KT decompress --db "$D/db.features.tsv" --out-dir "$W/keepdir" -i "$W/p.kt" >/dev/null 2>&1 && { echo "decompress --db of a plain track: EXIT 0"; fail=1; }
[ "$(cat "$W/keepdir/acrocentric.bed")" = PRECIOUS ] && ok "refused decompress --db: old BEDs intact" || no "refused decompress --db: OLD BEDS DESTROYED"
printf 'featureID\ts1\ts2\n1\tx\tp\n2\ty\tp\n' > "$W/two.tsv"
printf 's\t0\t10\tx\ns\t10\t20\ty\n' > "$W/n1.s1.bed"; printf 's\t0\t10\tp\ns\t10\t20\tp\n' > "$W/n1.s2.bed"
must_fail "touching same-label runs in one set" $KT compress --db "$W/two.tsv" -i "$W/n1.s1.bed" -i "$W/n1.s2.bed"
printf 's\t0\t10\t+1\ns\t10\t20\t2\n' > "$W/plus.bed"
must_fail "featureID +1 in a merged track" $KT compress --db "$W/two.tsv" -i "$W/plus.bed"
L201=$(head -c 201 /dev/zero | tr '\0' L); { printf 'featureID'; for i in $(seq 33); do printf '\ts%d' $i; done; printf '\n1'; for i in $(seq 33); do printf '\t%s%d' "$L201" $i; done; printf '\n'; } > "$W/t33.tsv"
f33=""; for i in $(seq 33); do printf 's\t0\t10\t%s%d\n' "$L201" $i > "$W/t33.$i.bed"; f33="$f33 -i $W/t33.$i.bed"; done
$KT compress --db "$W/t33.tsv" $f33 > "$W/t33.kt" 2>/dev/null && ok "33 sets with 201-byte labels: encoded" || no "33 sets with 201-byte labels: REFUSED"
L300=$(head -c 300 /dev/zero | tr '\0' L); printf "featureID\ts1\ts2\n1\t$L300\tp\n" > "$W/l300.tsv"
printf "s\t0\t10\t$L300\n" > "$W/l300.s1.bed"; printf 's\t0\t10\tp\n' > "$W/l300.s2.bed"; mkdir -p "$W/l300o"
$KT compress --db "$W/l300.tsv" -i "$W/l300.s1.bed" -i "$W/l300.s2.bed" > "$W/l300.kt" 2>/dev/null \
  && $KT decompress --db "$W/l300.tsv" --out-dir "$W/l300o" -i "$W/l300.kt" 2>/dev/null && cmp -s "$W/l300o/s1.bed" "$W/l300.s1.bed" \
  && ok "300-byte label in a set: round-trip" || no "300-byte label in a set: WRONG"
head -c -1 "$D/db.gene.bed" > "$W/nonl.gene.bed"
must_fail "a set BED without a final newline" $KT compress --db "$D/db.features.tsv" -i "$D/db.acrocentric.bed" -i "$D/db.chromosome.bed" -i "$W/nonl.gene.bed" -i "$D/db.region.bed" -i "$D/db.repeat.bed" -i "$D/db.subtelomeric.bed"

# round 23 findings, one case each
printf '\001\t0\t10\t1\n\001\t10\t20\t2\nb\t0\t10\t1\n' > "$W/c3.bed"
$KT compress --db "$W/two.tsv" -i "$W/c3.bed" > "$W/c3.kt" 2>/dev/null && $KT decompress --raw -i "$W/c3.kt" 2>/dev/null | cmp -s - "$W/c3.bed" \
  && ok "first sequence named 0x01 under --db: round-trip" || no "first sequence named 0x01 under --db: WRONG"
LD="$W"; while [ ${#LD} -lt 3970 ]; do LD="$LD/$(printf 'd%.0s' $(seq 1 99))"; done
LD="$LD/$(printf 'e%.0s' $(seq 1 $((4069 - ${#LD}))))"; mkdir -p "$LD"
must_fail "decompress --db --out-dir with a 4070-byte path" $KT decompress --db "$D/db.features.tsv" --out-dir "$LD" -i "$W/f.kt"
$KT --max-gap 18446744073709551615 decompress -i "$W/sm.kt" --hier "$W/sm.tsv" --smooth --set acrocentric 2>/dev/null | cmp -s - "$W/sm.1" && ok "MAXGAP of 2^64-1: same as default" || no "MAXGAP of 2^64-1: WRONG"
must_fail "region with END before START" $KT query -i "$W/p.kt" -r q1:150-20

# round 24 findings, one case each
rc=0; (ulimit -f 1; $KT compress < "$W/pbig.bed" > "$W/lim.kt" 2> "$W/lim.err") || rc=$?
[ $rc = 1 ] && ok "file-size limit: exit 1, not a signal" || no "file-size limit: rc=$rc"
sz=""; for e in high medium low; do sz="$sz $($KT -e $e compress < "$D/ls.bed" 2>/dev/null | wc -c)"; done
[ "$(echo $sz | tr ' ' '\n' | sort -n | tr '\n' ' ')" = "$(echo $sz) " ] && ok "effort high <= medium <= low:$sz B" || no "effort sizes out of order:$sz B"
# -e high: zpaq on both frames
awk 'BEGIN { srand(5); for (c = 0; c < 6; c++) { g = c % 4; nr = 20000 + int(rand() * 20000); p = 0; for (r = 0; r < nr; r++) { l = 1 + int(rand() * (g ? 9 * g : 37)); u = rand(); printf "ctg%d\t%d\t%d\t%s\n", c, p, p + l, u < 0.5 + g / 10 ? "chr" (1 + g) : "novel" int(rand() * 3); p += l } } }' > "$W/zt.bed"
$KT -t 4 -e high compress -i "$W/zt.bed" -o "$W/zt.kt" 2> "$W/zt.err"
grep -q "lengths zpaq-64MB, labels zpaq-64MB" "$W/zt.err" && ok "high: zpaq on both frames of six contigs" || no "high: zpaq not picked: $(grep layout "$W/zt.err")"
$KT -t 4 decompress -i "$W/zt.kt" 2>/dev/null | cmp -s - "$W/zt.bed" && ok "zpaq frames: round-trip" || no "zpaq frames: WRONG"
# -e high: interleaved zpaq model, short reads
awk 'BEGIN { srand(7); for (r = 0; r < 60000; r++) { k = rand() < 0.4 ? 1 : 2 + int(rand() * rand() * 12); p = 0; g = int(rand() * 4); for (j = 0; j < k; j++) { l = j == k - 1 ? 120 - p : 1 + int(rand() * (120 - p - (k - j - 1)) / (k - j)); printf "r%d\t%d\t%d\t%s\n", r, p, p + l, "chr" (1 + (g + j) % 5); p += l } } }' > "$W/zi.bed"
$KT -t 4 -e high compress -i "$W/zi.bed" -o "$W/zi.kt" 2> "$W/zi.err"
grep -q "lengths zpaq-64MB+ilv" "$W/zi.err" && ok "high: the interleaved zpaq model on short reads" || no "high: interleaved model not picked: $(grep layout "$W/zi.err")"
$KT -t 4 decompress -i "$W/zi.kt" 2>/dev/null | cmp -s - "$W/zi.bed" && ok "interleaved zpaq model: round-trip" || no "interleaved zpaq model: WRONG"
must_fail "--tmp-dir that does not exist" $KT --tmp-dir /nonexistent/dir compress < "$D/s.bed"
$KT --tmp-dir "$W" compress < "$D/s.bed" 2>/dev/null | cmp -s - "$W/p.kt" && ok "--tmp-dir honoured: same bytes" || no "--tmp-dir honoured: WRONG"

# round 25 findings, one case each
# environment variables are ignored
ok=1
for e in 1 3 16 100000; do
  OMP_NUM_THREADS=$e $KT compress < "$D/db.rall.bed" 2>/dev/null | cmp -s - "$W/omp.ref" || ok=0
  OMP_NUM_THREADS=$e $KT --threads 2 compress < "$D/db.rall.bed" 2>/dev/null | cmp -s - "$W/omp.ref" || ok=0
done
TMPDIR=/nonexistent $KT compress < "$D/s.bed" 2>/dev/null | cmp -s - "$W/p.kt" || ok=0
[ $ok = 1 ] && ok "OMP_NUM_THREADS, TMPDIR: ignored" || no "OMP_NUM_THREADS or TMPDIR: READ"
must_fail "--block 4abc" $KT compress --block 4abc < "$D/s.bed"
must_fail "-t 4x" $KT -t 4x compress < "$D/s.bed"

# round 26 findings, one case each
cp "$W/A.txt" "$W/A2.txt"; $KT --names "$W/A2.txt" index -i "$W/nx.kt" >/dev/null 2>&1; cp "$W/B.txt" "$W/A2.tmp"; mv -f "$W/A2.tmp" "$W/A2.txt"
[ "$($KT --names "$W/A2.txt" query -i "$W/nx.kt" -r r0 2>/dev/null | cut -f4)" = D ] && ok "index built before a name file swap: ignored" || no "index built before a name file swap: WRONG"
mkdir -p "$W/good"; $KT decompress --db "$D/db.features.tsv" --out-dir "$W/good" -i "$W/f.kt" >/dev/null 2>&1; m1=$(cat "$W/good"/*.bed | md5sum)
cp "$W/f.kt" "$W/fbad.kt"; sz=$(stat -c%s "$W/fbad.kt"); printf '\xff' | dd of="$W/fbad.kt" bs=1 seek=$((sz-100)) conv=notrunc status=none
$KT decompress --db "$D/db.features.tsv" --out-dir "$W/good" -i "$W/fbad.kt" >/dev/null 2>&1 && { echo "decompress of a damaged db track: EXIT 0"; fail=1; }
[ "$(cat "$W/good"/*.bed | md5sum)" = "$m1" ] && [ -z "$(ls "$W/good" | grep tmp)" ] && ok "failed decompress --db: previous output intact, no temporaries" || no "failed decompress --db: OUTPUT DAMAGED"
mkdir -p "$W/oo"; echo old > "$W/oo/k.kt"; $KT compress -i "$W/nonexist" -o "$W/oo/k.kt" >/dev/null 2>&1 || :
$KT decompress -i "$W/fbad.kt" -o "$W/oo/k.kt" >/dev/null 2>&1 || :; $KT compress -i "$W/nonexist" -o "$W/oo/n.kt" >/dev/null 2>&1 || :
[ "$(cat "$W/oo/k.kt")" = old ] && [ "$(ls "$W/oo")" = k.kt ] && ok "failed -o: previous file intact, no temporaries" || no "failed -o: OUTPUT DAMAGED"
$KT decompress -i "$W/nx.kt" -o "$W/oo/k.kt" >/dev/null 2>&1 || :; [ "$(ls "$W/oo")" = k.kt ] && [ -s "$W/oo/k.kt" ] && ok "-o over an old file: replaced" || no "-o over an old file: WRONG"

# round 27 findings, one case each
printf 'r0\nr1\nr2\nr3\n' > "$W/pn.txt"
$KT --names <(cat "$W/pn.txt") decompress -i "$W/nx.kt" 2>/dev/null | cmp -s - "$W/nx.bed" && ok "name file from a pipe: restored" || no "name file from a pipe: WRONG"
printf '@r0\nACGT\n+\nIIII\n@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIIII\n@r3\nACGT\n+\nIIII\n' > "$W/pn.fq"
$KT --names <(cat "$W/pn.fq") decompress -i "$W/nx.kt" 2>/dev/null | cmp -s - "$W/nx.bed" && ok "FASTQ from a pipe: restored" || no "FASTQ from a pipe: WRONG"
mkdir -p "$W/mix"; for f in acrocentric chromosome gene region repeat subtelomeric; do echo "OLD_$f" > "$W/mix/$f.bed"; done
mkdir -p "$W/mix/gene.bed.blocked" && mv "$W/mix/gene.bed" "$W/mix/gene.keep" && mkdir "$W/mix/gene.bed"
must_fail "decompress --db over a directory in --out" $KT decompress --db "$D/db.features.tsv" --out-dir "$W/mix" -i "$W/f.kt"
[ "$(cat "$W/mix/acrocentric.bed")" = OLD_acrocentric ] && [ -z "$(ls "$W/mix" | grep -E 'tmp|old')" ] \
  && ok "refused db decompress: all old files kept, no leftovers" || no "refused db decompress: HALF UPDATED"

# exclusive flag pairs and range limits
must_fail "--block junk" $KT compress --block xyz < "$D/s.bed"
must_fail "--threads 1025" $KT -t 1025 compress < "$D/s.bed"
$KT -t 1024 compress < "$D/s.bed" > /dev/null 2>&1 && ok "the documented limits: accepted" || no "the documented limits: REFUSED"

# flags work before or after the command
$KT --effort low compress < "$D/s.bed" > "$W/side1.kt" 2>/dev/null
$KT compress --effort low < "$D/s.bed" > "$W/side2.kt" 2>/dev/null
cmp -s "$W/side1.kt" "$W/side2.kt" && ok "a flag before the command: same bytes" || no "a flag before the command: DIFFERENT"

# a non-KST file: named as such
printf 'KTRACK02rest' > "$W/old.kst"
$KT decompress -i "$W/old.kst" 2>&1 | grep -q 'not "KST00001"' && ok "a foreign magic: named" || no "a foreign magic: WRONG MESSAGE"

must_fail "unknown flag" $KT --bogus compress
$KT --help > "$W/help.txt" 2>/dev/null && [ -s "$W/help.txt" ] && ok "--help: on stdout, exit 0" || no "--help: WRONG"
must_fail "compress --out" $KT compress --out-dir "$W" < "$D/s.bed"
must_fail "compress --names ordinal" $KT compress --names ordinal < "$D/s.bed"
must_fail "decompress two tracks" $KT decompress -i "$W/side1.kt" -i "$W/side2.kt"

# nothing from htscodecs reaches stderr
$KT compress < "$D/s.bed" 2>&1 >/dev/null | grep -qi "htscodecs" && { echo "compress stderr: HTSCODECS COMPLAINED"; fail=1; } || echo "compress stderr: clean"

# library API: same bytes as the CLI
if [ -x "$LT" ]; then
  lok=1
  for f in s.bed ls.bed db.rall.bed db.chromosome.bed ord.bed; do
    $KT compress < "$D/$f" > "$W/lc.kt" 2>/dev/null
    $LT "$D/$f" "$W/ll.kt" 2>/dev/null
    cmp -s "$W/lc.kt" "$W/ll.kt" || { lok=0; echo "  library differs on $f"; }
  done
  $KT --block 262144 compress < "$D/db.rall.bed" > "$W/lc.kt" 2>/dev/null
  $LT "$D/db.rall.bed" "$W/ll.kt" 262144 2>/dev/null
  cmp -s "$W/lc.kt" "$W/ll.kt" || { lok=0; echo "  library differs on --block 262144"; }
  $KT --block 4 compress < "$D/db.rall.bed" > "$W/lc.kt" 2>/dev/null
  $LT "$D/db.rall.bed" "$W/ll.kt" 4 2>/dev/null
  $LT -twice "$W/lt1.kt" "$W/lt2.kt" 2>/dev/null; printf 't1\t0\t7\tCCC\nt1\t7\t12\tAAA\n' | $KT compress 2>/dev/null | cmp -s - "$W/lt2.kt" \
    && ok "library API: a second writer starts clean" || no "library API: a second writer keeps the first's state"
  cmp -s "$W/lc.kt" "$W/ll.kt" || { lok=0; echo "  library differs on --block 4"; }
  [ $lok = 1 ] && ok "library API: same bytes as the CLI" || no "library API: DIFFERS from the CLI"
  # library refuses all 12 contract violations
  aok=1
  for c in 1 2 3 4 5 6 7 8 9 10 11 12; do
    rc=0; $LT -abuse $c >/dev/null 2>"$W/ab.err" || rc=$?
    [ $rc = 1 ] && grep -q '^karyozip: ' "$W/ab.err" || { aok=0; echo "  library abuse case $c: exit $rc, $(head -c 100 "$W/ab.err")"; }
  done
  [ $aok = 1 ] && ok "library API: every contract violation refused" || no "library API: A CONTRACT VIOLATION WAS ACCEPTED"
else
  echo "library API: FAILED, $LT not built"; fail=1
fi

# label cache collisions, prefix names, chunk boundaries
printf 'r1\t0\t10\tabcXYZ\nr1\t10\t20\tabqXYZ\nr1\t20\t30\tabcXYZ\n' > "$W/coll.bed"
$KT compress < "$W/coll.bed" > "$W/coll.kt" 2>/dev/null
$KT decompress -i "$W/coll.kt" 2>/dev/null | cmp -s - "$W/coll.bed" && ok "label cache: colliding keys round-trip" || no "label cache: COLLIDING KEYS WRONG"
python3 -c "
p=0
f=open('$W/many.bed','w')
for i in range(60000):
    f.write('r1\t%d\t%d\ta%04dz\n' % (p, p+10, i % 20000)); p += 10
"
$KT compress < "$W/many.bed" > "$W/many.kt" 2>/dev/null
$KT decompress -i "$W/many.kt" 2>/dev/null | cmp -s - "$W/many.bed" && ok "label cache: 20000 labels round-trip" || no "label cache: 20000 labels WRONG"
printf 'r1\t0\t10\tA\nr11\t0\t10\tB\nr111\t0\t10\tC\n' > "$W/pre.bed"
$KT compress < "$W/pre.bed" > "$W/pre.kt" 2>/dev/null
$KT decompress -i "$W/pre.kt" 2>/dev/null | cmp -s - "$W/pre.bed" && ok "prefix sequence names: round-trip" || no "prefix sequence names: WRONG"
python3 -c "
f=open('$W/chunk.bed','w')
for s in range(200000): f.write('seq%06d\t0\t10\tA\n' % s)
"
$KT -t 8 compress < "$W/chunk.bed" > "$W/chunk.kt" 2>/dev/null
$KT decompress -i "$W/chunk.kt" 2>/dev/null | cmp -s - "$W/chunk.bed" && ok "200000 one-row sequences at 8 threads: round-trip" || no "200000 one-row sequences: WRONG"
# a full temp filesystem must be fatal
mkdir -p "$W/tinytmp"
awk 'BEGIN{p=0; for(k=0;k<400000;k++){printf "r1\t%d\t%d\tchr%d\n", p, p+20, k%22; p+=20}}' > "$W/stage.bed"
# 256 KB cap, below the staging payload
( ulimit -f 512; trap "" XFSZ; $KT --tmp-dir "$W/tinytmp" compress < "$W/stage.bed" >/dev/null 2>&1 ) && { echo "a full staging filesystem: ACCEPTED"; fail=1; } \
  || echo "a full staging filesystem: refused"
[ -z "$(ls -A "$W/tinytmp" 2>/dev/null)" ] && ok "a full staging filesystem: no temp files left" || no "a full staging filesystem: TEMP FILES LEFT"

# whole-read blocks: layouts, queries, --db, --aln, --smooth
awk 'BEGIN { for (r = 0; r < 3000; r++) { n = 1 + r % 4; p = 0; for (i = 0; i < n; i++) { l = i + 1 == n ? 151 - p : 20 + (r * 7 + i * 13) % 40; printf "q%d\t%d\t%d\t%s\n", r, p, p + l, (r + i) % 5 ? "chr" (1 + (r + i) % 3) : "novel"; p += l } } }' > "$W/vb.bed"
for n in 1 7 500; do
  $KT compress --read-blocks $n -i "$W/vb.bed" -o "$W/vb.kt" 2> "$W/vb.err"
  $KT decompress -i "$W/vb.kt" 2>/dev/null | cmp -s - "$W/vb.bed" && ok "--read-blocks $n: byte-identical, $(grep -o 'layout: [a-z]*' "$W/vb.err"), $(grep -o 'blocks=[0-9]*' "$W/vb.err")" || no "--read-blocks $n: WRONG"
  for ix in 0 1; do [ $ix = 1 ] && $KT index -i "$W/vb.kt" 2>/dev/null; q=0
    for r in q0 q1 q499 q500 q1777 q2999; do $KT query -i "$W/vb.kt" -r $r 2>/dev/null | cmp -s - <(awk -v r=$r '$1 == r' "$W/vb.bed") && q=$((q+1)); done
    $KT query -i "$W/vb.kt" -r q1777:30-90 2>/dev/null | cmp -s - <(awk '$1 == "q1777" && $3 > 30 && $2 < 90' "$W/vb.bed") && q=$((q+1))
    [ $q = 7 ] && ok "--read-blocks $n query, index $ix: 7/7" || no "--read-blocks $n query, index $ix: only $q/7"; done
  rm -f "$W/vb.kt.kzi"
done
grep -q "layout: plain" "$W/vb.err" && no "--read-blocks 500: PLAIN, the layout search did not run" || ok "--read-blocks 500: a one-block layout on every block"
$KT compress --read-blocks 3 --db "$D/db.features.tsv" $in -o "$W/vbd.kt" 2>/dev/null; rm -rf "$W/vbo"; mkdir -p "$W/vbo"
$KT decompress --db "$D/db.features.tsv" --out-dir "$W/vbo" -i "$W/vbd.kt" 2>/dev/null
n=0; for f in $sets; do cmp -s "$W/vbo/$f.bed" "$D/db.$f.bed" && n=$((n+1)); done
[ "$n" = 6 ] && ok "--read-blocks 3 --db: 6/6 identical" || no "--read-blocks 3 --db: only $n/6 identical"
$KT compress --read-blocks 3 --aln "$D/p.paf" --ref "$D/ref.bed" -i "$D/s.bed" -o "$W/vba.kt" 2>/dev/null
$KT decompress --aln "$D/p.paf" --ref "$D/ref.bed" -i "$W/vba.kt" 2>/dev/null | cmp -s - "$D/s.bed" && ok "--read-blocks 3 --aln: byte-identical" || no "--read-blocks 3 --aln: WRONG"
{ printf 'feature_set\tchild\tparent\n'; for l in chr1 chr2 chr3 novel; do printf 'fs\t%s\tcategorized\n' $l; done; } > "$W/vbh.tsv"
$KT compress -i "$W/vb.bed" -o "$W/vb1.kt" 2>/dev/null; $KT compress --read-blocks 7 -i "$W/vb.bed" -o "$W/vb.kt" 2>/dev/null
for t in 1 4; do $KT decompress -t $t -i "$W/vb.kt" --smooth --hier "$W/vbh.tsv" 2>/dev/null | cmp -s - <($KT decompress -i "$W/vb1.kt" --smooth --hier "$W/vbh.tsv" 2>/dev/null) && ok "--read-blocks 7 --smooth -t $t: same as one block" || no "--read-blocks 7 --smooth -t $t: DIFFERS"; done

# decode ahead: same rows, same stop
awk 'BEGIN { x = 7; for (r = 0; r < 120000; r++) { n = 1 + r % 4; p = 0; for (i = 0; i < n; i++) { x = (x * 1103515245 + 12345) % 2147483648; l = i + 1 == n ? 151 - p : 1 + x % 37; printf "q%d\t%d\t%d\t%s\n", r, p, p + l, x % 7 ? "chr" (1 + x % 3) : "novel"; p += l } } }' > "$W/pf.bed"
for bl in "--read-blocks 70000" "--block 70000"; do
  $KT compress $bl -i "$W/pf.bed" -o "$W/pf.kt" 2>/dev/null
  $KT decompress -t 4 -i "$W/pf.kt" 2>/dev/null | cmp -s - "$W/pf.bed" && ok "decode ahead [$bl]: byte-identical" || no "decode ahead [$bl]: WRONG"
  cp "$W/pf.kt" "$W/pfb.kt"; sz=$(stat -c%s "$W/pfb.kt"); printf '\xff' | dd of="$W/pfb.kt" bs=1 seek=$((sz - 300)) conv=notrunc status=none
  a=0; $KT decompress -t 1 -i "$W/pfb.kt" > "$W/pf1.out" 2> "$W/pf1.err" || a=$?
  b=0; $KT decompress -t 4 -i "$W/pfb.kt" > "$W/pf4.out" 2> "$W/pf4.err" || b=$?
  cmp -s "$W/pf1.out" "$W/pf4.out" && cmp -s "$W/pf1.err" "$W/pf4.err" && [ $a = $b ] && [ $a = 1 ] && ok "decode ahead [$bl], damaged: stops at the same byte, $(wc -c < "$W/pf1.out") B" || no "decode ahead [$bl], damaged: DIFFERS from one thread"
done
if [ "$SAN" != "$KT" ]; then ASAN_OPTIONS=detect_leaks=1 $SAN decompress -t 4 -i "$W/pf.kt" 2> "$W/pfs.err" | cmp -s - "$W/pf.bed" && ! grep -qE "Sanitizer|runtime error" "$W/pfs.err" && ok "decode ahead: sanitizer clean" || no "decode ahead: SANITIZER REPORT"; fi

# pin the stored bytes to codec versions
want=55b6e9eadfe851a4374e094e045ee8ecb4ffac9a2fdd83de126ff076724e700b
got=$($KT --block 4 compress < "$D/db.rall.bed" 2>/dev/null | sha256sum | cut -d' ' -f1)
if [ "$got" = "$want" ]; then echo "codec versions: stored bytes unchanged"
else
  echo "codec versions: STORED BYTES MOVED"; fail=1
  echo "  want $want"; echo "  got  $got"; $KT --version | sed 's/^/  /'
fi

if bash "$D/fuzz.sh" "$KT" "$SAN"; then echo "fuzz: clean"; else echo "fuzz: FOUND SOMETHING"; fail=1; fi

exit $fail
