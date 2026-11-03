#!/bin/bash
# aln.expected cases, all computed by hand:
#   q1  forward, exact
#   q2  minus strand, maps to tp+n-1-pos
#   q3  50M20D50M, the deletion shifts the target
#   q4  50M20I50M, insertion bases have no target
#   q5  no alignment at all
#   run.sh BINARY LIBTEST [SANITIZER_BINARY]
# LIBTEST: C library test, compared to BINARY
# Every check prints one line, failures prefixed "FAIL: "; the last line is "run.sh: N passed, M failed, K skipped".
# CASE=REGEX: every command still runs (later checks reuse earlier files), but only checks whose description
#   matches are printed and counted; fuzz.sh runs only when CASE is unset or matches "fuzz".
# TEST_JOBS, KZ_TEST_THREADS, TEST_LIGHT, TIMEOUT_SCALE: see lib.sh. $KT and $SAN carry -t $KZ_TEST_THREADS;
#   $KZ and $KZSAN are the bare binaries, for the cases that test the default thread count.
set -eu
. "$(dirname "$0")/lib.sh"
need python3 zstd gzip awk sha256sum md5sum timeout realpath mktemp cmp diff dd stat tac cut sort uniq sed tr wc head tail seq xargs mkfifo
need_gnu
KZ=$(realpath -m "$1"); LT=$(realpath -m "$2"); KZSAN=$(realpath -m "${3:-$1}")
KT="$KZ -t $KZ_TEST_THREADS"; SAN="$KZSAN -t $KZ_TEST_THREADS"   # a later -t wins, in the C and the Rust
D=$(dirname "$0")
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
src=$D; D=$W/fx; mkdir -p "$D"; cp "$src"/fixtures.txt "$src"/*.py "$src"/*.sh "$D"/   # the fixtures, unpacked
awk '/^@@ /{f=d"/"$2; next} {print > f}' d="$D" "$D/fixtures.txt"
tscale "$KZSAN"
ok() { tpass "$1"; }
no() { tfail "$1"; }
must_fail() {  # DESC CMD...: expect exit 1, message
  local what=$1 rc=0; shift
  "$@" >/dev/null 2>"$W/mf.err" || rc=$?
  if [ $rc = 1 ] && grep -q '^karyozip: ' "$W/mf.err"; then ok "$what: refused ($(grep -m1 '^karyozip: ' "$W/mf.err" | cut -c11-90))"
  else no "$what: WRONG, exit $rc: $(head -c 200 "$W/mf.err")"; fi
}

$KT compress --aln "$D/p.paf" --ref "$D/ref.bed" < "$D/s.bed" > "$W/d.kt" 2>/dev/null
$KT decompress --raw -i "$W/d.kt" 2>/dev/null | cut -f2,3,4 > "$W/got"
cut -f2,3,4 "$D/aln.expected" > "$W/want"
cmp -s "$W/got" "$W/want" && ok "aln delta: correct" || { no "aln delta: WRONG"; diff "$W/got" "$W/want" | head; }

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
# several readers at once share no state: TEST_JOBS of them, at most 8, sharing KZ_TEST_THREADS
NR=$TEST_JOBS; [ "$NR" -gt 8 ] && NR=8
if stress "$NR readers at once"; then
$KT decompress -i "$W/rall.kst" > "$W/conc.want" 2>/dev/null
seq=$(head -1 "$D/db.rall.bed" | cut -f1)
$KT query -i "$W/rall.kst" -r "$seq" > "$W/conc.rwant" 2>/dev/null
KR="$KZ -t $(share $NR)"
for i in $(seq 1 $NR); do
  ( $KR decompress -i "$W/rall.kst" > "$W/c.d.$i" 2>/dev/null
    $KR query -i "$W/rall.kst" -r "$seq" > "$W/c.r.$i" 2>/dev/null
    $KR decompress --aln "$D/db.self.paf" --ref "$W/rall.kst" -i "$W/rall.aln.kst" > "$W/c.a.$i" 2>/dev/null ) &
done
wait
bad=0
for i in $(seq 1 $NR); do
  cmp -s "$W/c.d.$i" "$W/conc.want" || bad=1
  cmp -s "$W/c.r.$i" "$W/conc.rwant" || bad=1
  cmp -s "$W/c.a.$i" "$W/rall.want" || bad=1
done
[ $bad = 0 ] && ok "$NR readers at once: same bytes" || no "$NR readers at once: WRONG"
fi
# review cases of 2026-09-07, one each
python3 - "$W" <<'PY'
import sys; W = sys.argv[1]; L = "L" * 2100000
open(W + "/rv.tsv", "w").write("featureID\tA\tB\n1\tx\ty\n2\t%s\tz\n" % L)
open(W + "/rv1.bed", "w").write("s1\t0\t10\tx\ns1\t10\t20\t%s\n" % L); open(W + "/rv2.bed", "w").write("s1\t0\t10\ty\ns1\t10\t20\tz\n")
open(W + "/rv.big", "w").write("1125899906842625\t0\t10\tx\n"); open(W + "/rv.names", "w").write("q1\n\nq2\n")
PY
$KT --db "$W/rv.tsv" compress -i "$W/rv1.bed" -i "$W/rv2.bed" > "$W/rv.kt" 2>/dev/null
# 2100012 = "s1\t10\t20\t" (9 B) + the 2100000-byte label of set A + "\t" + set B's "z" + "\n"
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
  $KT compress $f < "$D/s.bed" > "$W/f.kt" 2>/dev/null || { no "$f: ENCODE FAILED"; continue; }
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
stress "decompress --db --aln, 4 threads: same as plain" && { $KT --db "$D/db.features.tsv" --aln "$D/db.self.paf" --ref "$D/db.rall.bed" -t 4 decompress -i "$W/da.kt" --hier "$W/inf.tsv" --smooth --set acrocentric 2>/dev/null | cmp -s - "$W/inf.want" && ok "decompress --db --aln, 4 threads: same as plain" || no "decompress --db --aln: DIFFERS"; }
$KT --aln "$D/p.paf" --ref "$D/ref.bed" decompress -i "$W/d.kt" --hier "$W/inf.tsv" --smooth --set acrocentric 2>/dev/null | cmp -s - "$W/inf.wants" && ok "decompress --aln: same as plain" || no "decompress --aln: DIFFERS"
sq=$(head -1 "$D/db.acrocentric.bed" | cut -f1)
$KT query -i "$W/pa.kt" -r "$sq":0-100000 > "$W/inf.want" 2>/dev/null
$KT --db "$D/db.features.tsv" query -i "$W/f.kt" -r "$sq":0-100000 2>/dev/null > "$W/inf.got"
[ "$(awk -F'\t' '{print NF}' "$W/inf.got" | sort -u)" = 9 ] && ok "decompress --db: 6 label columns" || no "decompress --db: wrong column count"
awk -F'\t' 'NR>1 && $1==p && $2==e && $4==l {e=$3; next} NR>1 {print p"\t"s"\t"e"\t"l} {p=$1; s=$2; e=$3; l=$4} END{print p"\t"s"\t"e"\t"l}' <(cut -f1-4 "$W/inf.got") | cmp -s - "$W/inf.want" && ok "decompress --db: acrocentric column, merged, same as plain" || no "decompress --db: DIFFERS"
$KT --db "$D/db.features.tsv" --aln "$D/db.self.paf" --ref "$D/db.rall.bed" query -i "$W/da.kt" -r "$sq":0-100000 2>/dev/null | cmp -s - "$W/inf.got" && ok "decompress --db --aln: same as --db" || no "decompress --db --aln: DIFFERS"
seqcmp $KT --aln "$D/p.paf" --ref "$D/ref.bed" query -i "$W/d.kt" -r q1 -- $KT query -i "$W/p.kt" -r q1 && ok "decompress --aln: same as plain" || no "decompress --aln: DIFFERS"
# --set: one set equals its own track
n=0; for f in $sets; do $KT --db "$D/db.features.tsv" --set $f decompress -i "$W/f.kt" 2>/dev/null | cmp -s - "$D/db.$f.bed" && n=$((n+1)); done
[ "$n" = 6 ] && ok "decompress --set: 6/6 equal to the set's track" || no "decompress --set: only $n/6"
$KT --db "$D/db.features.tsv" --aln "$D/db.self.paf" --ref "$D/db.rall.bed" --set gene decompress -i "$W/da.kt" 2>/dev/null | cmp -s - "$D/db.gene.bed" && ok "decompress --set on db+aln: same" || no "decompress --set on db+aln: DIFFERS"
n=0; for r in "-r $sq" "-r $sq:100-101" "-r $sq:0-100000"; do seqcmp $KT --db "$D/db.features.tsv" --set acrocentric query $r -i "$W/f.kt" -- $KT query $r -i "$W/pa.kt" && [ -n "$($KT query $r -i "$W/pa.kt" 2>/dev/null)" ] && n=$((n+1)); done
[ "$n" = 3 ] && ok "query --set: 3/3 forms equal the set's track" || no "query --set: only $n/3"
$KT --db "$D/db.features.tsv" --set gene,acrocentric decompress -i "$W/f.kt" 2>/dev/null > "$W/ga.bed"
[ "$(awk -F'\t' '{print NF}' "$W/ga.bed" | sort -u)" = 5 ] && awk -F'\t' 'NR>1 && $1==p && $2==e && $5==l {e=$3; next} NR>1 {print p"\t"s"\t"e"\t"l} {p=$1; s=$2; e=$3; l=$5} END{print p"\t"s"\t"e"\t"l}' "$W/ga.bed" | cmp -s - "$D/db.acrocentric.bed" && ok "decompress --set gene,acrocentric: two columns, merged jointly" || no "decompress --set gene,acrocentric: WRONG"
seqcmp $KT --db "$D/db.features.tsv" decompress -i "$W/f.kt" -- $KT --db "$D/db.features.tsv" --set all decompress -i "$W/f.kt" && ok "decompress --db on stdout: all sets, same as --set all" || no "decompress --db on stdout: DIFFERS from --set all"
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
# the --db tables are indexed by featureID: about 0.6 GB to compress and 1 GB to decompress this one row (measured, C)
stress "featureID 2^26-1, about 1 GB of --db tables: round trip" && { $KT compress --db "$W/big.tsv" -i "$W/big.bed" > "$W/big.kt" 2>/dev/null && $KT decompress --db "$W/big.tsv" --set acrocentric -i "$W/big.kt" 2>/dev/null | grep -q acro \
  && ok "featureID 2^26-1: round trip" || no "featureID 2^26-1: WRONG"; }
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
seqcmp $KT --db "$D/db.features.tsv" decompress -i "$W/hv.kt" --hier "$W/hv.tsv" --smooth --set repeat -- $KT decompress -i "$W/hvp.kt" --hier "$W/hv.tsv" --smooth --set repeat && ok "--hier: decompress --db same as the plain track" || no "--hier: decompress --db DIFFERS"
$KT --db "$D/db.features.tsv" --hier "$W/hv.tsv" --set repeat decompress -i "$W/hv.kt" 2>/dev/null | cmp -s - "$W/hv/repeat.bed" && ok "--hier: --set repeat gives the input" || no "--hier: --set repeat DIFFERS"
mkdir -p "$W/so"; $KT --db "$D/db.features.tsv" --set repeat,gene --out-dir "$W/so/new" decompress -i "$W/f.kt" 2>/dev/null; [ "$(ls "$W/so/new" | tr '\n' ' ')" = "gene.bed repeat.bed " ] && cmp -s "$W/so/new/gene.bed" "$D/db.gene.bed" && ok "decompress --set --out: only the picked sets, directory created" || no "decompress --set --out: WRONG"
# index needs no --aln or --db companions
cp "$W/d.kt" "$W/di.kt"; cp "$W/f.kt" "$W/fi.kt"
$KT index -i "$W/di.kt" >/dev/null 2>&1 && $KT index -i "$W/fi.kt" >/dev/null 2>&1 && seqcmp $KT --aln "$D/p.paf" --ref "$D/ref.bed" query -i "$W/di.kt" -r q1 -- $KT --aln "$D/p.paf" --ref "$D/ref.bed" query -i "$W/d.kt" -r q1 && ok "index on --aln and --db tracks: built, same answers" || no "index on --aln and --db tracks: FAILED"
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
for t in 2 4 8; do stress "decompress, $t threads: identical to serial" || continue
  $KT --threads $t decompress -i "$W/sm.kt" --hier "$W/sm.tsv" --smooth --set acrocentric > "$W/sm.$t" 2>/dev/null || { no "decompress, $t threads: FAILED"; continue; }
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
cp "$D/s.bed" "$W/keep.bed"; $KT compress "$W/keep.bed" < /dev/null > /dev/null 2>&1 && no "compress FILE: ACCEPTED"
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
rc=0; timeout "$(tmo 10)" $KT decompress -i "$W/p.kt" --hier "$W/cyc.tsv" --smooth --set acrocentric >/dev/null 2>&1 || rc=$?
[ $rc = 1 ] && ok "cyclic hierarchy: refused" || no "cyclic hierarchy: rc=$rc"
for c in "decompress -i $W/p.kt" "query -i $W/p.kt -r q1" "decompress -i $W/sm.kt --hier $W/sm.tsv --smooth --set acrocentric" "decompress --raw -i $W/p.kt"; do
  $KT $c > /dev/full 2>/dev/null && no "$c to /dev/full: EXIT 0" || ok "${c%% *} to /dev/full: exit 1"
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
# a ten-byte varint past 64 bits in the index table, re-signed: once read as 0, now refused; one that fits still reads
$KT --block 4 compress < "$D/s.bed" > "$W/v10.kt" 2>/dev/null
for t in 0x02 0x40 0x7f 0x81; do
  python3 "$D/kt.py" varint10 "$W/v10.kt" "$W/v10.$t.kt" $t "$W"
  must_fail "a ten-byte varint ending $t in the index table" $SAN decompress -i "$W/v10.$t.kt"
  grep -qx 'karyozip: corrupt varint' "$W/mf.err" && ok "a ten-byte varint ending $t: says corrupt varint" || no "a ten-byte varint ending $t: WRONG words: $(head -c 200 "$W/mf.err")"
done
python3 "$D/kt.py" varint10 "$W/v10.kt" "$W/v10.0.kt" 0x00 "$W"
$SAN decompress -i "$W/v10.0.kt" 2>/dev/null | cmp -s - "$D/s.bed" && ok "a ten-byte varint of 0 in the index table: still read" || no "a ten-byte varint of 0 in the index table: WRONG"

# a packed alignment needs the same reference
{ cat "$D/ref.bed"; printf 'chr2\t0\t50\tB\n'; } > "$W/ref2.bed"; { printf 'chr2\t0\t50\tB\n'; cat "$D/ref.bed"; } > "$W/ref.rev"
$KT aln-pack -a "$D/p.paf" --ref "$W/ref2.bed" -o "$W/p.kta" >/dev/null 2>&1
$KT compress --aln "$W/p.kta" --ref "$W/ref2.bed" < "$D/s.bed" > "$W/pk.kt" 2>/dev/null && $KT decompress --aln "$W/p.kta" --ref "$W/ref2.bed" -i "$W/pk.kt" 2>/dev/null | cmp -s - "$D/s.bed" \
  && ok "packed alignment: round-trip" || no "packed alignment: WRONG"
must_fail "packed alignment with a reordered reference" $KT decompress --aln "$W/p.kta" --ref "$W/ref.rev" -i "$W/pk.kt"
{ printf 'KSTAL2'; tail -c +7 "$W/p.kta"; } > "$W/old.kta"
must_fail "a packed alignment with another magic" $KT decompress --aln "$W/old.kta" --ref "$W/ref2.bed" -i "$W/pk.kt"
# a forged stream length is a clean refusal, not an oom; an empty stream's bytes are skipped, as the Rust
head -2 "$D/p.paf" > "$W/p2.paf"   # one gapless block a query: the deltas stream is empty
$KT aln-pack -a "$W/p2.paf" --ref "$W/ref2.bed" -o "$W/p2.kta" >/dev/null 2>&1
$KT compress --aln "$W/p2.kta" --ref "$W/ref2.bed" < "$D/s.bed" > "$W/pk2.kt" 2>/dev/null
python3 - "$D" "$W/p.kta" "$W/big.kta" "$W/p2.kta" "$W/skip.kta" <<'PY'
import sys
sys.path.insert(0, sys.argv[1]); import kt
def secs(p):
    d = open(p, "rb").read(); body = d[14:]; i = 0; sec = []
    for _ in range(10):
        cd, i = kt.vget(body, i); rl, i = kt.vget(body, i); cl, i = kt.vget(body, i); sec.append((cd, rl, body[i:i + cl])); i += cl
    return d[:14], sec
pack = lambda hd, s: hd + b"".join(kt.vput(c) + kt.vput(r) + kt.vput(len(x)) + x for c, r, x in s)
hd, sec = secs(sys.argv[2]); open(sys.argv[3], "wb").write(pack(hd, [(1, 1 << 40, b"")] + sec[1:]))
hd, sec = secs(sys.argv[4]); k = [j for j, x in enumerate(sec) if x[0] == 0]
if k: sec[k[0]] = (0, 0, b"junk")
open(sys.argv[5], "wb").write(pack(hd, sec)); print(k[0] if k else -1)
PY
must_fail "a packed alignment claiming a 1 TB stream" $KT decompress --aln "$W/big.kta" --ref "$W/ref2.bed" -i "$W/pk.kt"
grep -q 'impossible allocation size' "$W/mf.err" && ok "a packed alignment claiming a 1 TB stream: the allocation cap" || no "a packed alignment claiming a 1 TB stream: not the allocation cap"
$KT decompress --aln "$W/skip.kta" --ref "$W/ref2.bed" -i "$W/pk2.kt" 2>/dev/null | cmp -s - "$D/s.bed" && ! cmp -s "$W/skip.kta" "$W/p2.kta" \
  && ok "packed alignment, bytes behind an empty stream: skipped" || no "packed alignment, bytes behind an empty stream: WRONG"

# a header flag this reader does not know: a newer writer, refused before any size is read
for bit in 53 54; do
  python3 -c "import sys; d = bytearray(open(sys.argv[1], 'rb').read()); d[16 + $bit // 8] |= 1 << ($bit % 8); open(sys.argv[2], 'wb').write(d)" "$W/p.kt" "$W/fl$bit.kt"
  must_fail "header flag bit $bit" $KT decompress -i "$W/fl$bit.kt"
  grep -q 'written by a newer karyozip: unknown header flags' "$W/mf.err" && ok "header flag bit $bit: named as a newer writer" || no "header flag bit $bit: WRONG words"
done
# more than 64 -i: --db takes one BED per feature set
python3 - "$W" <<'PY'
import sys
W = sys.argv[1]; n = 65
open(W + "/f65.tsv", "w").write("featureID\t" + "\t".join("s%d" % k for k in range(n)) + "\n"
                                 + "".join("%d\t" % i + "\t".join("L%d" % ((i + k) % 3) for k in range(n)) + "\n" for i in range(1, 4)))
for k in range(n): open(W + "/f65.%d.bed" % k, "w").write("".join("c\t%d\t%d\tL%d\n" % (10 * j, 10 * j + 10, (j % 3 + k) % 3) for j in range(6)))
PY
in65=""; for k in $(seq 0 64); do in65="$in65 -i $W/f65.$k.bed"; done
mkdir -p "$W/o65"
$KT compress --db "$W/f65.tsv" $in65 > "$W/f65.kt" 2>/dev/null && $KT decompress --db "$W/f65.tsv" --out-dir "$W/o65" -i "$W/f65.kt" 2>/dev/null
n=0; for k in $(seq 0 64); do cmp -s "$W/o65/s$k.bed" "$W/f65.$k.bed" && n=$((n+1)); done
[ "$n" = 65 ] && ok "--db with 65 feature sets: 65/65 identical" || no "--db with 65 feature sets: only $n/65 identical"
# the default thread count is the CPUs this process may use, as the Rust's: the bare binary, pinned to one CPU
if command -v taskset >/dev/null; then
  taskset -c "$(python3 -c 'import os; print(min(os.sched_getaffinity(0)))')" $KZ --profile compress < "$D/s.bed" 2>&1 >/dev/null | grep -q '(1 allowed)' \
    && ok "default threads under a one-CPU affinity: 1" || no "default threads under a one-CPU affinity: WRONG"
fi

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
if stress "long temp dir, 80000 rows: round-trip, bystander intact"; then
awk 'BEGIN{for(i=0;i<80000;i++)printf "s%d\t0\t10\tA\n",i}' > "$W/long.bed"
$KT compress --tmp-dir "$L" --names drop < "$W/long.bed" > "$W/long.kt" 2>/dev/null \
  && $KT decompress --names <(cut -f1 "$W/long.bed") -i "$W/long.kt" 2>/dev/null | cmp -s - "$W/long.bed" && [ "$(cat "${L:0:511}" 2>/dev/null)" = keep ] \
  && ok "long temp dir: round-trip, bystander intact" || no "long temp dir: WRONG"
fi
$KT compress < "$D/s.bed" >&- 2>/dev/null && no "compress with stdout closed: EXIT 0" || ok "compress with stdout closed: refused"
$KT compress < "$W" > "$W/dir.kt" 2>/dev/null && no "compress from a directory: EXIT 0" || ok "compress from a directory: refused"
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
$KT decompress --raw -i "$W/f.kt" 2>/dev/null > "$W/merged.bed"; $KT compress --db "$D/db.features.tsv" -i "$W/merged.bed" > "$W/filedb.kt" 2>/dev/null
cat "$W/merged.bed" | $KT compress --db "$D/db.features.tsv" -i - > "$W/pipedb.kt" 2>/dev/null   # one karyozip at a time
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
# 301 = the 300 L's left once --v1 strips "_specific", plus cut's newline
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
$KT decompress --db "$D/db.features.tsv" --out-dir "$W/keepdir" -i "$W/p.kt" >/dev/null 2>&1 && no "decompress --db of a plain track: EXIT 0"
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
if stress "file-size limit, 100000 rows: exit 1, not a signal"; then
awk 'BEGIN{srand(3);for(i=0;i<20000;i++){p=0;for(j=0;j<5;j++){c=p+int(rand()*300)+1;printf "read%d\t%d\t%d\tL%d\n",i,p,c,int(rand()*6);p=c}}}' > "$W/pbig.bed"
rc=0; (ulimit -f 1; $KT compress < "$W/pbig.bed" > "$W/lim.kt" 2> "$W/lim.err") || rc=$?
[ $rc = 1 ] && ok "file-size limit: exit 1, not a signal" || no "file-size limit: rc=$rc"
fi
# ls.bed (a fixture) at each effort: a higher effort tries a superset of codecs, so it is never larger
sz=""; for e in high medium low; do sz="$sz $($KT -e $e compress < "$D/ls.bed" 2>/dev/null | wc -c)"; done
[ "$(echo $sz | tr ' ' '\n' | sort -n | tr '\n' ' ')" = "$(echo $sz) " ] && ok "effort high <= medium <= low:$sz B" || no "effort sizes out of order:$sz B"
# -e high: zpaq on both frames
if stress "-e high on six contigs of 20000 to 40000 rows: zpaq on both frames, round-trip"; then
awk 'BEGIN { srand(5); for (c = 0; c < 6; c++) { g = c % 4; nr = 20000 + int(rand() * 20000); p = 0; for (r = 0; r < nr; r++) { l = 1 + int(rand() * (g ? 9 * g : 37)); u = rand(); printf "ctg%d\t%d\t%d\t%s\n", c, p, p + l, u < 0.5 + g / 10 ? "chr" (1 + g) : "novel" int(rand() * 3); p += l } } }' > "$W/zt.bed"
$KT -e high compress -i "$W/zt.bed" -o "$W/zt.kt" 2> "$W/zt.err"
grep -q "lengths zpaq-64MB, labels zpaq-64MB" "$W/zt.err" && ok "high: zpaq on both frames of six contigs" || no "high: zpaq not picked: $(grep layout "$W/zt.err")"
$KT decompress -i "$W/zt.kt" 2>/dev/null | cmp -s - "$W/zt.bed" && ok "zpaq frames: round-trip" || no "zpaq frames: WRONG"
fi
# -e high: interleaved zpaq model, short reads
if stress "-e high on 60000 short reads: the interleaved zpaq model, round-trip"; then
awk 'BEGIN { srand(7); for (r = 0; r < 60000; r++) { k = rand() < 0.4 ? 1 : 2 + int(rand() * rand() * 12); p = 0; g = int(rand() * 4); for (j = 0; j < k; j++) { l = j == k - 1 ? 120 - p : 1 + int(rand() * (120 - p - (k - j - 1)) / (k - j)); printf "r%d\t%d\t%d\t%s\n", r, p, p + l, "chr" (1 + (g + j) % 5); p += l } } }' > "$W/zi.bed"
$KT -e high compress -i "$W/zi.bed" -o "$W/zi.kt" 2> "$W/zi.err"
# the layout line on stderr names the model the search picked; this input is built so the interleaved one wins
grep -q "lengths zpaq-64MB+ilv" "$W/zi.err" && ok "high: the interleaved zpaq model on short reads" || no "high: interleaved model not picked: $(grep layout "$W/zi.err")"
$KT decompress -i "$W/zi.kt" 2>/dev/null | cmp -s - "$W/zi.bed" && ok "interleaved zpaq model: round-trip" || no "interleaved zpaq model: WRONG"
fi
must_fail "--tmp-dir that does not exist" $KT --tmp-dir /nonexistent/dir compress < "$D/s.bed"
$KT --tmp-dir "$W" compress < "$D/s.bed" 2>/dev/null | cmp -s - "$W/p.kt" && ok "--tmp-dir honoured: same bytes" || no "--tmp-dir honoured: WRONG"

# round 25 findings, one case each
# environment variables are ignored; OMP_NUM_THREADS is a stress case: the bare binary at its default count, and -t 2
ok=1; stress "OMP_NUM_THREADS ignored, at the default thread count and at --threads 2" && dflt=1 || dflt=0
for e in 1 3 16 100000; do
  [ $dflt = 0 ] || OMP_NUM_THREADS=$e $KZ compress < "$D/db.rall.bed" 2>/dev/null | cmp -s - "$W/omp.ref" || ok=0
  [ $dflt = 0 ] || OMP_NUM_THREADS=$e $KT --threads 2 compress < "$D/db.rall.bed" 2>/dev/null | cmp -s - "$W/omp.ref" || ok=0
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
$KT decompress --db "$D/db.features.tsv" --out-dir "$W/good" -i "$W/fbad.kt" >/dev/null 2>&1 && no "decompress of a damaged db track: EXIT 0"
[ "$(cat "$W/good"/*.bed | md5sum)" = "$m1" ] && [ -z "$(ls "$W/good" | grep tmp)" ] && ok "failed decompress --db: previous output intact, no temporaries" || no "failed decompress --db: OUTPUT DAMAGED"
mkdir -p "$W/oo"; echo old > "$W/oo/k.kt"; $KT compress -i "$W/nonexist" -o "$W/oo/k.kt" >/dev/null 2>&1 || :
$KT decompress -i "$W/fbad.kt" -o "$W/oo/k.kt" >/dev/null 2>&1 || :; $KT compress -i "$W/nonexist" -o "$W/oo/n.kt" >/dev/null 2>&1 || :
[ "$(cat "$W/oo/k.kt")" = old ] && [ "$(ls "$W/oo")" = k.kt ] && ok "failed -o: previous file intact, no temporaries" || no "failed -o: OUTPUT DAMAGED"
$KT decompress -i "$W/nx.kt" -o "$W/oo/k.kt" >/dev/null 2>&1 || :; [ "$(ls "$W/oo")" = k.kt ] && [ -s "$W/oo/k.kt" ] && ok "-o over an old file: replaced" || no "-o over an old file: WRONG"

# output and temp-file safety: a failed -o, --out-dir, index or aln-pack never harms what was there
mkdir -p "$W/os"; cd "$W/os"
printf 'r0\tx\t5\tA\n' > bad.bed; echo old > k.kt; cp bad.bed same.bed
must_fail "-o, malformed input" $KT compress -i bad.bed -o k.kt
must_fail "-o names the input, malformed" $KT compress -i same.bed -o same.bed
[ "$(cat k.kt)" = old ] && cmp -s same.bed bad.bed && [ "$(ls | tr '\n' ' ')" = "bad.bed k.kt same.bed " ] \
  && ok "failed -o: sentinel and input intact, no temporaries" || no "failed -o: DAMAGED: $(ls | tr '\n' ' ')"
if [ "$(id -u)" != 0 ]; then   # root ignores the mode
  mkdir ro; echo old > ro/k.kt; cp "$D/s.bed" ro/s.bed; chmod 555 ro
  must_fail "-o in a directory that forbids new files" $KT compress -i "$D/s.bed" -o ro/k.kt
  must_fail "-o names the input, no staging" $KT compress -i ro/s.bed -o ro/s.bed
  must_fail "aln-pack in a directory that forbids new files" $KT aln-pack -a "$D/p.paf" --ref "$D/ref.bed" -o ro/k.kt
  [ "$(cat ro/k.kt)" = old ] && cmp -s ro/s.bed "$D/s.bed" && ok "-o unstageable: never truncated" || no "-o unstageable: TRUNCATED"
  chmod 755 ro
fi
# the old names were TARGET.tmp.PID, TRACK.kzi.PID and SET.bed.tmp.PID: a link planted there is not followed
echo victim > victim
( ln -s "$W/os/victim" "k.kt.tmp.$BASHPID"; exec $KT compress -i "$D/s.bed" -o k.kt ) 2>/dev/null
cp "$W/p.kt" ix.kt; ( ln -s "$W/os/victim" "ix.kt.kzi.$BASHPID"; exec $KT index -i ix.kt ) 2>/dev/null
mkdir -p sd; ( ln -s "$W/os/victim" "sd/gene.bed.tmp.$BASHPID"; exec $KT decompress --db "$D/db.features.tsv" --out-dir sd -i "$W/f.kt" ) 2>/dev/null
[ "$(cat victim)" = victim ] && cmp -s k.kt "$W/p.kt" && [ -s ix.kt.kzi ] && cmp -s sd/gene.bed "$D/db.gene.bed" \
  && [ "$(ls | grep -c 'k\.kt\.tmp\|kzi\.')" = 2 ] && [ "$(ls sd | grep -c tmp)" = 1 ] && ok "a link at an old temporary name: target untouched" || no "a link at an old temporary name: FOLLOWED"
# -o through a symlink: the target replaced, the link kept; a dangling link refused
echo old > tgt.kt; ln -s tgt.kt lnk.kt; $KT compress -i "$D/s.bed" -o lnk.kt 2>/dev/null
[ -L lnk.kt ] && cmp -s tgt.kt "$W/p.kt" && [ -z "$(ls | grep 'tgt.kt.tmp\|lnk.kt.tmp')" ] && ok "-o to a symlink: target replaced, link kept" || no "-o to a symlink: WRONG"
ln -s nowhere/x.kt dang.kt; must_fail "-o to a dangling symlink" $KT compress -i "$D/s.bed" -o dang.kt
[ -L dang.kt ] && [ ! -e nowhere ] && ok "dangling -o link: left alone" || no "dangling -o link: WRONG"
# aln-pack: staged as -o is
$KT aln-pack -a "$D/p.paf" --ref "$D/ref.bed" -o pk 2>/dev/null; cp pk pk.want
must_fail "aln-pack of a missing PAF" $KT aln-pack -a nosuch.paf --ref "$D/ref.bed" -o pk
cmp -s pk pk.want && [ -z "$(ls | grep 'pk.tmp')" ] && ok "failed aln-pack: old file intact" || no "failed aln-pack: OLD FILE DAMAGED"
echo old > pk2; $KT aln-pack -a "$D/p.paf" --ref "$D/ref.bed" -o pk2 2>/dev/null; cmp -s pk2 pk.want && ok "aln-pack over an old file: replaced" || no "aln-pack over an old file: WRONG"
echo old > od.bed
must_fail "-o with --out-dir" $KT decompress --db "$D/db.features.tsv" --out-dir od -o od.bed -i "$W/f.kt"
must_fail "--out-dir ''" $KT decompress --db "$D/db.features.tsv" --out-dir "" -i "$W/f.kt"
must_fail "--smooth --out-dir ''" $KT --db "$D/db.features.tsv" decompress -i "$W/f.kt" --hier "$W/full.tsv" --smooth --out-dir ""
[ "$(cat od.bed)" = old ] && [ ! -e od ] && ok "-o with --out-dir: refused before writing" || no "-o with --out-dir: WROTE"
cd - >/dev/null

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
stress "the documented limits: accepted" && { $KT -t 1024 compress < "$D/s.bed" > /dev/null 2>&1 && ok "the documented limits: accepted" || no "the documented limits: REFUSED"; }

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
$KT compress < "$D/s.bed" 2>&1 >/dev/null | grep -qi "htscodecs" && no "compress stderr: HTSCODECS COMPLAINED" || ok "compress stderr: clean"

# library API: same bytes as the CLI; the library has no -t, so its test runs on KZ_TEST_THREADS CPUs
if [ -x "$LT" ]; then
  LTB=$LT; LT="$(pin "$KZ_TEST_THREADS") $LTB"
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
  # the first through the one-block layout search, the second many blocks: the CLI's bytes and its stderr, totals included
  $KT --block 4 compress < "$D/db.rall.bed" > "$W/sq.want" 2> "$W/sq.wanterr"
  $LT -seq "$D/s.bed" "$W/sq1.kt" "$D/db.rall.bed" "$W/sq2.kt" 4 2> "$W/sq.err"
  cmp -s "$W/sq2.kt" "$W/sq.want" && grep -q '^layout: interleaved' "$W/sq.err" \
    && tail -n "$(wc -l < "$W/sq.wanterr")" "$W/sq.err" | cmp -s - "$W/sq.wanterr" \
    && ok "library API: a many-block writer after a one-block one: the CLI's bytes and messages" \
    || no "library API: a many-block writer after a one-block one keeps the first's state"
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
  no "library API: FAILED, $LT not built"
fi

# label cache collisions, prefix names, chunk boundaries
printf 'r1\t0\t10\tabcXYZ\nr1\t10\t20\tabqXYZ\nr1\t20\t30\tabcXYZ\n' > "$W/coll.bed"
$KT compress < "$W/coll.bed" > "$W/coll.kt" 2>/dev/null
$KT decompress -i "$W/coll.kt" 2>/dev/null | cmp -s - "$W/coll.bed" && ok "label cache: colliding keys round-trip" || no "label cache: COLLIDING KEYS WRONG"
if stress "label cache, 20000 labels on 60000 rows: round-trip"; then
python3 -c "
p=0
f=open('$W/many.bed','w')
for i in range(60000):
    f.write('r1\t%d\t%d\ta%04dz\n' % (p, p+10, i % 20000)); p += 10
"
$KT compress < "$W/many.bed" > "$W/many.kt" 2>/dev/null
$KT decompress -i "$W/many.kt" 2>/dev/null | cmp -s - "$W/many.bed" && ok "label cache: 20000 labels round-trip" || no "label cache: 20000 labels WRONG"
fi
printf 'r1\t0\t10\tA\nr11\t0\t10\tB\nr111\t0\t10\tC\n' > "$W/pre.bed"
$KT compress < "$W/pre.bed" > "$W/pre.kt" 2>/dev/null
$KT decompress -i "$W/pre.kt" 2>/dev/null | cmp -s - "$W/pre.bed" && ok "prefix sequence names: round-trip" || no "prefix sequence names: WRONG"
if stress "200000 one-row sequences at 8 threads: round-trip"; then
python3 -c "
f=open('$W/chunk.bed','w')
for s in range(200000): f.write('seq%06d\t0\t10\tA\n' % s)
"
$KT -t 8 compress < "$W/chunk.bed" > "$W/chunk.kt" 2>/dev/null
$KT decompress -i "$W/chunk.kt" 2>/dev/null | cmp -s - "$W/chunk.bed" && ok "200000 one-row sequences at 8 threads: round-trip" || no "200000 one-row sequences: WRONG"
fi
# a full temp filesystem must be fatal
mkdir -p "$W/tinytmp"
if stress "a full staging filesystem, 400000 rows: refused, no temp files left"; then
awk 'BEGIN{p=0; for(k=0;k<400000;k++){printf "r1\t%d\t%d\tchr%d\n", p, p+20, k%22; p+=20}}' > "$W/stage.bed"
# 256 KB cap, below the staging payload
( ulimit -f 512; trap "" XFSZ; $KT --tmp-dir "$W/tinytmp" compress < "$W/stage.bed" >/dev/null 2>&1 ) && no "a full staging filesystem: ACCEPTED" \
  || ok "a full staging filesystem: refused"
[ -z "$(ls -A "$W/tinytmp" 2>/dev/null)" ] && ok "a full staging filesystem: no temp files left" || no "a full staging filesystem: TEMP FILES LEFT"
fi
# the staging fits, the compressed body does not: one-row blocks make it three times the staged bytes; -t 2, the ordered loop
if stress "a full body temp file at 2 threads"; then
awk 'BEGIN{p=0; for(k=0;k<20000;k++){printf "r1\t%d\t%d\tchr%d\n", p, p+20, k%22; p+=20}}' > "$W/body.bed"
rc=0; ( ulimit -f 64; trap "" XFSZ; timeout 60 $KT -t 2 --tmp-dir "$W/tinytmp" --block 1 compress < "$W/body.bed" >/dev/null 2>"$W/body.err" ) || rc=$?
[ $rc = 1 ] && grep -qx 'karyozip: tmp write' "$W/body.err" && ok "a full body temp file: refused at once (tmp write)" || no "a full body temp file: WRONG, exit $rc: $(head -c 200 "$W/body.err")"
[ -z "$(ls -A "$W/tinytmp" 2>/dev/null)" ] && ok "a full body temp file: no temp files left" || no "a full body temp file: TEMP FILES LEFT"
fi

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
for t in 1 4; do [ $t = 1 ] || stress "--read-blocks 7 --smooth -t $t: same as one block" || continue; seqcmp $KT decompress -t $t -i "$W/vb.kt" --smooth --hier "$W/vbh.tsv" -- $KT decompress -i "$W/vb1.kt" --smooth --hier "$W/vbh.tsv" && ok "--read-blocks 7 --smooth -t $t: same as one block" || no "--read-blocks 7 --smooth -t $t: DIFFERS"; done

# decode ahead: same rows, same stop
if stress "decode ahead at -t 4: same rows as -t 1, same stop, sanitizer clean"; then
awk 'BEGIN { x = 7; for (r = 0; r < 120000; r++) { n = 1 + r % 4; p = 0; for (i = 0; i < n; i++) { x = (x * 1103515245 + 12345) % 2147483648; l = i + 1 == n ? 151 - p : 1 + x % 37; printf "q%d\t%d\t%d\t%s\n", r, p, p + l, x % 7 ? "chr" (1 + x % 3) : "novel"; p += l } } }' > "$W/pf.bed"
for bl in "--read-blocks 70000" "--block 70000"; do
  $KT compress $bl -i "$W/pf.bed" -o "$W/pf.kt" 2>/dev/null
  $KT decompress -t 4 -i "$W/pf.kt" 2>/dev/null | cmp -s - "$W/pf.bed" && ok "decode ahead [$bl]: byte-identical" || no "decode ahead [$bl]: WRONG"
  cp "$W/pf.kt" "$W/pfb.kt"; sz=$(stat -c%s "$W/pfb.kt"); printf '\xff' | dd of="$W/pfb.kt" bs=1 seek=$((sz - 300)) conv=notrunc status=none
  a=0; $KT decompress -t 1 -i "$W/pfb.kt" > "$W/pf1.out" 2> "$W/pf1.err" || a=$?
  b=0; $KT decompress -t 4 -i "$W/pfb.kt" > "$W/pf4.out" 2> "$W/pf4.err" || b=$?
  cmp -s "$W/pf1.out" "$W/pf4.out" && cmp -s "$W/pf1.err" "$W/pf4.err" && [ $a = $b ] && [ $a = 1 ] && ok "decode ahead [$bl], damaged: stops at the same byte, $(wc -c < "$W/pf1.out") B" || no "decode ahead [$bl], damaged: DIFFERS from one thread"
done
if [ "$KZSAN" != "$KZ" ]; then ASAN_OPTIONS=detect_leaks=1 $SAN decompress -t 4 -i "$W/pf.kt" 2> "$W/pfs.err" | cmp -s - "$W/pf.bed" && ! grep -qE "Sanitizer|runtime error" "$W/pfs.err" && ok "decode ahead: sanitizer clean" || no "decode ahead: SANITIZER REPORT"; fi
fi

# forged BGZF, under the sanitizer binary: ISIZE bounds, BC after another subfield, cuts inside the extra field
X=$W/bx; mkdir -p "$X"; python3 "$D/kt.py" bgzfx "$X"
bx_ok() {  # DESC GZ BED: file and pipe give the plain BED's bytes, sanitizer clean
  $KT compress -i "$3" > "$W/bx.want" 2>/dev/null
  ASAN_OPTIONS=detect_leaks=1 $SAN compress -i "$2" > "$W/bx.f" 2> "$W/bx.err"; ASAN_OPTIONS=detect_leaks=1 $SAN compress < "$2" > "$W/bx.p" 2>> "$W/bx.err"
  cmp -s "$W/bx.f" "$W/bx.want" && cmp -s "$W/bx.p" "$W/bx.want" && ! grep -qE "Sanitizer|runtime error" "$W/bx.err" && ok "$1: same bytes" || no "$1: WRONG"
}
bx_no() {  # DESC GZ: refused from file and pipe, sanitizer clean
  local rc=0 rp=0
  local A=detect_leaks=0   # a refusal exits through die, leaks are not the point here
  ASAN_OPTIONS=$A $SAN compress -i "$2" > /dev/null 2> "$W/bx.err" || rc=$?; ASAN_OPTIONS=$A $SAN compress < "$2" > /dev/null 2>> "$W/bx.err" || rp=$?
  [ $rc = 1 ] && [ $rp = 1 ] && [ "$(grep -c '^karyozip: ' "$W/bx.err")" = 2 ] && ! grep -qE "Sanitizer|runtime error" "$W/bx.err" \
    && ok "$1: refused ($(grep -m1 '^karyozip: ' "$W/bx.err" | cut -c11-90))" || no "$1: WRONG, exit $rc/$rp: $(head -c 300 "$W/bx.err")"
}
bx_ok "BGZF ISIZE 65535" "$X/isz65535.gz" "$X/isz65535.bed"
bx_ok "BGZF ISIZE 65536" "$X/isz65536.gz" "$X/isz65536.bed"
bx_no "BGZF ISIZE 65537" "$X/isz65537.gz"
bx_no "BGZF BSIZE ending inside its header" "$X/bsize.gz"
bx_ok "BGZF BC after another subfield" "$X/bcafter.gz" "$X/bc.bed"
bx_no "BGZF BC payload past the end" "$X/bcend1.gz"; bx_no "BGZF BC payload past the end, second block" "$X/bcend2.gz"
cb=0; for f in "$X"/cut[12].*.gz; do bx_no "cut" "$f" > "$W/bx.one" || true; grep -q WRONG "$W/bx.one" && { cat "$W/bx.one"; cb=1; }; done
[ $cb = 0 ] && ok "BGZF cut at every byte of the extra field, one block and two: refused, sanitizer clean" || no "BGZF cut inside the extra field: WRONG"
if stress "BGZF blocks claiming 4 GiB each, 300000 blocks claiming 64 KiB each, under a 1.5 GB address limit at -t 2"; then
bx_no "BGZF blocks claiming 4 GiB each" "$X/huge.gz"; bx_no "BGZF 300000 blocks claiming 64 KiB each" "$X/many.gz"
for f in huge many; do   # refused before allocating what the trailers claim
  rc=0; ( ulimit -v 1500000; exec timeout 30 $KT compress -t 2 -i "$X/$f.gz" ) > /dev/null 2> "$W/bx.err" || rc=$?
  [ $rc = 1 ] && grep -q '^karyozip: ' "$W/bx.err" && ! grep -q '^karyozip: oom' "$W/bx.err" && ok "BGZF $f.gz under a 1.5 GB address limit: refused, not out of memory" || no "BGZF $f.gz under a 1.5 GB address limit: WRONG, exit $rc: $(head -c 200 "$W/bx.err")"
done
fi

# tests/kt.py's flag mask is kst.rs's FLAGS: forged fixtures re-sign with it
KST=$src/../rust/src/kst.rs
if [ -f "$KST" ]; then
  python3 - "$KST" "$D" <<'PY' && ok "kt.py FLAGS: same bits as rust/src/kst.rs" || no "kt.py FLAGS: DIFFERS from rust/src/kst.rs"
import re, sys; sys.path.insert(0, sys.argv[2]); import kt
s = open(sys.argv[1]).read()
bit = {n: 1 << int(b) for n, b in re.findall(r'pub const (\w+): u64 = 1 << (\d+);', s)}
expr = re.search(r'pub const FLAGS: u64 = ([^;]+);', s).group(1)
mask = 0
for n in re.findall(r'\w+', expr): mask |= bit[n]
sys.exit(0 if mask == kt.FLAGS else 1)
PY
else echo "kt.py FLAGS: skipped, no $KST"; fi

# existing behaviour with no test before: a name with ':', --max-gap 0, --tmp-dir on readers, -o to a FIFO, signals
printf 'HLA:A*01\t0\t10\tA\nHLA:A*01\t10\t20\tB\nx:5-6\t0\t10\tC\n' > "$W/colon.bed"
$KT compress < "$W/colon.bed" > "$W/colon.kt" 2>/dev/null
$KT query -i "$W/colon.kt" -r 'HLA:A*01' 2>/dev/null | cmp -s - <(head -2 "$W/colon.bed") \
  && $KT query -i "$W/colon.kt" -r 'HLA:A*01:12-13' 2>/dev/null | cmp -s - <(sed -n 2p "$W/colon.bed") \
  && ok "-r with a name holding ':': the name, and a final :START-END a range" || no "-r with a name holding ':': WRONG"
must_fail "-r NAME:5-6 is a range on NAME, not the name 'NAME:5-6'" $KT query -i "$W/colon.kt" -r x:5-6
$KT query -i "$W/colon.kt" -r 'x:5-6:0-3' 2>/dev/null | cmp -s - <(sed -n 3p "$W/colon.bed") && ok "-r 'x:5-6:0-3': a range on the name x:5-6" || no "-r 'x:5-6:0-3': WRONG"
# --max-gap 0: accepted, smooths less than the default, rows still tile the input
$KT --max-gap 0 decompress -i "$W/ch.kst" --hier "$W/ch.tsv" --smooth --set fs > "$W/mg0" 2>/dev/null && rc=0 || rc=$?
awk -F'\t' -v e="$(tail -1 "$W/ch.bed" | cut -f3)" 'NR == 1 && $2 != 0 {bad = 1} NR > 1 && $2 != p {bad = 1} {p = $3} END {exit bad || p != e}' "$W/mg0" \
  && [ $rc = 0 ] && [ "$(wc -l < "$W/mg0")" -gt "$(wc -l < "$W/ch.whole")" ] && ! cmp -s "$W/mg0" "$W/ch.whole" \
  && ok "--max-gap 0: accepted, $(wc -l < "$W/mg0") rows tiling the input, more than the default's $(wc -l < "$W/ch.whole")" || no "--max-gap 0: WRONG, rc=$rc"
# readers make no temporary file: --tmp-dir is accepted and ignored, even a missing one
mkdir -p "$W/rtd"; rok=1
for c in "decompress -i $W/p.kt" "query -i $W/p.kt -r q1" "decompress -i $W/ch.kst --hier $W/ch.tsv --smooth --set fs"; do
  $KT $c > "$W/rtd.want" 2>/dev/null
  $KT --tmp-dir "$W/rtd" $c 2>/dev/null | cmp -s - "$W/rtd.want" || rok=0
  $KT --tmp-dir /nonexistent $c 2>/dev/null | cmp -s - "$W/rtd.want" || rok=0
done
[ $rok = 1 ] && [ -z "$(ls -A "$W/rtd")" ] && ok "--tmp-dir on readers: same rows, nothing written, a missing directory ignored" || no "--tmp-dir on readers: WRONG"
# -o to a FIFO: written in place, no temporary beside it
mkdir -p "$W/ff"; mkfifo "$W/ff/o.bed" "$W/ff/o.kt"
timeout "$(tmo 30)" cat "$W/ff/o.bed" > "$W/ff.got" & cpid=$!
$KT decompress -i "$W/p.kt" -o "$W/ff/o.bed" 2>/dev/null && rc=0 || rc=$?; wait $cpid || :
timeout "$(tmo 30)" cat "$W/ff/o.kt" > "$W/ff.kt" & cpid=$!
$KT compress -i "$D/s.bed" -o "$W/ff/o.kt" 2>/dev/null && rc2=0 || rc2=$?; wait $cpid || :
[ $rc = 0 ] && [ $rc2 = 0 ] && cmp -s "$W/ff.got" "$D/s.bed" && cmp -s "$W/ff.kt" "$W/p.kt" && [ "$(ls "$W/ff" | tr '\n' ' ')" = "o.bed o.kt " ] \
  && ok "-o to a FIFO: decompress and compress write through it, no temporaries" || no "-o to a FIFO: WRONG, rc $rc/$rc2"
# SIGTERM while blocked: the name file is a FIFO nobody writes, so the reader waits with its temporaries open
# sigwait DIR GLOB PID: up to 30 s (scaled) for GLOB to show in DIR while PID lives
sigwait() { local n=0 lim; lim=$(awk -v t="$(tmo 30)" 'BEGIN { print int(t * 20) }')
  while [ $n -lt $lim ] && kill -0 "$3" 2>/dev/null; do compgen -G "$1/$2" >/dev/null && return 0; sleep 0.05; n=$((n+1)); done; return 1; }
mkdir -p "$W/sg"; $KT compress --names drop < "$W/ch.bed" > "$W/sg.kt" 2>/dev/null; mkfifo "$W/sg.names"
$KT --names "$W/sg.names" decompress -i "$W/sg.kt" -o "$W/sg/out.bed" 2>/dev/null & pid=$!
if sigwait "$W/sg" 'out.bed.tmp.*' $pid; then kill -TERM $pid 2>/dev/null || :; wait $pid && rc=0 || rc=$?
  [ $rc = 1 ] && [ -z "$(ls -A "$W/sg")" ] && ok "SIGTERM during decompress -o: exit 1, no temporary, no output" || no "SIGTERM during decompress -o: rc=$rc, left: $(ls "$W/sg")"
else kill -KILL $pid 2>/dev/null || :; wait $pid 2>/dev/null || :; no "SIGTERM during decompress -o: no temporary appeared"; fi
$KT compress --db "$D/db.features.tsv" --names drop $in > "$W/sgd.kt" 2>/dev/null
$KT --db "$D/db.features.tsv" --hier "$W/full.tsv" --names "$W/sg.names" decompress -i "$W/sgd.kt" --smooth --out-dir "$W/sgd" 2>/dev/null & pid=$!
if sigwait "$W/sgd" '*.bed.tmp.*' $pid; then kill -TERM $pid 2>/dev/null || :; wait $pid && rc=0 || rc=$?
  [ $rc = 1 ] && [ -z "$(ls -A "$W/sgd")" ] && ok "SIGTERM during decompress --out-dir: exit 1, no temporaries, no BEDs" || no "SIGTERM during decompress --out-dir: rc=$rc, left: $(ls "$W/sgd")"
else kill -KILL $pid 2>/dev/null || :; wait $pid 2>/dev/null || :; no "SIGTERM during decompress --out-dir: no temporary appeared"; fi

# delayed write errors (tests/failclose.c, LD_PRELOAD): a staged output whose fsync or close fails is refused, every
# old output kept byte for byte, no temporaries left; -o, --out-dir, --smooth --out-dir, index and aln-pack
FC=$W/failclose.so
if "${CC:-cc}" -shared -fPIC -O2 -o "$FC" "$src/failclose.c" -ldl 2>"$W/fc.err"; then
  fcfail() {  # DESC MESSAGE MATCH CMD...: exit 1 and MESSAGE alone, the shim failing $op on MATCH's descriptors
    local what=$1 msg=$2 m=$3 rc=0; shift 3
    env KZ_FAIL_FD_OP=$op KZ_FAIL_FD_MATCH="$m" LD_PRELOAD="$FC" "$@" >/dev/null 2>"$W/fc.err" || rc=$?
    if [ $rc = 1 ] && [ "$(grep '^karyozip: ' "$W/fc.err")" = "karyozip: $msg" ]; then ok "$what: refused ($msg)"
    else no "$what: WRONG, exit $rc: $(head -c 200 "$W/fc.err")"; fi
  }
  mkdir -p "$W/fc"; cd "$W/fc"
  for op in sync close; do
    rm -rf ./*; mkdir od sm
    for s in gene repeat; do echo "OLD_$s" > od/$s.bed; echo "OLD_$s" > sm/$s.bed; done
    echo stale > od/gene.bed.old.AAAAAA; echo old > k.kt; echo old > pk; cp "$W/p.kt" ix.kt; echo oldidx > ix.kt.kzi   # not an index: ignored when read
    fcfail "failed $op, decompress --db --out-dir" "write failed" .bed.tmp. $KT decompress --db "$D/db.features.tsv" --set gene,repeat --out-dir od -i "$W/f.kt"
    fcfail "failed $op, --smooth --out-dir" "write failed" .bed.tmp. $KT --db "$D/db.features.tsv" --set gene,repeat decompress -i "$W/f.kt" --hier "$W/full.tsv" --smooth --out-dir sm
    fcfail "failed $op, index" "index write failed" .kzi. $KT index -i ix.kt
    fcfail "failed $op, -o" "write failed" k.kt.tmp. $KT compress -i "$D/s.bed" -o k.kt
    fcfail "failed $op, aln-pack" "write failed" pk.tmp. $KT aln-pack -a "$D/p.paf" --ref "$D/ref.bed" -o pk
    [ "$(cat od/gene.bed od/repeat.bed od/gene.bed.old.AAAAAA sm/gene.bed sm/repeat.bed k.kt pk ix.kt.kzi | tr '\n' ' ')" = "OLD_gene OLD_repeat stale OLD_gene OLD_repeat old old oldidx " ] \
      && cmp -s ix.kt "$W/p.kt" && [ "$(ls -A . od sm | tr '\n' ' ')" = ".: ix.kt ix.kt.kzi k.kt od pk sm  od: gene.bed gene.bed.old.AAAAAA repeat.bed  sm: gene.bed repeat.bed " ] \
      && ok "failed $op: every old output intact, no temporaries" || no "failed $op: DAMAGED: $(ls -A . od sm | tr '\n' ' ')"
  done
  # the shim loaded, nothing failed: replaced, the stale backup kept
  env KZ_FAIL_FD_MATCH=nomatch LD_PRELOAD="$FC" $KT decompress --db "$D/db.features.tsv" --set gene,repeat --out-dir od -i "$W/f.kt" 2>/dev/null
  cmp -s od/gene.bed "$D/db.gene.bed" && cmp -s od/repeat.bed "$D/db.repeat.bed" && [ "$(cat od/gene.bed.old.AAAAAA)" = stale ] \
    && [ "$(ls -A od | tr '\n' ' ')" = "gene.bed gene.bed.old.AAAAAA repeat.bed " ] && ok "the shim, nothing failing: replaced, a stale backup kept" || no "the shim, nothing failing: WRONG: $(ls -A od | tr '\n' ' ')"
  cd - >/dev/null
else no "tests/failclose.c: does not build with ${CC:-cc}: $(head -c 200 "$W/fc.err")"; fi

# --out-dir's commit, all or nothing: KZ_TEST_SWAP_PAUSE_MS pauses after each of its renames, so the test acts in
# the window. Two old BEDs and two stale backups (of an earlier run, one with the old PID naming) beside them
PAUSE=$(awk -v s="${TIMEOUT_SCALE:-1}" 'BEGIN { print int(600 * s) }')
swdir() { rm -rf "$1"; mkdir -p "$1"; echo OLD_gene > "$1/gene.bed"; echo OLD_repeat > "$1/repeat.bed"
  echo stale > "$1/gene.bed.old.AAAAAA"; echo stale > "$1/repeat.bed.old.12345"; }
# swgo DIR: exec, so $! is karyozip itself
swgo() { KZ_TEST_SWAP_PAUSE_MS=$PAUSE exec $KT decompress --db "$D/db.features.tsv" --set gene,repeat --out-dir "$1" -i "$W/f.kt" >/dev/null 2>"$W/sw.err"; }
# swgen DIR: old, new or mixed; swclean DIR: nothing but the BEDs and the untouched stale backups
swgen() { if [ "$(cat "$1/gene.bed" "$1/repeat.bed" 2>/dev/null | tr '\n' ' ')" = "OLD_gene OLD_repeat " ]; then echo old
  elif cmp -s "$1/gene.bed" "$D/db.gene.bed" && cmp -s "$1/repeat.bed" "$D/db.repeat.bed"; then echo new; else echo mixed; fi; }
swclean() { [ "$(ls -A "$1" | tr '\n' ' ')" = "gene.bed gene.bed.old.AAAAAA repeat.bed repeat.bed.old.12345 " ] \
  && [ "$(cat "$1/gene.bed.old.AAAAAA" "$1/repeat.bed.old.12345" | tr '\n' ' ')" = "stale stale " ]; }
# swwait DIR COND PID: up to 30 s (scaled) for the shell test COND, run in DIR, while PID lives
swwait() { local n=0 lim; lim=$(awk -v t="$(tmo 30)" 'BEGIN { print int(t * 50) }')
  while [ $n -lt $lim ] && kill -0 "$3" 2>/dev/null; do (cd "$1" && eval "$2") && return 0; sleep 0.02; n=$((n+1)); done; return 1; }
nold='[ $(ls | grep -c "\.bed\.old\.") -gt 2 ]'; ntmp='[ $(ls | grep -c "\.bed\.tmp\.") -lt 2 ]'
swdir "$W/sw"; ( swgo "$W/sw" ) && [ "$(swgen "$W/sw")" = new ] && swclean "$W/sw" && ok "--out-dir over old BEDs and stale backups: replaced, the stale ones kept" || no "--out-dir over stale backups: WRONG: $(ls -A "$W/sw" | tr '\n' ' ')"
for at in backup install; do
  swdir "$W/sw"; ( swgo "$W/sw" ) & pid=$!
  if swwait "$W/sw" "$nold" $pid && { [ $at = backup ] || swwait "$W/sw" "$ntmp" $pid; }; then
    kill -TERM $pid 2>/dev/null || :; wait $pid && rc=0 || rc=$?; g=$(swgen "$W/sw")
    { [ $rc = 0 ] || [ $rc = 1 ]; } && [ $g != mixed ] && swclean "$W/sw" && ok "SIGTERM in --out-dir's $at renames: exit $rc, all $g, no leftovers" \
      || no "SIGTERM in --out-dir's $at renames: rc=$rc, $g, left: $(ls -A "$W/sw" | tr '\n' ' ')"
  else kill -KILL $pid 2>/dev/null || :; wait $pid 2>/dev/null || :; no "SIGTERM in --out-dir's $at renames: the window never showed"; fi
done
# a rename fails: a directory takes gene.bed's place inside the window; with no old gene.bed the rollback succeeds
swdir "$W/sw"; rm "$W/sw/gene.bed"; ( swgo "$W/sw" ) & pid=$!
if swwait "$W/sw" '[ ! -e repeat.bed ]' $pid; then mkdir -p "$W/sw/gene.bed/x"; wait $pid && rc=0 || rc=$?
  [ $rc = 1 ] && [ "$(grep '^karyozip: ' "$W/sw.err")" = "karyozip: cannot rename a derived BED into place" ] && [ "$(cat "$W/sw/repeat.bed")" = OLD_repeat ] \
    && rmdir "$W/sw/gene.bed/x" "$W/sw/gene.bed" && echo OLD_gene > "$W/sw/gene.bed" && swclean "$W/sw" \
    && ok "a failed --out-dir rename: rolled back, all old, no leftovers" || no "a failed --out-dir rename: rc=$rc, $(grep -m1 '^karyozip: ' "$W/sw.err" | head -c 200), left: $(ls -A "$W/sw" | tr '\n' ' ')"
else kill -KILL $pid 2>/dev/null || :; wait $pid 2>/dev/null || :; no "a failed --out-dir rename: the window never showed"; fi
# the rollback fails too: gene.bed's old file cannot go back over the directory, and the refusal says so
swdir "$W/sw"; ( swgo "$W/sw" ) & pid=$!
if swwait "$W/sw" '[ ! -e gene.bed ]' $pid; then mkdir -p "$W/sw/gene.bed/x"; wait $pid && rc=0 || rc=$?
  kept=$(ls "$W/sw" | grep -v AAAAAA | grep '^gene\.bed\.old\.' || :)
  [ $rc = 1 ] && [ "$(grep '^karyozip: ' "$W/sw.err")" = "karyozip: cannot rename a derived BED into place, nor put every old one back: old BEDs remain as *.bed.old.* under $W/sw" ] \
    && [ -n "$kept" ] && [ "$(cat "$W/sw/$kept")" = OLD_gene ] && [ "$(cat "$W/sw/repeat.bed")" = OLD_repeat ] \
    && rmdir "$W/sw/gene.bed/x" "$W/sw/gene.bed" && mv "$W/sw/$kept" "$W/sw/gene.bed" && swclean "$W/sw" \
    && ok "a failed --out-dir rollback: refused, the old BED kept as $kept" || no "a failed --out-dir rollback: rc=$rc, $(grep -m1 '^karyozip: ' "$W/sw.err" | head -c 200), left: $(ls -A "$W/sw" | tr '\n' ' ')"
else kill -KILL $pid 2>/dev/null || :; wait $pid 2>/dev/null || :; no "a failed --out-dir rollback: the window never showed"; fi

# a handler that has decided to exit races the commit's start: KZ_TEST_SIG_PAUSE_MS holds the handler right after its
# decision, on a worker thread (signalled by its thread id, the thread Linux then prefers), while KZ_TEST_SWAP_PAUSE_MS
# holds the main thread just before the window, then lets it go on into the commit. --smooth --out-dir over 100
# contigs, so the C too keeps a team of workers. Signalled once both BEDs are staged and closed, before the window:
# the handler exits and the commit must not start, all old (all new only if the window opened first). Signalled
# inside the window: the handler defers, all new, exit 1
mkdir -p "$W/rc"; printf 'featureID\tgene\trepeat\n1\tCDS\tLINE\n2\tCDS\tSINE\n3\tintergenic\tLINE\n4\tintergenic\tSINE\n' > "$W/rc/f.tsv"
printf 'feature_set\tchild\tparent\ngene\tCDS\tcategorized\ngene\tintergenic\tcategorized\nrepeat\tLINE\tcategorized\nrepeat\tSINE\tcategorized\n' > "$W/rc/h.tsv"
awk -v d="$W/rc" 'BEGIN { for (c = 0; c < 100; c++) for (j = 0; j < 4; j++) {
  printf "c%03d\t%d\t%d\t%s\n", c, 10 * j, 10 * j + 10, (c + j) % 2 ? "CDS" : "intergenic" > (d "/gene.bed")
  printf "c%03d\t%d\t%d\t%s\n", c, 10 * j, 10 * j + 10, j % 2 ? "LINE" : "SINE" > (d "/repeat.bed") } }'
$KT compress --db "$W/rc/f.tsv" -i "$W/rc/gene.bed" -i "$W/rc/repeat.bed" > "$W/rc/r.kt" 2>/dev/null
rcgo() { exec $KZ -t 4 --db "$W/rc/f.tsv" decompress -i "$W/rc/r.kt" --hier "$W/rc/h.tsv" --smooth --out-dir "$1" >/dev/null 2>"$W/sw.err"; }
rcrace() { KZ_TEST_SWAP_PAUSE_MS=$PAUSE KZ_TEST_SIG_PAUSE_MS=$((2 * PAUSE)) rcgo "$1"; }
rcgen() { if [ "$(cat "$1/gene.bed" "$1/repeat.bed" 2>/dev/null | tr '\n' ' ')" = "OLD_gene OLD_repeat " ]; then echo old
  elif cmp -s "$1/gene.bed" "$W/rc/new/gene.bed" && cmp -s "$1/repeat.bed" "$W/rc/new/repeat.bed"; then echo new; else echo mixed; fi; }
staged='[ $(ls | grep -c "\.bed\.tmp\.") = 2 ] && ! ls -l /proc/$pid/fd 2>/dev/null | grep -q "\.bed\.tmp\."'
( rcgo "$W/rc/new" ) || :
for at in before inside; do
  swdir "$W/sw"; ( rcrace "$W/sw" ) & pid=$!; [ $at = before ] && cond=$staged || cond=$nold; tid=
  swwait "$W/sw" "$cond" $pid && tid=$(ls /proc/$pid/task 2>/dev/null | grep -vx "$pid" | head -n1) || :
  if [ -n "$tid" ]; then
    kill -TERM "$tid" 2>/dev/null || :; wait $pid && rc=0 || rc=$?; g=$(rcgen "$W/sw")
    [ $rc = 1 ] && { [ $g = new ] || { [ $at = before ] && [ $g = old ]; }; } && swclean "$W/sw" \
      && ok "SIGTERM on a worker, its handler held, $at the commit window: exit 1, all $g, no leftovers" \
      || no "SIGTERM on a worker, its handler held, $at the commit window: rc=$rc, $g, left: $(ls -A "$W/sw" | tr '\n' ' ')"
  else kill -KILL $pid 2>/dev/null || :; wait $pid 2>/dev/null || :; no "SIGTERM on a worker $at the commit window: no window, or no worker thread"; fi
done

# pin the stored bytes to codec versions: `karyozip --block 4 compress` of the fixture db.rall.bed
# (tests/fixtures.txt), default effort (medium), no other flag but $KT's -t (threads never change the bytes).
# It moves when a vendored codec, zstd or the format changes; then check diff.sh still passes and update the pin on purpose.
want=55b6e9eadfe851a4374e094e045ee8ecb4ffac9a2fdd83de126ff076724e700b
got=$($KT --block 4 compress < "$D/db.rall.bed" 2>/dev/null | sha256sum | cut -d' ' -f1)
if [ "$got" = "$want" ]; then ok "codec versions: stored bytes unchanged"
else
  no "codec versions: STORED BYTES MOVED"
  echo "  want $want"; echo "  got  $got"; $KT --version | sed 's/^/  /'
fi

# fuzz.sh reads TEST_LIGHT itself: under TEST_LIGHT=1 it runs its light subset (10 cases per generator, no forged-index)
if [ -n "${CASE:-}" ] && ! [[ fuzz =~ $CASE ]]; then echo "fuzz: skipped (CASE=$CASE)"
elif CASE= bash "$D/fuzz.sh" "$KZ" "$KZSAN"; then ok "fuzz: clean"; else no "fuzz: FOUND SOMETHING"; fi

summary run.sh
exit $fail
