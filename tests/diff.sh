#!/bin/bash
# C and Rust: same bytes, same words
#   diff.sh [C_BINARY] [RUST_BINARY]
# Every compared form counts once; failures are prefixed "FAIL: "; the last lines are
# "diff.sh: N passed, M failed, K skipped" and DIFF_PASS or DIFF_FAIL.
# KZ_TEST_THREADS, TEST_LIGHT: see lib.sh. "$C" and "$R" are wrappers that add
#   -t $KZ_TEST_THREADS; $CB and $RB are the bare binaries, for the cli section's exact argument lists.
# CASE=REGEX runs only the sections whose name matches: versions compress decompress cross region index
#   reader-flags smooth db aln refusals awkward gzip cli tail zpaq damaged forged edge bodytmp forged2 varint finalize
set -u
. "$(dirname "$0")/lib.sh"
need python3 zstd gzip awk realpath mktemp cmp diff dd stat cut sort uniq sed tr head tail
SCASE=${CASE:-}; CASE=   # sections, not lines, are filtered here
sec() { [ -z "$SCASE" ] || [[ $1 =~ $SCASE ]]; }
C=${1:-$(dirname "$0")/../c/build/karyozip}
R=${2:-$(dirname "$0")/../rust/target/release/karyozip}
D=$(cd "$(dirname "$0")" && pwd)
[ -x "$C" ] && [ -x "$R" ] || { echo "build both first: make all rust"; exit 1; }
C=$(realpath "$C"); R=$(realpath "$R")
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
CB=$C; RB=$R; mkdir -p "$W/bin"
for x in C R; do printf '#!/bin/sh\nexec "%s" -t %s "$@"\n' "${!x}" "$KZ_TEST_THREADS" > "$W/bin/$x"; chmod +x "$W/bin/$x"; done
C=$W/bin/C; R=$W/bin/R   # a later -t wins, in the C and the Rust
src=$D; D=$W/fx; mkdir -p "$D"; cp "$src"/fixtures.txt "$src"/*.py "$src"/*.sh "$D"/   # the fixtures, unpacked
awk '/^@@ /{f=d"/"$2; next} {print > f}' d="$D" "$D/fixtures.txt"

cv=$("$C" --version | awk '/zstd|libbz3|liblzma/{print $1"="$2}' | tr '\n' ' ')
rv=$("$R" --version | awk '/zstd|libbz3|liblzma/{print $1"="$2}' | tr '\n' ' ')
if [ "$cv" = "$rv" ]; then tpass "codec versions: $cv"
else tfail "codec versions DIFFER: c=[$cv] rust=[$rv]"; fi

# both: same stdout, stderr, exit code
ok=0; bad=0; EXPECT=ok
both() {  # a shared refusal passes if EXPECT=fail
  local l=$1; shift
  "$C" "$@" < "${IN:-/dev/null}" > "$W/c.out" 2>"$W/c.err"; local a=$?
  "$R" "$@" < "${IN:-/dev/null}" > "$W/r.out" 2>"$W/r.err"; local b=$?
  if cmp -s "$W/c.out" "$W/r.out" && cmp -s "$W/c.err" "$W/r.err" && [ $a = $b ] && { [ $a = 0 ] || [ $EXPECT = fail ]; }; then ok=$((ok+1)); tpass
  else bad=$((bad+1)); tfail "differs: $l (rc $a/$b)"; diff "$W/c.err" "$W/r.err" | head -4; fi
}
tally() { [ $bad = 0 ] && echo "$1: $ok forms identical" || echo "$1: $bad of $((ok+bad)) DIFFER"; ok=0; bad=0; }

# the --chunk fixture: a two-level hierarchy, one long sequence (sections smooth and edge)
chfix() {
python3 -c "
import random; random.seed(5); p=0
with open('$W/ch.bed','w') as f:
    for i in range(6000):
        l=random.randint(1,40); f.write('s1\t%d\t%d\t%s\n'%(p,p+l,random.choice('a b c d G1 G2 a a c novel'.split()))); p+=l
"
printf 'feature_set\tchild\tparent\nfs\tG1\tcategorized\nfs\tG2\tcategorized\nfs\ta\tG1\nfs\tb\tG1\nfs\tc\tG2\nfs\td\tG2\n' > "$W/ch.tsv"
"$C" compress < "$W/ch.bed" > "$W/ch.kst" 2>/dev/null
}

# 2. encode, every fixture and flag
if sec compress; then
for f in "$D"/*.bed; do
  for fl in "" "--block 4" "--block 262144" \
            "--effort low" "--effort low --block 4" "--names drop" "--effort high" \
            "--read-blocks 1" "--read-blocks 3" "--read-blocks 40" "--effort low --read-blocks 3"; do
    IN=$f both "compress $(basename "$f") [$fl]" $fl compress
  done
done
tally compress
fi

# 3. Rust decodes what C encoded
if sec decompress; then
dok=0; dbad=0
for f in "$D"/*.bed; do
  for fl in "" "--block 4" "--block 262144" "--read-blocks 3" "--read-blocks 40"; do
    "$C" $fl compress < "$f" > "$W/t.kst" 2>/dev/null
    "$R" decompress -i "$W/t.kst" 2>/dev/null | cmp -s - "$f" && { dok=$((dok+1)); tpass; } \
      || { dbad=$((dbad+1)); tfail "decompress differs: $(basename "$f") [$fl]"; }
  done
done
[ $dbad = 0 ] && echo "decompress: $dok forms give the input back" \
  || echo "decompress: $dbad of $((dok+dbad)) DIFFER"
fi

# 4. C decodes what Rust encoded
if sec cross; then
x0=$NFAIL
for f in "$D"/ls.bed "$D"/db.rall.bed; do
  "$R" compress < "$f" > "$W/r.kst" 2>/dev/null
  "$C" decompress -i "$W/r.kst" 2>/dev/null | cmp -s - "$f" \
    && tpass || tfail "the C cannot decompress the Rust's $(basename "$f")"
done
[ $NFAIL = $x0 ] && echo "cross-decompress: the C reads the Rust's files"
fi

# 5. region: same rows, every query shape
if sec region; then
for f in ls.bed s.bed db.rall.bed ord.bed; do
  for bl in "" "--block 4" "--read-blocks 3"; do
    "$C" $bl compress < "$D/$f" > "$W/t.kst" 2>/dev/null
    for seq in $(awk -F'\t' '{print $1}' "$D/$f" | uniq | head -3); do
      for q in "-r $seq" "-r $seq:50-51" "-r $seq:0-100" "-r $seq:999999-1000000"; do
        both "region $f [$bl] [$q]" query $q -i "$W/t.kst"
      done
    done
  done
done
tally region
fi

# 6. index: same trailer, cross-read
if sec index; then
"$C" compress < "$D/db.rall.bed" > "$W/ci.kst" 2>/dev/null
cp "$W/ci.kst" "$W/ri.kst"
"$C" index -i "$W/ci.kst" >/dev/null 2>&1; irc=$?
"$R" index -i "$W/ri.kst" >/dev/null 2>&1; irr=$?
if [ $irc = 0 ] && [ $irr = 0 ] && cmp -s "$W/ci.kst" "$W/ri.kst"; then tpass "index: byte-identical trailer"
else tfail "index: TRAILERS DIFFER (rc $irc/$irr)"; fi
both "index cross-read" query -i "$W/ci.kst" -r seq0005; tally "index cross-read"
fi

# 7. reader flags: dropped names, --raw
if sec reader-flags; then
for f in ls.bed s.bed db.rall.bed ord.bed; do
  "$C" compress --names drop < "$D/$f" > "$W/d.kst" 2>/dev/null
  cut -f1 "$D/$f" | uniq > "$W/nm"
  for dfl in "--names $W/nm" "--raw"; do both "decompress $dfl on $f" $dfl decompress -i "$W/d.kst"; done
  EXPECT=fail both "decompress without names on $f" decompress -i "$W/d.kst"; EXPECT=ok
  "$C" compress < "$D/$f" > "$W/p.kst" 2>/dev/null
  both "--raw $f" --raw decompress -i "$W/p.kst"
done
tally "reader flags"
fi

# 8. smooth over a real hierarchy shape
if sec smooth; then
for f in ls.bed s.bed db.rall.bed; do
  { printf 'feature_set\tchild\tparent\n'; cut -f4 "$D/$f" | sort -u | sed 's/^/fs\t/; s/$/\tcategorized/'; } > "$W/h.tsv"
  "$C" compress < "$D/$f" > "$W/t.kst" 2>/dev/null
  for q in "" "--v1" "--max-gap 0" "--max-gap 5" "--max-gap 100000"; do both "decompress $f [$q]" $q decompress -i "$W/t.kst" --hier "$W/h.tsv" --smooth --set fs; done
  for seq in $(awk -F'\t' '{print $1}' "$D/$f" | uniq | head -2); do
    for q in "-r $seq" "-r $seq:0-100"; do both "query $f [$q]" query -i "$W/t.kst" --hier "$W/h.tsv" --smooth --set fs $q; done
  done
  EXPECT=fail both "query --smooth without -r $f" query -i "$W/t.kst" --hier "$W/h.tsv" --smooth --set fs; EXPECT=ok
done
# --chunk: two-level hierarchy, one long sequence
chfix
for q in "--chunk 1" "--chunk 7" "--chunk 64" "--chunk 64 --v1" "--chunk 100000"; do both "decompress --chunk [$q]" decompress -i "$W/ch.kst" --hier "$W/ch.tsv" --smooth --set fs $q; done
both "decompress --smooth, set inferred" decompress -i "$W/ch.kst" --hier "$W/ch.tsv" --smooth
both "query --smooth -r, set inferred" query -i "$W/ch.kst" --hier "$W/ch.tsv" --smooth -r s1:10-20
both "query --chunk -r" query -i "$W/ch.kst" --hier "$W/ch.tsv" --smooth --set fs --chunk 64 -r s1:1000-2000
tally decompress
fi

# 9. --db, both directions.
if sec db; then
DB="-i $D/db.acrocentric.bed -i $D/db.chromosome.bed -i $D/db.gene.bed -i $D/db.region.bed -i $D/db.repeat.bed -i $D/db.subtelomeric.bed"
"$C" compress --db "$D/db.features.tsv" $DB > "$W/cdb.kst" 2>/dev/null
"$R" compress --db "$D/db.features.tsv" $DB > "$W/rdb.kst" 2>/dev/null
if cmp -s "$W/cdb.kst" "$W/rdb.kst"; then tpass "--db compress: identical"
else tfail "--db compress: DIFFERS"; fi
rm -rf "$W/co" "$W/ro"; mkdir -p "$W/co" "$W/ro"
"$C" decompress --db "$D/db.features.tsv" --out-dir "$W/co" -i "$W/cdb.kst" >/dev/null 2>&1
"$R" decompress --db "$D/db.features.tsv" --out-dir "$W/ro" -i "$W/cdb.kst" >/dev/null 2>&1
if diff -r "$W/co" "$W/ro" >/dev/null 2>&1; then tpass "--db decompress: identical"
else tfail "--db decompress: DIFFERS"; fi
# outside tuple with --hier: 2^26 + rank
{ printf 'feature_set\tchild\tparent\nrepeat\trepeat\tcategorized\nrepeat\tLINE\trepeat\n'; } > "$W/hv.tsv"
awk -F'\t' 'NR==5{$4="repeat"} {print}' OFS='\t' "$D/db.repeat.bed" > "$W/hv.repeat.bed"
HV=${DB/$D\/db.repeat.bed/$W\/hv.repeat.bed}
both "--db --hier" compress --db "$D/db.features.tsv" --hier "$W/hv.tsv" $HV; cp "$W/c.out" "$W/hv.kst"
{ printf 'feature_set\tchild\tparent\n'; for s in $(head -1 "$D/db.features.tsv" | cut -f2-); do cut -f4 "$D/db.$s.bed" | sort -u | grep -vx novel | sed "s/^/$s\t/; s/\$/\tcategorized/"; done; } > "$W/full.tsv"
"$C" --db "$D/db.features.tsv" compress $DB > "$W/fd.kst" 2>/dev/null
both "--db decompress --smooth --out-dir" --db "$D/db.features.tsv" --hier "$W/full.tsv" decompress -i "$W/fd.kst" --smooth --out-dir "$W/sod"
both "--db decompress --smooth --set gene,repeat --out-dir" --db "$D/db.features.tsv" --hier "$W/full.tsv" decompress -i "$W/fd.kst" --smooth --set gene,repeat --out-dir "$W/sod1"
EXPECT=fail both "--db --hier decompress --smooth --out-dir, a set not in the hierarchy" --db "$D/db.features.tsv" --hier "$W/hv.tsv" decompress -i "$W/hv.kst" --smooth --out-dir "$W/sod2"; EXPECT=ok
both "--db --hier decompress" --db "$D/db.features.tsv" --hier "$W/hv.tsv" --set repeat decompress -i "$W/hv.kst"; tally "--db --hier"
dn=0
for b in $DB; do [ "$b" = -i ] && continue;
  n=$(basename "$b" .bed); n=${n#db.}
  cmp -s "$W/ro/$n.bed" "$b" || dn=$((dn+1))
done
[ $dn = 0 ] && tpass "--db decompress: restores all six inputs" \
  || tfail "--db decompress: $dn inputs not restored"
fi

# 10. --aln, both directions, and aln-pack.
if sec aln; then
aok=0; abad=0
for pair in "s.bed p.paf ref.bed" "ls.bed lp.paf lref.bed"; do
  set -- $pair
  for fl in "" "--block 4" "--effort low"; do
    "$C" $fl compress --aln "$D/$2" --ref "$D/$3" < "$D/$1" > "$W/c.kst" 2>/dev/null; c1=$?
    "$R" $fl compress --aln "$D/$2" --ref "$D/$3" < "$D/$1" > "$W/r.kst" 2>/dev/null; c2=$?
    if cmp -s "$W/c.kst" "$W/r.kst" && [ $c1 = $c2 ]; then aok=$((aok+1)); tpass
    else abad=$((abad+1)); tfail "--aln compress differs: $1 [$fl]"; fi
  done
  "$C" compress --aln "$D/$2" --ref "$D/$3" < "$D/$1" > "$W/a.kst" 2>/dev/null
  "$C" decompress --aln "$D/$2" --ref "$D/$3" -i "$W/a.kst" > "$W/c.out" 2>/dev/null; c1=$?
  "$R" decompress --aln "$D/$2" --ref "$D/$3" -i "$W/a.kst" > "$W/r.out" 2>/dev/null; c2=$?
  if cmp -s "$W/c.out" "$W/r.out" && [ $c1 = $c2 ] && cmp -s "$W/r.out" "$D/$1"; then aok=$((aok+1)); tpass
  else abad=$((abad+1)); tfail "--aln decompress differs: $1"; fi
  # a .kst reference: same bytes and rows
  "$C" compress < "$D/$3" > "$W/c.ref.kst" 2>/dev/null; "$R" compress < "$D/$3" > "$W/r.ref.kst" 2>/dev/null
  "$C" compress --aln "$D/$2" --ref "$W/r.ref.kst" < "$D/$1" > "$W/c.kst" 2>/dev/null; c1=$?
  "$R" compress --aln "$D/$2" --ref "$W/c.ref.kst" < "$D/$1" > "$W/r.kst" 2>/dev/null; c2=$?
  "$R" decompress --aln "$D/$2" --ref "$W/c.ref.kst" -i "$W/a.kst" > "$W/r.out" 2>/dev/null; c3=$?
  if cmp -s "$W/c.ref.kst" "$W/r.ref.kst" && cmp -s "$W/c.kst" "$W/a.kst" && cmp -s "$W/r.kst" "$W/a.kst" && [ $c1 = 0 ] && [ $c2 = 0 ] && [ $c3 = 0 ] && cmp -s "$W/r.out" "$D/$1"; then aok=$((aok+1)); tpass
  else abad=$((abad+1)); tfail ".kst reference differs: $1"; fi
done
# one-block query: empty deltas stream
printf 's\t0\t1\tA\n' > "$W/one.ref.bed"; printf 'q\t1\t0\t1\t+\ts\t1\t0\t1\t1\t1\t60\tcg:Z:1M\n' > "$W/one.paf"
for pair in "p.paf ref.bed" "lp.paf lref.bed" "db.self.paf db.rall.bed" "$W/one.paf $W/one.ref.bed"; do
  set -- $pair; a=$1; b=$2; [ -e "$a" ] || a=$D/$1; [ -e "$b" ] || b=$D/$2
  "$C" aln-pack -a "$a" --ref "$b" -o "$W/c.pk" >/dev/null 2>"$W/c.pe"; c1=$?
  "$R" aln-pack -a "$a" --ref "$b" -o "$W/r.pk" >/dev/null 2>"$W/r.pe"; c2=$?
  if cmp -s "$W/c.pk" "$W/r.pk" && cmp -s "$W/c.pe" "$W/r.pe" && [ $c1 = $c2 ]; then aok=$((aok+1)); tpass
  else abad=$((abad+1)); tfail "aln-pack differs: $1"; fi
done
[ $abad = 0 ] && echo "--aln and aln-pack: $aok forms identical" \
  || echo "--aln: $abad DIFFER"
fi

# 11. same refusals, same words
if sec refusals; then
x0=$NFAIL
EXPECT=fail
for bad_in in 'r1\t0\t10' 'r1\t5\t10\tA' 'r1\t0\t0\tA' '\t0\t10\tA' 'r1\tx\t10\tA' 'r1\t00\t10\tA' 'r1\t0\t010\tA'; do
  printf "$bad_in\n" | "$C" compress 2> "$W/ce" >/dev/null; c1=$?; cm=$(tail -1 "$W/ce")
  printf "$bad_in\n" | "$R" compress 2> "$W/re" >/dev/null; r1=$?; rm2=$(tail -1 "$W/re")
  [ "$cm" = "$rm2" ] && [ "$c1" = 1 ] && [ "$r1" = 1 ] && tpass || tfail "refusal differs on [$bad_in]: c=[$cm] rc $c1, rust=[$rm2] rc $r1"
done
# refusals on companion files
printf 'q1\n\nq2\n' > "$W/blank"; gzip -c "$D/ref.bed" | head -c 40 > "$W/cut.gz"
printf '\ns\t0\t10\t1\n' > "$W/dbblank.bed"; printf 's\t00\t10\t1\n' > "$W/dbzero.bed"
"$C" aln-pack -a "$D/p.paf" --ref "$D/ref.bed" -o "$W/over.pk" >/dev/null 2>&1; printf '\0' >> "$W/over.pk"
printf 'KSTAL1abc' > "$W/short.pk"; printf 'KST' > "$W/tiny.paf"   # the magic is read before the body
"$C" compress --aln "$D/p.paf" --ref "$D/ref.bed" < "$D/s.bed" > "$W/a.kst" 2>/dev/null
for form in "compress --aln $D/p.paf --ref $D/ref.bed --aln-names $W/blank" "compress --aln $D/p.paf --ref $W/cut.gz" "decompress --aln $W/over.pk --ref $D/ref.bed -i $W/a.kst" \
            "decompress --aln $W/short.pk --ref $D/ref.bed -i $W/a.kst" "decompress --aln $W/tiny.paf --ref $D/ref.bed -i $W/a.kst" \
            "--db $D/db.features.tsv compress -i $W/dbblank.bed" "--db $D/db.features.tsv compress -i $W/dbzero.bed"; do
  "$C" $form < "$D/s.bed" 2> "$W/ce" >/dev/null; c1=$?; cm=$(tail -1 "$W/ce")
  "$R" $form < "$D/s.bed" 2> "$W/re" >/dev/null; r1=$?; rm2=$(tail -1 "$W/re")
  [ "$cm" = "$rm2" ] && [ "$c1" = 1 ] && [ "$r1" = 1 ] && tpass || tfail "refusal differs on [$form]: c=[$cm] rc $c1, rust=[$rm2] rc $r1"
done
[ $NFAIL = $x0 ] && echo "refusals: same message from both"
EXPECT=ok
fi

# 12. awkward valid input: same bytes
if sec awkward; then
python3 "$D/kt.py" gen "$W" 20
gok=0; gbad=0
for f in "$W"/ok*.bed; do
  for fl in "" "--block 4" "--block 262144" "--names drop" "--effort low --block 3"; do
    "$C" $fl compress < "$f" > "$W/c.kst" 2>/dev/null; c1=$?
    "$R" $fl compress < "$f" > "$W/r.kst" 2>/dev/null; c2=$?
    if cmp -s "$W/c.kst" "$W/r.kst" && [ $c1 = $c2 ]; then gok=$((gok+1)); tpass
    else gbad=$((gbad+1)); tfail "awkward input differs: $(basename "$f") [$fl]"; fi
  done
done
[ $gbad = 0 ] && echo "awkward input: $gok forms identical" || echo "awkward input: $gbad DIFFER"
fi

EXPECT=ok
# 11b. gzip input: same track as plain
if sec gzip; then
gok=0; gbad=0
for f in "$D"/s.bed "$D"/ls.bed "$D"/db.rall.bed "$D"/ord.bed; do
  gzip -c "$f" > "$W/in.gz"
  "$C" compress < "$f" > "$W/c.plain" 2>/dev/null
  for b in "$C" "$R"; do
    "$b" compress < "$W/in.gz" > "$W/g.kst" 2>/dev/null
    if cmp -s "$W/c.plain" "$W/g.kst"; then gok=$((gok+1)); tpass
    else gbad=$((gbad+1)); tfail "gz differs: $(basename "$f") $b"; fi
  done
done
gzip -c "$D/ref.bed" > "$W/ref.gz"; gzip -c "$D/p.paf" > "$W/paf.gz"
gzip -c "$D/s.bed" > "$W/s.gz"
"$C" compress --aln "$D/p.paf" --ref "$D/ref.bed" < "$D/s.bed" > "$W/a.plain" 2>/dev/null
for b in "$C" "$R"; do
  "$b" compress --aln "$W/paf.gz" --ref "$W/ref.gz" < "$W/s.gz" > "$W/a.gz.kst" 2>/dev/null
  if cmp -s "$W/a.plain" "$W/a.gz.kst"; then gok=$((gok+1)); tpass
  else gbad=$((gbad+1)); tfail "gz --aln differs: $b"; fi
done
# many BGZF blocks, file and pipe
awk 'BEGIN{for(i=0;i<60000;i++){p=0;for(j=0;j<5;j++){printf "read%d\t%d\t%d\tL%d\n",i,p,p+j+7,j%3;p+=j+7}}}' > "$W/bz.bed"
python3 "$D/kt.py" bgzf "$W/bz.bed" "$W/bz.bgz"
"$C" compress < "$W/bz.bed" > "$W/bz.plain" 2>/dev/null
for b in "$C" "$R"; do
  for how in file pipe; do
    if [ $how = file ]; then "$b" compress -i "$W/bz.bgz" > "$W/bz.kst" 2>/dev/null; else "$b" compress < "$W/bz.bgz" > "$W/bz.kst" 2>/dev/null; fi
    if cmp -s "$W/bz.plain" "$W/bz.kst"; then gok=$((gok+1)); tpass; else gbad=$((gbad+1)); tfail "bgzf differs: $b $how"; fi
  done
done
sz=$(stat -c%s "$W/bz.bgz"); head -c $((sz / 2)) "$W/bz.bgz" > "$W/bzcut.bgz"   # cut and byte-flipped BGZF
cp "$W/bz.bgz" "$W/bzflip.bgz"; printf 'U' | dd of="$W/bzflip.bgz" bs=1 seek=$((sz / 3)) conv=notrunc status=none
for f in bzcut bzflip; do
  "$C" compress -i "$W/$f.bgz" 2> "$W/ce" >/dev/null; c1=$?; cm=$(tail -1 "$W/ce")
  "$R" compress -i "$W/$f.bgz" 2> "$W/re" >/dev/null; r1=$?; rm2=$(tail -1 "$W/re")
  [ "$cm" = "$rm2" ] && [ "$c1" = 1 ] && [ "$r1" = 1 ] && { gok=$((gok+1)); tpass; } || { gbad=$((gbad+1)); tfail "$f refusal differs: c=[$cm] rc $c1, rust=[$rm2] rc $r1"; }
done
# forged BGZF: ISIZE bounds, BC after another subfield, cuts inside the extra field; same bytes or same refusal
mkdir -p "$W/bx"; python3 "$D/kt.py" bgzfx "$W/bx"
for f in "$W"/bx/*.gz; do
  for how in file pipe; do
    if [ $how = file ]; then "$C" compress -i "$f" > "$W/bxc" 2> "$W/ce"; c1=$?; "$R" compress -i "$f" > "$W/bxr" 2> "$W/re"; r1=$?
    else "$C" compress < "$f" > "$W/bxc" 2> "$W/ce"; c1=$?; "$R" compress < "$f" > "$W/bxr" 2> "$W/re"; r1=$?; fi
    cmp -s "$W/bxc" "$W/bxr" && cmp -s "$W/ce" "$W/re" && [ $c1 = $r1 ] && { gok=$((gok+1)); tpass; } \
      || { gbad=$((gbad+1)); tfail "forged bgzf differs: $(basename "$f") $how, c=[$(tail -1 "$W/ce")] rc $c1, rust=[$(tail -1 "$W/re")] rc $r1"; }
  done
done
[ $gbad = 0 ] && echo "gzip input: $gok forms match the plain track" \
  || echo "gzip input: $gbad DIFFER"
fi

# 12. CLI: same words, same exit code
if sec cli; then
cok=0; cbad=0
cli() {
  (cd "$W" && "$CB" "$@" </dev/null > cliout.c 2>&1); local cr=$?   # in $W, since forms may write files
  (cd "$W" && "$RB" "$@" </dev/null > cliout.r 2>&1); local rr=$?
  if cmp -s "$W/cliout.c" "$W/cliout.r" && [ "$cr" = "$rr" ]; then cok=$((cok+1)); tpass; return; fi
  cbad=$((cbad+1)); tfail "cli differs: [$*] rc $cr/$rr"; echo "    c=[$(head -c 200 "$W/cliout.c" | tr -d '\0' | head -1)]"
  echo "    rust=[$(head -c 200 "$W/cliout.r" | tr -d '\0' | head -1)]"
}
cli; cli --help; cli -h; cli --version; cli --bogus; cli --bogus compress
cli decompress; cli compress extra; cli --profile; cli --threads; cli --threads 4
cli badcmd; cli region; cli decompress -i f; cli index; cli decompress -i f; cli aln-pack
cli compress --names ordinal; cli index --effort low
cli compress --block; cli compress --effort max; cli compress --block 0
cli -t; cli -t 0; cli -t 99999; cli --tmp-dir; cli decompress --raw x --raw
cli index --block 0; cli aln-pack -a a --ref b; cli decompress -i f --hier h; cli query -i f -r s:5-5
cli compress --names foo; cli --effort low; cli -r x
cli aln-pack --names foo -a a --ref b -o c; cli compress --db; cli --version compress
cli index -i f x; cli decompress -i f x; cli query -i f --hier h --smooth --set fs -r s:9-3; cli aln-pack -a a --ref b -o c d
cli query -i f --hier h --smooth --set fs -r s:1-2; cli query -i f -r s --halo 3; cli query -i f -r s:1-2
cli decompress --max-gap x -i f --hier h --smooth --set fs; cli query -i f --hier h --smooth --set fs; cli query -i f --smooth --out-dir d
cli decompress --out-dir D; cli decompress -i f --halo 2; cli query -i f -r s:+1-5; cli query -i f -r s:99999999999999999999-1
cli compress --block +5; cli compress --aln p; cli compress --ref r
cli --db $D/db.features.tsv compress -i nosuch.bed; cli --db $D/db.features.tsv compress -i $D/db.gene.bed -i nosuch.bed -i x -i y -i z -i w
cli decompress --aln-names f; cli index --names ordinal
cli decompress --aln p; cli decompress --ref r; cli decompress --aln p --ref r
cli compress --names drop --aln p --ref r; cli decompress --db f; cli compress --db f
cli decompress -i a -i b; cli compress -i nosuch; cli index -i f -o x; cli query -i f -r s:2-1
cli query -i f -r s --out-dir d; cli query -i f -r s --raw; cli decompress -i f --halo 2; cli region -i f -r s
cli --tmp-dir /nonexistent compress; cli --tmp-dir /nonexistent --db $D/db.features.tsv compress -i $D/db.gene.bed
cli --block=4 compress; cli -b4 compress; cli -t4 compress; cli --effort=low compress; cli -V; cli foo --raw
cli decompress -i f --max-gap -1; cli query -i f -r -x; cli compress -- -x; cli -o -i compress
cli decompress --db f --out-dir d -o x -i f; cli decompress --db f --out-dir "" -i f; cli decompress -i f --hier h --smooth --out-dir ""
ln -s nowhere/x "$W/dang"; cli compress -o dang; cli compress -o nosuchdir/x; cli aln-pack -a a --ref b -o nosuchdir/c
[ $cbad = 0 ] && echo "cli: $cok forms give the same words and code" \
  || echo "cli: $cbad of $((cok+cbad)) DIFFER"
fi

# 13. tail predictor, past TAILSKIP=64 sequences
if sec tail; then
python3 - "$W" <<'PY'
import sys, random
random.seed(5); w = sys.argv[1]
for name, lens in (("tail1", [100]), ("tail2", [100, 150])):
    with open(f"{w}/{name}.bed", "w") as f:
        for r in range(200):
            tot = random.choice(lens); cuts = sorted(random.sample(range(1, tot), random.randint(0, 3)))
            p = 0
            for k, e in enumerate(cuts + [tot]): f.write(f"q{r}\t{p}\t{e}\t{'ABCD'[(r + k) % 4]}\n"); p = e
PY
tok=0; tbad=0
for f in tail1 tail2; do for flags in "" "-e low" "-e high" "--block 7"; do
  "$C" $flags compress -i "$W/$f.bed" > "$W/c.kst" 2>/dev/null; "$R" $flags compress -i "$W/$f.bed" > "$W/r.kst" 2>/dev/null
  if cmp -s "$W/c.kst" "$W/r.kst" && "$R" decompress -i "$W/c.kst" 2>/dev/null | cmp -s - "$W/$f.bed"; then tok=$((tok+1)); tpass; else tbad=$((tbad+1)); tfail "tail differs: $f $flags"; fi
done; done
[ $tbad = 0 ] && echo "tail predictor: $tok forms identical" || echo "tail predictor: $tbad DIFFER"
fi

# 13b. -e high's zpaq frames, both ways
if sec zpaq; then
awk 'BEGIN { srand(5); for (c = 0; c < 6; c++) { g = c % 4; nr = 20000 + int(rand() * 20000); p = 0; for (r = 0; r < nr; r++) { l = 1 + int(rand() * (g ? 9 * g : 37)); u = rand(); printf "ctg%d\t%d\t%d\t%s\n", c, p, p + l, u < 0.5 + g / 10 ? "chr" (1 + g) : "novel" int(rand() * 3); p += l } } }' > "$W/zt.bed"
"$C" -e high compress -i "$W/zt.bed" > "$W/c.kst" 2>/dev/null; "$R" -e high compress -i "$W/zt.bed" > "$W/r.kst" 2>/dev/null
if cmp -s "$W/c.kst" "$W/r.kst" && "$R" decompress -i "$W/c.kst" 2>/dev/null | cmp -s - "$W/zt.bed" && "$C" decompress -i "$W/r.kst" 2>/dev/null | cmp -s - "$W/zt.bed"
then tpass "zpaq frames: identical, each reads the other"; else tfail "zpaq frames: DIFFER"; fi
awk 'BEGIN { srand(7); for (r = 0; r < 60000; r++) { k = rand() < 0.4 ? 1 : 2 + int(rand() * rand() * 12); p = 0; g = int(rand() * 4); for (j = 0; j < k; j++) { l = j == k - 1 ? 120 - p : 1 + int(rand() * (120 - p - (k - j - 1)) / (k - j)); printf "r%d\t%d\t%d\t%s\n", r, p, p + l, "chr" (1 + (g + j) % 5); p += l } } }' > "$W/zi.bed"
"$C" -e high compress -i "$W/zi.bed" > "$W/c.kst" 2>/dev/null; "$R" -e high compress -i "$W/zi.bed" > "$W/r.kst" 2>/dev/null
if cmp -s "$W/c.kst" "$W/r.kst" && "$R" decompress -i "$W/c.kst" 2>/dev/null | cmp -s - "$W/zi.bed" && "$C" decompress -i "$W/r.kst" 2>/dev/null | cmp -s - "$W/zi.bed"
then tpass "interleaved zpaq model: identical, each reads the other"; else tfail "interleaved zpaq model: DIFFER"; fi
fi

# 14. damaged file: same rows before refusal
if sec damaged; then
python3 -c "
import random; random.seed(9)
with open('$W/many.bed','w') as f:
    for r in range(400000): f.write(f'read{r:07d}\t0\t{random.randint(50,150)}\t{random.choice(\"ABCD\")}\n')"
"$C" -e low --block 4096 compress -i "$W/many.bed" > "$W/many.kst" 2>/dev/null
dok=0; dbad=0
for pct in 75 90; do
  python3 -c "import sys; d=bytearray(open(sys.argv[1],'rb').read()); d[len(d)*$pct//100]^=0xff; open(sys.argv[2],'wb').write(d)" "$W/many.kst" "$W/cut.kst"
  "$C" decompress -i "$W/cut.kst" > "$W/c.out" 2>"$W/c.err"; c1=$?; "$R" decompress -i "$W/cut.kst" > "$W/r.out" 2>"$W/r.err"; r1=$?
  if cmp -s "$W/c.out" "$W/r.out" && cmp -s "$W/c.err" "$W/r.err" && [ $c1 = $r1 ] && [ $c1 != 0 ] && [ -s "$W/c.out" ]; then dok=$((dok+1)); tpass; else dbad=$((dbad+1)); tfail "damaged decompress differs at $pct%: $(wc -c < "$W/c.out") vs $(wc -c < "$W/r.out") B, rc $c1/$r1"; fi
done
[ $dbad = 0 ] && echo "damaged decompress: $dok flips stop at the same byte" || echo "damaged decompress: $dbad DIFFER"
fi

# 15. forged chunk tables, re-signed
if sec forged; then
python3 -c "
with open('$W/ck.bed','w') as f:
    for i in range(65537): f.write(f'r{i}\t0\t1\t{\"B\" if i & 1 else \"A\"}\n')"
"$C" -e low --block 1000 compress -i "$W/ck.bed" > "$W/ck.kst" 2>/dev/null
python3 - "$W" "$D" <<'PY'
import sys
sys.path.insert(0, sys.argv[2])
import kt
w = sys.argv[1]
src, h = kt.load(f"{w}/ck.kst")
off = 136 + h[6] + (h[1] & ~kt.FLAGS) + h[2] + h[3]
tab = kt.unzst(src[off:off + h[11]]); pairs = []; p = 0
while p < len(tab):
    cl, p = kt.vget(tab, p); sp, p = kt.vget(tab, p); pairs.append([cl, sp])
for mode in ("wrap", "shift"):
    q = [x[:] for x in pairs]
    if mode == "wrap": q[0][0], q[1][0] = (1 << 64) - 1, h[2] + 1
    else: q[0][1], q[1][1] = 1, 65537
    z = kt.zst(w, b"".join(kt.vput(a) + kt.vput(b) for a, b in q))
    out = bytearray(src); out[off:off + h[11]] = z; hh = h[:]; hh[11] = len(z); kt.sign(out, hh)
    open(f"{w}/{mode}.kst", "wb").write(out)
PY
EXPECT=fail
for m in wrap shift; do both "forged chunk table $m" query -i "$W/$m.kst" -r r0; done
EXPECT=ok
tally "forged chunk tables"
fi

# 16. edge cases: empty input, a name holding ':', --tmp-dir on readers, --max-gap 0
if sec edge; then
IN=/dev/null both "compress of empty input" compress; cp "$W/c.out" "$W/empty.kst"
: > "$W/empty.bed"; gzip -c < /dev/null > "$W/empty.gz"
both "compress -i an empty file" compress -i "$W/empty.bed"
both "compress -i an empty gzip" compress -i "$W/empty.gz"
both "compress --names drop of empty input" compress --names drop
both "compress --read-blocks 3 of empty input" --read-blocks 3 compress
both "decompress of an empty track" decompress -i "$W/empty.kst"
both "decompress --raw of an empty track" --raw decompress -i "$W/empty.kst"
EXPECT=fail both "query of an empty track" query -i "$W/empty.kst" -r x
printf 'HLA:A*01\t0\t10\tA\nHLA:A*01\t10\t20\tB\nx:5-6\t0\t10\tC\n' > "$W/colon.bed"
"$C" compress < "$W/colon.bed" > "$W/colon.kst" 2>/dev/null
for q in 'HLA:A*01' 'HLA:A*01:12-13' 'x:5-6:0-3'; do both "query a name holding ':' [-r $q]" query -i "$W/colon.kst" -r "$q"; done
EXPECT=fail both "query -r x:5-6, a range on a missing x" query -i "$W/colon.kst" -r x:5-6; EXPECT=ok
"$C" compress < "$D/s.bed" > "$W/td.kst" 2>/dev/null
for t in /nonexistent "$W"; do
  both "--tmp-dir $t decompress" --tmp-dir "$t" decompress -i "$W/td.kst"
  both "--tmp-dir $t query" --tmp-dir "$t" query -i "$W/td.kst" -r q1
done
[ -s "$W/ch.kst" ] || chfix
for q in "--max-gap 0" "--max-gap 0 --chunk 64" "--max-gap 0 --v1"; do both "decompress --smooth [$q]" decompress -i "$W/ch.kst" --hier "$W/ch.tsv" --smooth --set fs $q; done
tally edge
fi


# 17. a body temp file that fills: the staging fits, the one-row blocks' frames do not; the Rust used to hang here
if sec bodytmp && stress "full body temp file at 2 threads: both refuse alike"; then
awk 'BEGIN{p=0; for(k=0;k<20000;k++){printf "r1\t%d\t%d\tchr%d\n", p, p+20, k%22; p+=20}}' > "$W/body.bed"
mkdir -p "$W/btmp"
for x in C R; do
  rc=0; ( ulimit -f 64; trap "" XFSZ; timeout 60 "${!x}" -t 2 --tmp-dir "$W/btmp" --block 1 compress < "$W/body.bed" >/dev/null 2>"$W/$x.berr" ) || rc=$?
  eval "rc$x=$rc"
done
if [ $rcC = 1 ] && [ $rcR = 1 ] && cmp -s "$W/C.berr" "$W/R.berr" && [ -z "$(ls -A "$W/btmp")" ]; then tpass "full body temp file: both refuse alike, $(cat "$W/R.berr" | tail -1)"
else tfail "full body temp file: DIFFER (rc $rcC/$rcR)"; diff "$W/C.berr" "$W/R.berr" | head -4; ls "$W/btmp"; fi
fi

# 18. forged inputs the readers once told apart: same words, or the same rows
if sec forged2; then
"$C" compress < "$D/s.bed" > "$W/pl.kst" 2>/dev/null
"$C" aln-pack -a "$D/p.paf" --ref "$D/ref.bed" -o "$W/fg.pk" >/dev/null 2>&1
head -2 "$D/p.paf" > "$W/p2.paf"   # one gapless block a query: the deltas stream is empty
"$C" aln-pack -a "$W/p2.paf" --ref "$D/ref.bed" -o "$W/fg2.pk" >/dev/null 2>&1
python3 - "$D" "$W" <<'PY'
import sys
sys.path.insert(0, sys.argv[1]); import kt
W = sys.argv[2]
for bit in (53, 54):   # a flag bit no reader knows
    d = bytearray(open(W + "/pl.kst", "rb").read()); d[16 + bit // 8] |= 1 << (bit % 8)
    open(W + "/fl%d.kst" % bit, "wb").write(d)
def secs(p):
    d = open(p, "rb").read(); body = d[14:]; i = 0; sec = []
    for _ in range(10):
        cd, i = kt.vget(body, i); rl, i = kt.vget(body, i); cl, i = kt.vget(body, i); sec.append((cd, rl, body[i:i + cl])); i += cl
    return d[:14], sec
pack = lambda hd, s: hd + b"".join(kt.vput(c) + kt.vput(r) + kt.vput(len(x)) + x for c, r, x in s)
hd, sec = secs(W + "/fg.pk"); open(W + "/big.pk", "wb").write(pack(hd, [(1, 1 << 40, b"")] + sec[1:]))   # a 1 TB stream
hd, sec = secs(W + "/fg2.pk"); k = [j for j, x in enumerate(sec) if x[0] == 0]
if k: sec[k[0]] = (0, 0, b"junk")   # bytes behind an empty stream
open(W + "/skip.pk", "wb").write(pack(hd, sec))
PY
"$C" compress --aln "$D/p.paf" --ref "$D/ref.bed" < "$D/s.bed" > "$W/fa.kst" 2>/dev/null
"$C" compress --aln "$W/fg2.pk" --ref "$D/ref.bed" < "$D/s.bed" > "$W/fa2.kst" 2>/dev/null
EXPECT=fail
for bit in 53 54; do both "header flag bit $bit" decompress -i "$W/fl$bit.kst"; done
both "a packed alignment claiming a 1 TB stream" decompress --aln "$W/big.pk" --ref "$D/ref.bed" -i "$W/fa.kst"
EXPECT=ok
both "a packed alignment with bytes behind an empty stream" decompress --aln "$W/skip.pk" --ref "$D/ref.bed" -i "$W/fa2.kst"
python3 - "$W" <<'PY'
import sys
W = sys.argv[1]; n = 65
open(W + "/f65.tsv", "w").write("featureID\t" + "\t".join("s%d" % k for k in range(n)) + "\n"
                                 + "".join("%d\t" % i + "\t".join("L%d" % ((i + k) % 3) for k in range(n)) + "\n" for i in range(1, 4)))
for k in range(n): open(W + "/f65.%d.bed" % k, "w").write("".join("c\t%d\t%d\tL%d\n" % (10 * j, 10 * j + 10, (j % 3 + k) % 3) for j in range(6)))
PY
in65=""; for k in $(seq 0 64); do in65="$in65 -i $W/f65.$k.bed"; done
both "--db with 65 -i" --db "$W/f65.tsv" compress $in65
tally "forged headers, packed alignments, 65 inputs"
fi

# 19. a ten-byte varint past 64 bits in the index table, re-signed: once read as its low 64 bits, now refused alike
if sec varint; then
"$C" --block 4 compress < "$D/s.bed" > "$W/v10.kst" 2>/dev/null
for t in 0x00 0x01 0x02 0x40 0x7f 0x80 0x81 0xff; do python3 "$D/kt.py" varint10 "$W/v10.kst" "$W/v10.$t.kst" $t "$W"; done
EXPECT=fail
for t in 0x02 0x40 0x7f 0x80 0x81 0xff; do
  both "a ten-byte varint ending $t, decompress" decompress -i "$W/v10.$t.kst"
  if grep -qx 'karyozip: corrupt varint' "$W/c.err"; then tpass; else tfail "a ten-byte varint ending $t: C says $(head -c 200 "$W/c.err")"; fi
  both "a ten-byte varint ending $t, query" query -i "$W/v10.$t.kst" -r q1
done
both "a ten-byte varint of 2^63 where 0 belongs" decompress -i "$W/v10.0x01.kst"   # it fits, and a later check refuses it
EXPECT=ok
both "a ten-byte varint of 0, not the shortest but it fits" decompress -i "$W/v10.0x00.kst"
if cmp -s "$W/c.out" "$D/s.bed"; then tpass; else tfail "a ten-byte varint of 0: rows WRONG"; fi
tally "ten-byte varints"
fi

# 20. a staged output whose fsync or close fails (tests/failclose.c, LD_PRELOAD): the same refusal, every old output kept
if sec finalize; then
FC=$W/failclose.so
if "${CC:-cc}" -shared -fPIC -O2 -o "$FC" "$src/failclose.c" -ldl 2>"$W/fc.err"; then
  fin=""; for f in acrocentric chromosome gene region repeat subtelomeric; do fin="$fin -i $D/db.$f.bed"; done
  "$C" compress --db "$D/db.features.tsv" $fin > "$W/fz.kt" 2>/dev/null; "$C" compress < "$D/s.bed" > "$W/fz.p.kt" 2>/dev/null
  { printf 'feature_set\tchild\tparent\n'; for s in $(head -1 "$D/db.features.tsv" | cut -f2-); do cut -f4 "$D/db.$s.bed" | sort -u | grep -vx novel | sed "s/^/$s\t/; s/\$/\tcategorized/"; done; } > "$W/fz.hier"
  mkdir -p "$W/fz/od" "$W/fz/sm"; for s in gene repeat; do echo "OLD_$s" > "$W/fz/od/$s.bed"; echo "OLD_$s" > "$W/fz/sm/$s.bed"; done
  echo old > "$W/fz/k"; echo old > "$W/fz/pk"; echo oldidx > "$W/fz.p.kt.kzi"
  EXPECT=fail
  for op in sync close; do
    KZ_FAIL_FD_OP=$op KZ_FAIL_FD_MATCH=.bed.tmp. LD_PRELOAD=$FC both "failed $op, decompress --db --out-dir" \
      decompress --db "$D/db.features.tsv" --set gene,repeat --out-dir "$W/fz/od" -i "$W/fz.kt"
    KZ_FAIL_FD_OP=$op KZ_FAIL_FD_MATCH=.bed.tmp. LD_PRELOAD=$FC both "failed $op, --smooth --out-dir" \
      --db "$D/db.features.tsv" --set gene,repeat decompress -i "$W/fz.kt" --hier "$W/fz.hier" --smooth --out-dir "$W/fz/sm"
    KZ_FAIL_FD_OP=$op KZ_FAIL_FD_MATCH=.kzi. LD_PRELOAD=$FC both "failed $op, index" index -i "$W/fz.p.kt"
    KZ_FAIL_FD_OP=$op KZ_FAIL_FD_MATCH=fz/k.tmp. LD_PRELOAD=$FC both "failed $op, -o" compress -i "$D/s.bed" -o "$W/fz/k"
    KZ_FAIL_FD_OP=$op KZ_FAIL_FD_MATCH=fz/pk.tmp. LD_PRELOAD=$FC both "failed $op, aln-pack" aln-pack -a "$D/p.paf" --ref "$D/ref.bed" -o "$W/fz/pk"
  done
  EXPECT=ok
  grep -q '^karyozip: write failed$' "$W/r.err" && [ "$(cat "$W"/fz/od/*.bed "$W"/fz/sm/*.bed "$W/fz/k" "$W/fz/pk" "$W/fz.p.kt.kzi" | tr '\n' ' ')" = "OLD_gene OLD_repeat OLD_gene OLD_repeat old old oldidx " ] \
    && [ "$(ls -A "$W/fz" "$W/fz/od" "$W/fz/sm" | tr '\n' ' ')" = "$W/fz: k od pk sm  $W/fz/od: gene.bed repeat.bed  $W/fz/sm: gene.bed repeat.bed " ] \
    && tpass || tfail "failed sync or close: an old output DAMAGED or a temporary left: $(ls -A "$W/fz" "$W/fz/od" "$W/fz/sm" | tr '\n' ' ')"
else tfail "tests/failclose.c: does not build with ${CC:-cc}: $(head -c 200 "$W/fc.err")"; fi
tally "sync and close failures"
fi

summary diff.sh
[ $fail = 0 ] && echo "DIFF_PASS" || echo "DIFF_FAIL"
exit $fail
