#!/bin/bash
# C and Rust: same bytes, same words
set -u
C=${1:-$(dirname "$0")/../c/build/karyozip}
R=${2:-$(dirname "$0")/../rust/target/release/karyozip}
D=$(cd "$(dirname "$0")" && pwd)
[ -x "$C" ] && [ -x "$R" ] || { echo "build both first: make all rust"; exit 1; }
C=$(realpath "$C"); R=$(realpath "$R")
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
fail=0

cv=$("$C" --version | awk '/zstd|libbz3|liblzma/{print $1"="$2}' | tr '\n' ' ')
rv=$("$R" --version | awk '/zstd|libbz3|liblzma/{print $1"="$2}' | tr '\n' ' ')
if [ "$cv" = "$rv" ]; then echo "codec versions: $cv"
else echo "codec versions DIFFER: c=[$cv] rust=[$rv]"; fail=1; fi

# both: same stdout, stderr, exit code
ok=0; bad=0; EXPECT=ok
both() {  # a shared refusal passes if EXPECT=fail
  local l=$1; shift
  "$C" "$@" < "${IN:-/dev/null}" > "$W/c.out" 2>"$W/c.err"; local a=$?
  "$R" "$@" < "${IN:-/dev/null}" > "$W/r.out" 2>"$W/r.err"; local b=$?
  if cmp -s "$W/c.out" "$W/r.out" && cmp -s "$W/c.err" "$W/r.err" && [ $a = $b ] && { [ $a = 0 ] || [ $EXPECT = fail ]; }; then ok=$((ok+1))
  else bad=$((bad+1)); echo "  differs: $l (rc $a/$b)"; diff "$W/c.err" "$W/r.err" | head -4; fi
}
tally() { [ $bad = 0 ] && echo "$1: $ok forms identical" || { echo "$1: $bad of $((ok+bad)) DIFFER"; fail=1; }; ok=0; bad=0; }

# 2. encode, every fixture and flag
for f in "$D"/*.bed; do
  for fl in "" "--block 4" "--block 262144" \
            "--effort low" "--effort low --block 4" "--names drop" "--effort high" \
            "--read-blocks 1" "--read-blocks 3" "--read-blocks 40" "--effort low --read-blocks 3"; do
    IN=$f both "compress $(basename "$f") [$fl]" $fl compress
  done
done
tally compress

# 3. Rust decodes what C encoded
dok=0; dbad=0
for f in "$D"/*.bed; do
  for fl in "" "--block 4" "--block 262144" "--read-blocks 3" "--read-blocks 40"; do
    "$C" $fl compress < "$f" > "$W/t.kst" 2>/dev/null
    "$R" decompress -i "$W/t.kst" 2>/dev/null | cmp -s - "$f" && dok=$((dok+1)) \
      || { dbad=$((dbad+1)); echo "  decompress differs: $(basename "$f") [$fl]"; }
  done
done
[ $dbad = 0 ] && echo "decompress: $dok forms give the input back" \
  || { echo "decompress: $dbad of $((dok+dbad)) DIFFER"; fail=1; }

# 4. C decodes what Rust encoded
for f in "$D"/ls.bed "$D"/db.rall.bed; do
  "$R" compress < "$f" > "$W/r.kst" 2>/dev/null
  "$C" decompress -i "$W/r.kst" 2>/dev/null | cmp -s - "$f" \
    || { echo "  the C cannot decompress the Rust's $(basename "$f")"; fail=1; }
done
echo "cross-decompress: the C reads the Rust's files"

# 5. region: same rows, every query shape
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

# 6. index: same trailer, cross-read
"$C" compress < "$D/db.rall.bed" > "$W/ci.kst" 2>/dev/null
cp "$W/ci.kst" "$W/ri.kst"
"$C" index -i "$W/ci.kst" >/dev/null 2>&1; irc=$?
"$R" index -i "$W/ri.kst" >/dev/null 2>&1; irr=$?
if [ $irc = 0 ] && [ $irr = 0 ] && cmp -s "$W/ci.kst" "$W/ri.kst"; then echo "index: byte-identical trailer"
else echo "index: TRAILERS DIFFER (rc $irc/$irr)"; fail=1; fi
both "index cross-read" query -i "$W/ci.kst" -r seq0005; tally "index cross-read"

# 7. reader flags: dropped names, --raw
for f in ls.bed s.bed db.rall.bed ord.bed; do
  "$C" compress --names drop < "$D/$f" > "$W/d.kst" 2>/dev/null
  cut -f1 "$D/$f" | uniq > "$W/nm"
  for dfl in "--names $W/nm" "--raw"; do both "decompress $dfl on $f" $dfl decompress -i "$W/d.kst"; done
  EXPECT=fail both "decompress without names on $f" decompress -i "$W/d.kst"; EXPECT=ok
  "$C" compress < "$D/$f" > "$W/p.kst" 2>/dev/null
  both "--raw $f" --raw decompress -i "$W/p.kst"
done
tally "reader flags"

# 8. smooth over a real hierarchy shape
for f in ls.bed s.bed db.rall.bed; do
  { printf 'feature_set\tchild\tparent\n'; cut -f4 "$D/$f" | sort -u | sed 's/^/fs\t/; s/$/\tcategorized/'; } > "$W/h.tsv"
  "$C" compress < "$D/$f" > "$W/t.kst" 2>/dev/null
  for q in "" "--v1" "--max-gap 5" "--max-gap 100000"; do both "decompress $f [$q]" $q decompress -i "$W/t.kst" --hier "$W/h.tsv" --smooth --set fs; done
  for seq in $(awk -F'\t' '{print $1}' "$D/$f" | uniq | head -2); do
    for q in "-r $seq" "-r $seq:0-100"; do both "query $f [$q]" query -i "$W/t.kst" --hier "$W/h.tsv" --smooth --set fs $q; done
  done
done
# --chunk on a two-level hierarchy, one long sequence
python3 -c "
import random; random.seed(5); p=0
with open('$W/ch.bed','w') as f:
    for i in range(6000):
        l=random.randint(1,40); f.write('s1\t%d\t%d\t%s\n'%(p,p+l,random.choice('a b c d G1 G2 a a c novel'.split()))); p+=l
"
printf 'feature_set\tchild\tparent\nfs\tG1\tcategorized\nfs\tG2\tcategorized\nfs\ta\tG1\nfs\tb\tG1\nfs\tc\tG2\nfs\td\tG2\n' > "$W/ch.tsv"
"$C" compress < "$W/ch.bed" > "$W/ch.kst" 2>/dev/null
for q in "--chunk 1" "--chunk 7" "--chunk 64" "--chunk 64 --v1" "--chunk 100000"; do both "decompress --chunk [$q]" decompress -i "$W/ch.kst" --hier "$W/ch.tsv" --smooth --set fs $q; done
both "decompress --smooth, set inferred" decompress -i "$W/ch.kst" --hier "$W/ch.tsv" --smooth
both "query --smooth -r, set inferred" query -i "$W/ch.kst" --hier "$W/ch.tsv" --smooth -r s1:10-20
both "query --chunk -r" query -i "$W/ch.kst" --hier "$W/ch.tsv" --smooth --set fs --chunk 64 -r s1:1000-2000
tally decompress

# 9. --db, both directions.
DB="-i $D/db.acrocentric.bed -i $D/db.chromosome.bed -i $D/db.gene.bed -i $D/db.region.bed -i $D/db.repeat.bed -i $D/db.subtelomeric.bed"
"$C" compress --db "$D/db.features.tsv" $DB > "$W/cdb.kst" 2>/dev/null
"$R" compress --db "$D/db.features.tsv" $DB > "$W/rdb.kst" 2>/dev/null
if cmp -s "$W/cdb.kst" "$W/rdb.kst"; then echo "--db compress: identical"
else echo "--db compress: DIFFERS"; fail=1; fi
rm -rf "$W/co" "$W/ro"; mkdir -p "$W/co" "$W/ro"
"$C" decompress --db "$D/db.features.tsv" --out-dir "$W/co" -i "$W/cdb.kst" >/dev/null 2>&1
"$R" decompress --db "$D/db.features.tsv" --out-dir "$W/ro" -i "$W/cdb.kst" >/dev/null 2>&1
if diff -r "$W/co" "$W/ro" >/dev/null 2>&1; then echo "--db decompress: identical"
else echo "--db decompress: DIFFERS"; fail=1; fi
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
[ $dn = 0 ] && echo "--db decompress: restores all six inputs" \
  || { echo "--db decompress: $dn inputs not restored"; fail=1; }

# 10. --aln, both directions, and aln-pack.
aok=0; abad=0
for pair in "s.bed p.paf ref.bed" "ls.bed lp.paf lref.bed"; do
  set -- $pair
  for fl in "" "--block 4" "--effort low"; do
    "$C" $fl compress --aln "$D/$2" --ref "$D/$3" < "$D/$1" > "$W/c.kst" 2>/dev/null; c1=$?
    "$R" $fl compress --aln "$D/$2" --ref "$D/$3" < "$D/$1" > "$W/r.kst" 2>/dev/null; c2=$?
    if cmp -s "$W/c.kst" "$W/r.kst" && [ $c1 = $c2 ]; then aok=$((aok+1))
    else abad=$((abad+1)); echo "  --aln compress differs: $1 [$fl]"; fi
  done
  "$C" compress --aln "$D/$2" --ref "$D/$3" < "$D/$1" > "$W/a.kst" 2>/dev/null
  "$C" decompress --aln "$D/$2" --ref "$D/$3" -i "$W/a.kst" > "$W/c.out" 2>/dev/null; c1=$?
  "$R" decompress --aln "$D/$2" --ref "$D/$3" -i "$W/a.kst" > "$W/r.out" 2>/dev/null; c2=$?
  if cmp -s "$W/c.out" "$W/r.out" && [ $c1 = $c2 ] && cmp -s "$W/r.out" "$D/$1"; then aok=$((aok+1))
  else abad=$((abad+1)); echo "  --aln decompress differs: $1"; fi
  # a .kst reference: same bytes and rows
  "$C" compress < "$D/$3" > "$W/c.ref.kst" 2>/dev/null; "$R" compress < "$D/$3" > "$W/r.ref.kst" 2>/dev/null
  "$C" compress --aln "$D/$2" --ref "$W/r.ref.kst" < "$D/$1" > "$W/c.kst" 2>/dev/null; c1=$?
  "$R" compress --aln "$D/$2" --ref "$W/c.ref.kst" < "$D/$1" > "$W/r.kst" 2>/dev/null; c2=$?
  "$R" decompress --aln "$D/$2" --ref "$W/c.ref.kst" -i "$W/a.kst" > "$W/r.out" 2>/dev/null; c3=$?
  if cmp -s "$W/c.ref.kst" "$W/r.ref.kst" && cmp -s "$W/c.kst" "$W/a.kst" && cmp -s "$W/r.kst" "$W/a.kst" && [ $c1 = 0 ] && [ $c2 = 0 ] && [ $c3 = 0 ] && cmp -s "$W/r.out" "$D/$1"; then aok=$((aok+1))
  else abad=$((abad+1)); echo "  .kst reference differs: $1"; fi
done
# one-block query: empty deltas stream
printf 's\t0\t1\tA\n' > "$W/one.ref.bed"; printf 'q\t1\t0\t1\t+\ts\t1\t0\t1\t1\t1\t60\tcg:Z:1M\n' > "$W/one.paf"
for pair in "p.paf ref.bed" "lp.paf lref.bed" "db.self.paf db.rall.bed" "$W/one.paf $W/one.ref.bed"; do
  set -- $pair; a=$1; b=$2; [ -e "$a" ] || a=$D/$1; [ -e "$b" ] || b=$D/$2
  "$C" aln-pack -a "$a" --ref "$b" -o "$W/c.pk" >/dev/null 2>"$W/c.pe"; c1=$?
  "$R" aln-pack -a "$a" --ref "$b" -o "$W/r.pk" >/dev/null 2>"$W/r.pe"; c2=$?
  if cmp -s "$W/c.pk" "$W/r.pk" && cmp -s "$W/c.pe" "$W/r.pe" && [ $c1 = $c2 ]; then aok=$((aok+1))
  else abad=$((abad+1)); echo "  aln-pack differs: $1"; fi
done
[ $abad = 0 ] && echo "--aln and aln-pack: $aok forms identical" \
  || { echo "--aln: $abad DIFFER"; fail=1; }

# 11. same refusals, same words
EXPECT=fail
for bad_in in 'r1\t0\t10' 'r1\t5\t10\tA' 'r1\t0\t0\tA' '\t0\t10\tA' 'r1\tx\t10\tA' 'r1\t00\t10\tA' 'r1\t0\t010\tA'; do
  printf "$bad_in\n" | "$C" compress 2> "$W/ce" >/dev/null; c1=$?; cm=$(tail -1 "$W/ce")
  printf "$bad_in\n" | "$R" compress 2> "$W/re" >/dev/null; r1=$?; rm2=$(tail -1 "$W/re")
  [ "$cm" = "$rm2" ] && [ "$c1" = 1 ] && [ "$r1" = 1 ] || { echo "  refusal differs on [$bad_in]: c=[$cm] rc $c1, rust=[$rm2] rc $r1"; fail=1; }
done
# refusals on companion files
printf 'q1\n\nq2\n' > "$W/blank"; gzip -c "$D/ref.bed" | head -c 40 > "$W/cut.gz"
printf '\ns\t0\t10\t1\n' > "$W/dbblank.bed"; printf 's\t00\t10\t1\n' > "$W/dbzero.bed"
"$C" aln-pack -a "$D/p.paf" --ref "$D/ref.bed" -o "$W/over.pk" >/dev/null 2>&1; printf '\0' >> "$W/over.pk"
"$C" compress --aln "$D/p.paf" --ref "$D/ref.bed" < "$D/s.bed" > "$W/a.kst" 2>/dev/null
for form in "compress --aln $D/p.paf --ref $D/ref.bed --aln-names $W/blank" "compress --aln $D/p.paf --ref $W/cut.gz" "decompress --aln $W/over.pk --ref $D/ref.bed -i $W/a.kst" \
            "--db $D/db.features.tsv compress -i $W/dbblank.bed" "--db $D/db.features.tsv compress -i $W/dbzero.bed"; do
  "$C" $form < "$D/s.bed" 2> "$W/ce" >/dev/null; c1=$?; cm=$(tail -1 "$W/ce")
  "$R" $form < "$D/s.bed" 2> "$W/re" >/dev/null; r1=$?; rm2=$(tail -1 "$W/re")
  [ "$cm" = "$rm2" ] && [ "$c1" = 1 ] && [ "$r1" = 1 ] || { echo "  refusal differs on [$form]: c=[$cm] rc $c1, rust=[$rm2] rc $r1"; fail=1; }
done
echo "refusals: same message from both"

# 12. awkward valid input: same bytes
python3 "$D/kt.py" gen "$W" 20
gok=0; gbad=0
for f in "$W"/ok*.bed; do
  for fl in "" "--block 4" "--block 262144" "--names drop" "--effort low --block 3"; do
    "$C" $fl compress < "$f" > "$W/c.kst" 2>/dev/null; c1=$?
    "$R" $fl compress < "$f" > "$W/r.kst" 2>/dev/null; c2=$?
    if cmp -s "$W/c.kst" "$W/r.kst" && [ $c1 = $c2 ]; then gok=$((gok+1))
    else gbad=$((gbad+1)); echo "  awkward input differs: $(basename "$f") [$fl]"; fi
  done
done
[ $gbad = 0 ] && echo "awkward input: $gok forms identical" || { echo "awkward input: $gbad DIFFER"; fail=1; }

EXPECT=ok
# 11b. gzip input: same track as plain
gok=0; gbad=0
for f in "$D"/s.bed "$D"/ls.bed "$D"/db.rall.bed "$D"/ord.bed; do
  gzip -c "$f" > "$W/in.gz"
  "$C" compress < "$f" > "$W/c.plain" 2>/dev/null
  for b in "$C" "$R"; do
    "$b" compress < "$W/in.gz" > "$W/g.kst" 2>/dev/null
    if cmp -s "$W/c.plain" "$W/g.kst"; then gok=$((gok+1))
    else gbad=$((gbad+1)); echo "  gz differs: $(basename "$f") $b"; fi
  done
done
gzip -c "$D/ref.bed" > "$W/ref.gz"; gzip -c "$D/p.paf" > "$W/paf.gz"
gzip -c "$D/s.bed" > "$W/s.gz"
"$C" compress --aln "$D/p.paf" --ref "$D/ref.bed" < "$D/s.bed" > "$W/a.plain" 2>/dev/null
for b in "$C" "$R"; do
  "$b" compress --aln "$W/paf.gz" --ref "$W/ref.gz" < "$W/s.gz" > "$W/a.gz.kst" 2>/dev/null
  if cmp -s "$W/a.plain" "$W/a.gz.kst"; then gok=$((gok+1))
  else gbad=$((gbad+1)); echo "  gz --aln differs: $b"; fi
done
# many BGZF blocks, file and pipe
awk 'BEGIN{for(i=0;i<60000;i++){p=0;for(j=0;j<5;j++){printf "read%d\t%d\t%d\tL%d\n",i,p,p+j+7,j%3;p+=j+7}}}' > "$W/bz.bed"
python3 "$D/kt.py" bgzf "$W/bz.bed" "$W/bz.bgz"
"$C" compress < "$W/bz.bed" > "$W/bz.plain" 2>/dev/null
for b in "$C" "$R"; do
  for how in file pipe; do
    if [ $how = file ]; then "$b" compress -i "$W/bz.bgz" > "$W/bz.kst" 2>/dev/null; else "$b" compress < "$W/bz.bgz" > "$W/bz.kst" 2>/dev/null; fi
    if cmp -s "$W/bz.plain" "$W/bz.kst"; then gok=$((gok+1)); else gbad=$((gbad+1)); echo "  bgzf differs: $b $how"; fi
  done
done
sz=$(stat -c%s "$W/bz.bgz"); head -c $((sz / 2)) "$W/bz.bgz" > "$W/bzcut.bgz"   # cut and byte-flipped BGZF
cp "$W/bz.bgz" "$W/bzflip.bgz"; printf 'U' | dd of="$W/bzflip.bgz" bs=1 seek=$((sz / 3)) conv=notrunc status=none
for f in bzcut bzflip; do
  "$C" compress -i "$W/$f.bgz" 2> "$W/ce" >/dev/null; c1=$?; cm=$(tail -1 "$W/ce")
  "$R" compress -i "$W/$f.bgz" 2> "$W/re" >/dev/null; r1=$?; rm2=$(tail -1 "$W/re")
  [ "$cm" = "$rm2" ] && [ "$c1" = 1 ] && [ "$r1" = 1 ] && gok=$((gok+1)) || { gbad=$((gbad+1)); echo "  $f refusal differs: c=[$cm] rc $c1, rust=[$rm2] rc $r1"; }
done
[ $gbad = 0 ] && echo "gzip input: $gok forms match the plain track" \
  || { echo "gzip input: $gbad DIFFER"; fail=1; }

# 12. CLI: same words, same exit code
cok=0; cbad=0
cli() {
  (cd "$W" && "$C" "$@" </dev/null > cliout.c 2>&1); local cr=$?   # in $W, since forms may write files
  (cd "$W" && "$R" "$@" </dev/null > cliout.r 2>&1); local rr=$?
  if cmp -s "$W/cliout.c" "$W/cliout.r" && [ "$cr" = "$rr" ]; then cok=$((cok+1)); return; fi
  cbad=$((cbad+1)); echo "  cli differs: [$*] rc $cr/$rr"; echo "    c=[$(head -c 200 "$W/cliout.c" | tr -d '\0' | head -1)]"
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
cli decompress --max-gap x -i f --hier h --smooth --set fs
cli decompress --out-dir D; cli decompress -i f --halo 2; cli query -i f -r s:+1-5; cli query -i f -r s:99999999999999999999-1
cli compress --block +5; cli compress --aln p; cli compress --ref r
cli decompress --aln-names f; cli index --names ordinal
cli decompress --aln p; cli decompress --ref r; cli decompress --aln p --ref r
cli compress --names drop --aln p --ref r; cli decompress --db f; cli compress --db f
cli decompress -i a -i b; cli compress -i nosuch; cli index -i f -o x; cli query -i f -r s:2-1
cli query -i f -r s --out-dir d; cli query -i f -r s --raw; cli decompress -i f --halo 2; cli region -i f -r s
cli --tmp-dir /nonexistent compress; cli --tmp-dir /nonexistent --db $D/db.features.tsv compress -i $D/db.gene.bed
cli --block=4 compress; cli -b4 compress; cli -t4 compress; cli --effort=low compress; cli -V; cli foo --raw
cli decompress -i f --max-gap -1; cli query -i f -r -x; cli compress -- -x; cli -o -i compress
[ $cbad = 0 ] && echo "cli: $cok forms give the same words and code" \
  || { echo "cli: $cbad of $((cok+cbad)) DIFFER"; fail=1; }

# 13. tail predictor, past TAILSKIP=64 sequences
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
  if cmp -s "$W/c.kst" "$W/r.kst" && "$R" decompress -i "$W/c.kst" 2>/dev/null | cmp -s - "$W/$f.bed"; then tok=$((tok+1)); else tbad=$((tbad+1)); echo "  tail differs: $f $flags"; fi
done; done
[ $tbad = 0 ] && echo "tail predictor: $tok forms identical" || { echo "tail predictor: $tbad DIFFER"; fail=1; }

# 13b. -e high's zpaq frames, both ways
awk 'BEGIN { srand(5); for (c = 0; c < 6; c++) { g = c % 4; nr = 20000 + int(rand() * 20000); p = 0; for (r = 0; r < nr; r++) { l = 1 + int(rand() * (g ? 9 * g : 37)); u = rand(); printf "ctg%d\t%d\t%d\t%s\n", c, p, p + l, u < 0.5 + g / 10 ? "chr" (1 + g) : "novel" int(rand() * 3); p += l } } }' > "$W/zt.bed"
"$C" -t 4 -e high compress -i "$W/zt.bed" > "$W/c.kst" 2>/dev/null; "$R" -t 4 -e high compress -i "$W/zt.bed" > "$W/r.kst" 2>/dev/null
if cmp -s "$W/c.kst" "$W/r.kst" && "$R" decompress -i "$W/c.kst" 2>/dev/null | cmp -s - "$W/zt.bed" && "$C" decompress -i "$W/r.kst" 2>/dev/null | cmp -s - "$W/zt.bed"
then echo "zpaq frames: identical, each reads the other"; else echo "zpaq frames: DIFFER"; fail=1; fi
awk 'BEGIN { srand(7); for (r = 0; r < 60000; r++) { k = rand() < 0.4 ? 1 : 2 + int(rand() * rand() * 12); p = 0; g = int(rand() * 4); for (j = 0; j < k; j++) { l = j == k - 1 ? 120 - p : 1 + int(rand() * (120 - p - (k - j - 1)) / (k - j)); printf "r%d\t%d\t%d\t%s\n", r, p, p + l, "chr" (1 + (g + j) % 5); p += l } } }' > "$W/zi.bed"
"$C" -t 4 -e high compress -i "$W/zi.bed" > "$W/c.kst" 2>/dev/null; "$R" -t 4 -e high compress -i "$W/zi.bed" > "$W/r.kst" 2>/dev/null
if cmp -s "$W/c.kst" "$W/r.kst" && "$R" decompress -i "$W/c.kst" 2>/dev/null | cmp -s - "$W/zi.bed" && "$C" decompress -i "$W/r.kst" 2>/dev/null | cmp -s - "$W/zi.bed"
then echo "interleaved zpaq model: identical, each reads the other"; else echo "interleaved zpaq model: DIFFER"; fail=1; fi

# 14. damaged file: same rows before refusal
python3 -c "
import random; random.seed(9)
with open('$W/many.bed','w') as f:
    for r in range(400000): f.write(f'read{r:07d}\t0\t{random.randint(50,150)}\t{random.choice(\"ABCD\")}\n')"
"$C" -e low --block 4096 compress -i "$W/many.bed" > "$W/many.kst" 2>/dev/null
dok=0; dbad=0
for pct in 75 90; do
  python3 -c "import sys; d=bytearray(open(sys.argv[1],'rb').read()); d[len(d)*$pct//100]^=0xff; open(sys.argv[2],'wb').write(d)" "$W/many.kst" "$W/cut.kst"
  "$C" decompress -i "$W/cut.kst" > "$W/c.out" 2>"$W/c.err"; c1=$?; "$R" decompress -i "$W/cut.kst" > "$W/r.out" 2>"$W/r.err"; r1=$?
  if cmp -s "$W/c.out" "$W/r.out" && cmp -s "$W/c.err" "$W/r.err" && [ $c1 = $r1 ] && [ $c1 != 0 ] && [ -s "$W/c.out" ]; then dok=$((dok+1)); else dbad=$((dbad+1)); echo "  damaged decompress differs at $pct%: $(wc -c < "$W/c.out") vs $(wc -c < "$W/r.out") B, rc $c1/$r1"; fi
done
[ $dbad = 0 ] && echo "damaged decompress: $dok flips stop at the same byte" || { echo "damaged decompress: $dbad DIFFER"; fail=1; }

# 15. forged chunk tables, re-signed
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

[ $fail = 0 ] && echo "DIFF_PASS" || echo "DIFF_FAIL"
exit $fail
