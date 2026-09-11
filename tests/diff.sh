#!/bin/bash
# The equivalence gate: c/ and rust/ must produce the same bytes.
set -u
C=${1:-$(dirname "$0")/../c/build/karyozip}
R=${2:-$(dirname "$0")/../rust/target/release/karyozip}
D=$(cd "$(dirname "$0")/../c/tests" && pwd)
[ -x "$C" ] && [ -x "$R" ] || { echo "build both first: make all rust"; exit 1; }
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
fail=0

cv=$("$C" --version | awk '/zstd|libbz3/{print $1"="$2}' | tr '\n' ' ')
rv=$("$R" --version | awk '/zstd|libbz3/{print $1"="$2}' | tr '\n' ' ')
if [ "$cv" = "$rv" ]; then echo "codec versions: $cv"
else echo "codec versions DIFFER: c=[$cv] rust=[$rv]"; fail=1; fi

# Run both on the same input and arguments: stdout bytes and exit code must
# match. IN names the stdin file. Every section below tallies with these.
ok=0; bad=0; EXPECT=ok
both() {  # both LABEL args...  a shared refusal passes only where EXPECT=fail
  local l=$1; shift
  "$C" "$@" < "${IN:-/dev/null}" > "$W/c.out" 2>/dev/null; local a=$?
  "$R" "$@" < "${IN:-/dev/null}" > "$W/r.out" 2>/dev/null; local b=$?
  if cmp -s "$W/c.out" "$W/r.out" && [ $a = $b ] && { [ $a = 0 ] || [ $EXPECT = fail ]; }; then ok=$((ok+1))
  else bad=$((bad+1)); echo "  differs: $l (rc $a/$b)"; fi
}
tally() { [ $bad = 0 ] && echo "$1: $ok forms identical" || { echo "$1: $bad of $((ok+bad)) DIFFER"; fail=1; }; ok=0; bad=0; }

# 1. every codec gives the same bytes from both, and round-trips.
for f in db.rall.bed db.chromosome.bed ls.bed; do
  "$R" codecs "$D/$f" 2>/dev/null | grep -q 'roundtrip=false' \
    && { echo "codecs: a round trip failed on $f"; fail=1; }
done
echo "codecs: every candidate round-trips"

# 2. encode must be byte-identical, over every fixture and flag.
for f in "$D"/*.bed; do
  for fl in "" "--block 4" "--block 262144" \
            "--effort low" "--effort low --block 4" "--names drop" "--effort high" \
            "--sort" "--sort --names drop" "--sort=1M" \
            "--sort --block 4"; do
    IN=$f both "encode $(basename "$f") [$fl]" $fl encode
  done
done
tally encode

# 3. decode: the Rust must give back exactly what was encoded.
dok=0; dbad=0
for f in "$D"/*.bed; do
  for fl in "" "--block 4" "--block 262144"; do
    "$C" $fl encode < "$f" > "$W/t.kst" 2>/dev/null
    "$R" decode -i "$W/t.kst" 2>/dev/null | cmp -s - "$f" && dok=$((dok+1)) \
      || { dbad=$((dbad+1)); echo "  decode differs: $(basename "$f") [$fl]"; }
  done
done
[ $dbad = 0 ] && echo "decode: $dok forms give the input back" \
  || { echo "decode: $dbad of $((dok+dbad)) DIFFER"; fail=1; }

# 4. the C must read what the Rust wrote.
for f in "$D"/ls.bed "$D"/db.rall.bed; do
  "$R" encode < "$f" > "$W/r.kst" 2>/dev/null
  "$C" decode -i "$W/r.kst" 2>/dev/null | cmp -s - "$f" \
    || { echo "  the C cannot decode the Rust's $(basename "$f")"; fail=1; }
done
echo "cross-decode: the C reads the Rust's files"

# 5. region: the same rows for every shape of query.
for f in ls.bed s.bed db.rall.bed ord.bed; do
  for bl in "" "--block 4"; do
    "$C" $bl encode < "$D/$f" > "$W/t.kst" 2>/dev/null
    for seq in $(awk -F'\t' '{print $1}' "$D/$f" | uniq | head -3); do
      for q in "-r $seq" "-r $seq:50-51" "-r $seq:0-100" "-r $seq:0-100 --halo 2" "-r $seq:999999-1000000"; do
        both "region $f [$bl] [$q]" region $q -i "$W/t.kst"
      done
    done
  done
done
tally region

# 6. index: the same trailer, and each reads the other's.
"$C" encode < "$D/db.rall.bed" > "$W/ci.kst" 2>/dev/null
cp "$W/ci.kst" "$W/ri.kst"
"$C" index -i "$W/ci.kst" >/dev/null 2>&1; irc=$?
"$R" index -i "$W/ri.kst" >/dev/null 2>&1; irr=$?
if [ $irc = 0 ] && [ $irr = 0 ] && cmp -s "$W/ci.kst" "$W/ri.kst"; then echo "index: byte-identical trailer"
else echo "index: TRAILERS DIFFER (rc $irc/$irr)"; fail=1; fi
both "index cross-read" region -i "$W/ci.kst" -r seq0005; tally "index cross-read"

# 6b. a sorted track: the Rust must resolve the ranks the same way.
for f in ls.bed s.bed db.rall.bed ord.bed; do
  for fl in "--sort"; do
    "$C" $fl encode < "$D/$f" > "$W/t.kst" 2>/dev/null
    both "sorted decode $f [$fl]" decode -i "$W/t.kst"
  done
  "$C" --sort --names drop encode < "$D/$f" > "$W/t.kst" 2>/dev/null
  both "sorted+drop $f" --names ordinal decode -i "$W/t.kst"
done
tally "sorted decode"

# 7. the reader's flags: dropped names, a sidecar, --raw.
for f in ls.bed s.bed db.rall.bed ord.bed; do
  "$C" encode --names drop < "$D/$f" > "$W/d.kst" 2>/dev/null
  for dfl in "--names ordinal" "--raw"; do both "decode $dfl on $f" $dfl decode -i "$W/d.kst"; done
  # an ordinal-named input has no names to put beside the track: refused by design
  if "$C" encode --names "drop:$W/side" < "$D/$f" > "$W/s.kst" 2>/dev/null; then both "sidecar $f" --names "$W/side" decode -i "$W/s.kst"
  else EXPECT=fail both "sidecar $f refused" encode --names "drop:$W/side2" < "$D/$f"; EXPECT=ok; fi
  "$C" encode < "$D/$f" > "$W/p.kst" 2>/dev/null
  both "--raw $f" --raw decode -i "$W/p.kst"
done
tally "reader flags"

# 8. smooth: the same rows, over a real hierarchy shape.
printf 'feature_set\tchild\tparent\nfs\tA\tcategorized\nfs\tB\tcategorized\nfs\tC\tcategorized\nfs\tD\tcategorized\n' > "$W/h.tsv"
for f in ls.bed s.bed db.rall.bed; do
  "$C" encode < "$D/$f" > "$W/t.kst" 2>/dev/null
  for q in "" "--v1" "--max-gap 5" "--max-gap 100000"; do both "smooth $f [$q]" $q smooth -i "$W/t.kst" --hier "$W/h.tsv" -f fs; done
  for seq in $(awk -F'\t' '{print $1}' "$D/$f" | uniq | head -2); do
    for q in "-r $seq" "-r $seq:0-100"; do both "smooth $f [$q]" smooth -i "$W/t.kst" --hier "$W/h.tsv" -f fs $q; done
  done
done
tally smooth

# 9. --db, both directions.
DB="-i $D/db.acrocentric.bed -i $D/db.chromosome.bed -i $D/db.gene.bed -i $D/db.region.bed -i $D/db.repeat.bed -i $D/db.subtelomeric.bed"
"$C" encode --db "$D/db.features.tsv" $DB > "$W/cdb.kst" 2>/dev/null
"$R" encode --db "$D/db.features.tsv" $DB > "$W/rdb.kst" 2>/dev/null
if cmp -s "$W/cdb.kst" "$W/rdb.kst"; then echo "--db encode: identical"
else echo "--db encode: DIFFERS"; fail=1; fi
rm -rf "$W/co" "$W/ro"; mkdir -p "$W/co" "$W/ro"
"$C" decode --db "$D/db.features.tsv" --out-dir "$W/co" -i "$W/cdb.kst" >/dev/null 2>&1
"$R" decode --db "$D/db.features.tsv" --out-dir "$W/ro" -i "$W/cdb.kst" >/dev/null 2>&1
if diff -r "$W/co" "$W/ro" >/dev/null 2>&1; then echo "--db decode: identical"
else echo "--db decode: DIFFERS"; fail=1; fi
# a tuple outside features.tsv, given the hierarchy: 2^26 + its rank, on both
{ printf 'feature_set\tchild\tparent\nrepeat\trepeat\tcategorized\nrepeat\tLINE\trepeat\n'; } > "$W/hv.tsv"
awk -F'\t' 'NR==5{$4="repeat"} {print}' OFS='\t' "$D/db.repeat.bed" > "$W/hv.repeat.bed"
HV=${DB/$D\/db.repeat.bed/$W\/hv.repeat.bed}
both "--db --hier" encode --db "$D/db.features.tsv" --hier "$W/hv.tsv" $HV; cp "$W/c.out" "$W/hv.kst"
both "--db --hier decode" --db "$D/db.features.tsv" --hier "$W/hv.tsv" --set repeat decode -i "$W/hv.kst"; tally "--db --hier"
dn=0
for b in $DB; do [ "$b" = -i ] && continue;
  n=$(basename "$b" .bed); n=${n#db.}
  cmp -s "$W/ro/$n.bed" "$b" || dn=$((dn+1))
done
[ $dn = 0 ] && echo "--db decode: restores all six inputs" \
  || { echo "--db decode: $dn inputs not restored"; fail=1; }

# 10. --aln, both directions, and aln-pack.
aok=0; abad=0
for pair in "s.bed p.paf ref.bed" "ls.bed lp.paf lref.bed"; do
  set -- $pair
  for fl in "" "--block 4" "--effort low"; do
    "$C" $fl encode --aln "$D/$2" --ref "$D/$3" < "$D/$1" > "$W/c.kst" 2>/dev/null; c1=$?
    "$R" $fl encode --aln "$D/$2" --ref "$D/$3" < "$D/$1" > "$W/r.kst" 2>/dev/null; c2=$?
    if cmp -s "$W/c.kst" "$W/r.kst" && [ $c1 = $c2 ]; then aok=$((aok+1))
    else abad=$((abad+1)); echo "  --aln encode differs: $1 [$fl]"; fi
  done
  "$C" encode --aln "$D/$2" --ref "$D/$3" < "$D/$1" > "$W/a.kst" 2>/dev/null
  "$C" decode --aln "$D/$2" --ref "$D/$3" -i "$W/a.kst" > "$W/c.out" 2>/dev/null; c1=$?
  "$R" decode --aln "$D/$2" --ref "$D/$3" -i "$W/a.kst" > "$W/r.out" 2>/dev/null; c2=$?
  if cmp -s "$W/c.out" "$W/r.out" && [ $c1 = $c2 ] && cmp -s "$W/r.out" "$D/$1"; then aok=$((aok+1))
  else abad=$((abad+1)); echo "  --aln decode differs: $1"; fi
  # the reference as a track, from either binary: the same bytes and rows as the text gives
  "$C" encode < "$D/$3" > "$W/c.ref.kst" 2>/dev/null; "$R" encode < "$D/$3" > "$W/r.ref.kst" 2>/dev/null
  "$C" encode --aln "$D/$2" --ref "$W/r.ref.kst" < "$D/$1" > "$W/c.kst" 2>/dev/null; c1=$?
  "$R" encode --aln "$D/$2" --ref "$W/c.ref.kst" < "$D/$1" > "$W/r.kst" 2>/dev/null; c2=$?
  "$R" decode --aln "$D/$2" --ref "$W/c.ref.kst" -i "$W/a.kst" > "$W/r.out" 2>/dev/null; c3=$?
  if cmp -s "$W/c.ref.kst" "$W/r.ref.kst" && cmp -s "$W/c.kst" "$W/a.kst" && cmp -s "$W/r.kst" "$W/a.kst" && [ $c1 = 0 ] && [ $c2 = 0 ] && [ $c3 = 0 ] && cmp -s "$W/r.out" "$D/$1"; then aok=$((aok+1))
  else abad=$((abad+1)); echo "  .kst reference differs: $1"; fi
done
for pair in "p.paf ref.bed" "lp.paf lref.bed" "db.self.paf db.rall.bed"; do
  set -- $pair
  "$C" aln-pack -a "$D/$1" --ref "$D/$2" -o "$W/c.pk" >/dev/null 2>&1; c1=$?
  "$R" aln-pack -a "$D/$1" --ref "$D/$2" -o "$W/r.pk" >/dev/null 2>&1; c2=$?
  if cmp -s "$W/c.pk" "$W/r.pk" && [ $c1 = $c2 ]; then aok=$((aok+1))
  else abad=$((abad+1)); echo "  aln-pack differs: $1"; fi
done
[ $abad = 0 ] && echo "--aln and aln-pack: $aok forms identical" \
  || { echo "--aln: $abad DIFFER"; fail=1; }

# 11. the same refusals, with the same words.
EXPECT=fail
for bad_in in 'r1\t0\t10' 'r1\t5\t10\tA' 'r1\t0\t0\tA' '\t0\t10\tA' 'r1\tx\t10\tA'; do
  printf "$bad_in\n" | "$C" encode 2> "$W/ce" >/dev/null; c1=$?; cm=$(tail -1 "$W/ce")
  printf "$bad_in\n" | "$R" encode 2> "$W/re" >/dev/null; r1=$?; rm2=$(tail -1 "$W/re")
  [ "$cm" = "$rm2" ] && [ "$c1" = 1 ] && [ "$r1" = 1 ] || { echo "  refusal differs on [$bad_in]: c=[$cm] rc $c1, rust=[$rm2] rc $r1"; fail=1; }
done
# the same refusals on companions: a blank --aln-names line, a cut .gz reference, a packed alignment with a byte over
printf 'q1\n\nq2\n' > "$W/blank"; gzip -c "$D/ref.bed" | head -c 40 > "$W/cut.gz"
"$C" aln-pack -a "$D/p.paf" --ref "$D/ref.bed" -o "$W/over.pk" >/dev/null 2>&1; printf '\0' >> "$W/over.pk"
"$C" encode --aln "$D/p.paf" --ref "$D/ref.bed" < "$D/s.bed" > "$W/a.kst" 2>/dev/null
for form in "encode --aln $D/p.paf --ref $D/ref.bed --aln-names $W/blank" "encode --aln $D/p.paf --ref $W/cut.gz" "decode --aln $W/over.pk --ref $D/ref.bed -i $W/a.kst"; do
  "$C" $form < "$D/s.bed" 2> "$W/ce" >/dev/null; c1=$?; cm=$(tail -1 "$W/ce")
  "$R" $form < "$D/s.bed" 2> "$W/re" >/dev/null; r1=$?; rm2=$(tail -1 "$W/re")
  [ "$cm" = "$rm2" ] && [ "$c1" = 1 ] && [ "$r1" = 1 ] || { echo "  refusal differs on [$form]: c=[$cm] rc $c1, rust=[$rm2] rc $r1"; fail=1; }
done
echo "refusals: same message from both"

# 12. valid input of every awkward shape, the fuzz's generator: both must
# write the same bytes under every flag form.
python3 "$D/kt.py" gen "$W" 20
gok=0; gbad=0
for f in "$W"/ok*.bed; do
  for fl in "" "--block 4" "--block 262144" "--sort" "--names drop" "--effort low --block 3"; do
    "$C" $fl encode < "$f" > "$W/c.kst" 2>/dev/null; c1=$?
    "$R" $fl encode < "$f" > "$W/r.kst" 2>/dev/null; c2=$?
    if cmp -s "$W/c.kst" "$W/r.kst" && [ $c1 = $c2 ]; then gok=$((gok+1))
    else gbad=$((gbad+1)); echo "  awkward input differs: $(basename "$f") [$fl]"; fi
  done
done
[ $gbad = 0 ] && echo "awkward input: $gok forms identical" || { echo "awkward input: $gbad DIFFER"; fail=1; }

EXPECT=ok
# 11b. gzipped input, which is how the real tracks are stored. zlib and
# zlib-rs must inflate to the same bytes, so the track must be identical to
# the one built from the plain file.
gok=0; gbad=0
for f in "$D"/s.bed "$D"/ls.bed "$D"/db.rall.bed "$D"/ord.bed; do
  gzip -c "$f" > "$W/in.gz"
  "$C" encode < "$f" > "$W/c.plain" 2>/dev/null
  for b in "$C" "$R"; do
    "$b" encode < "$W/in.gz" > "$W/g.kst" 2>/dev/null
    if cmp -s "$W/c.plain" "$W/g.kst"; then gok=$((gok+1))
    else gbad=$((gbad+1)); echo "  gz differs: $(basename "$f") $b"; fi
  done
done
gzip -c "$D/ref.bed" > "$W/ref.gz"; gzip -c "$D/p.paf" > "$W/paf.gz"
gzip -c "$D/s.bed" > "$W/s.gz"
"$C" encode --aln "$D/p.paf" --ref "$D/ref.bed" < "$D/s.bed" > "$W/a.plain" 2>/dev/null
for b in "$C" "$R"; do
  "$b" encode --aln "$W/paf.gz" --ref "$W/ref.gz" < "$W/s.gz" > "$W/a.gz.kst" 2>/dev/null
  if cmp -s "$W/a.plain" "$W/a.gz.kst"; then gok=$((gok+1))
  else gbad=$((gbad+1)); echo "  gz --aln differs: $b"; fi
done
[ $gbad = 0 ] && echo "gzip input: $gok forms match the plain track" \
  || { echo "gzip input: $gbad DIFFER"; fail=1; }

# 12. the command line itself: the same words and the same exit code for
# every malformed form, including --help and --version.
cok=0; cbad=0
cli() {
  "$C" "$@" </dev/null > "$W/cliout.c" 2>&1; local cr=$?
  "$R" "$@" </dev/null > "$W/cliout.r" 2>&1; local rr=$?
  if cmp -s "$W/cliout.c" "$W/cliout.r" && [ "$cr" = "$rr" ]; then cok=$((cok+1)); return; fi
  cbad=$((cbad+1)); echo "  cli differs: [$*] rc $cr/$rr"; echo "    c=[$(head -c 200 "$W/cliout.c" | tr -d '\0' | head -1)]"
  echo "    rust=[$(head -c 200 "$W/cliout.r" | tr -d '\0' | head -1)]"
}
cli; cli --help; cli -h; cli --version; cli --bogus; cli --bogus encode
cli decode; cli encode extra; cli --profile; cli --threads; cli --threads 4
cli badcmd; cli region; cli region -i f; cli index; cli smooth -i f; cli aln-pack
cli encode --names ordinal; cli decode --sort; cli index --effort low
cli encode --block; cli encode --effort max; cli encode --block 0
cli -t; cli -t 0; cli -t 99999; cli --tmp-dir; cli decode --raw x --raw
cli index --block 0; cli aln-pack -a a --ref b; cli smooth -i f --hier h; cli region -i f -r s:5-5
cli encode --sort=0; cli encode --names foo; cli decode --names keep
cli aln-pack --names foo -a a --ref b -o c; cli encode --db; cli --version encode
cli index -i f x; cli decode -i f x; cli smooth -i f --hier h -f fs -r s:9-3; cli aln-pack -a a --ref b -o c d
cli smooth -i f --hier h -f fs -r s:1-2; cli region -i f -r s --halo 3; cli region -i f -r s:1-2
cli smooth --max-gap x -i f --hier h -f fs
cli decode --out-dir D; cli region -i f --halo 2; cli region -i f -r s:+1-5; cli region -i f -r s:99999999999999999999-1
cli encode --block +5; cli encode --aln p; cli encode --ref r
cli decode --aln-names f; cli index --names ordinal; cli encode --sort --sort=1M
cli decode --aln p; cli decode --ref r; cli decode --aln p --ref r
cli encode --names drop --aln p --ref r; cli decode --db f; cli encode --db f
cli encode --names drop --names drop:x; cli encode --names drop:
cli encode --sort --names drop --aln p --ref r --aln-names q
cli decode -i a -i b; cli encode -i nosuch; cli index -i f -o x; cli region -i f -r s:2-1
[ $cbad = 0 ] && echo "cli: $cok forms give the same words and code" \
  || { echo "cli: $cbad of $((cok+cbad)) DIFFER"; fail=1; }

[ $fail = 0 ] && echo "DIFF_PASS" || echo "DIFF_FAIL"
exit $fail
