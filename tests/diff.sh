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

# 1. every codec gives the same bytes from both, and round-trips.
for f in db.rall.bed db.chromosome.bed ls.bed; do
  "$R" codecs "$D/$f" 2>/dev/null | grep -q 'roundtrip=false' \
    && { echo "codecs: a round trip failed on $f"; fail=1; }
done
echo "codecs: every candidate round-trips"

# 2. encode must be byte-identical, over every fixture and flag.
ok=0; bad=0
for f in "$D"/*.bed; do
  for fl in "" "--block 4" "--block all" "--level 3" "--level 22" "--names keep" \
            "--slow" "--slow --block 4" "--tail 0" "--tail 1" "--names drop" \
            "--tok-chunk 250000" "--slow --level 22" \
            "--sort" "--sort --names keep" "--sort --names drop" "--sort=1M" \
            "--sort --block 4"; do
    "$C" $fl encode < "$f" > "$W/c.kst" 2>/dev/null
    "$R" $fl encode < "$f" > "$W/r.kst" 2>/dev/null
    if cmp -s "$W/c.kst" "$W/r.kst"; then ok=$((ok+1))
    else bad=$((bad+1)); echo "  encode differs: $(basename "$f") [$fl]"; fi
  done
done
[ $bad = 0 ] && echo "encode: $ok forms byte-identical" \
  || { echo "encode: $bad of $((ok+bad)) DIFFER"; fail=1; }

# 3. decode: the Rust must give back exactly what was encoded.
dok=0; dbad=0
for f in "$D"/*.bed; do
  for fl in "" "--block 4" "--block all" "--level 3" "--names keep"; do
    "$C" $fl encode < "$f" > "$W/t.kst" 2>/dev/null
    "$R" decode "$W/t.kst" 2>/dev/null | cmp -s - "$f" && dok=$((dok+1)) \
      || { dbad=$((dbad+1)); echo "  decode differs: $(basename "$f") [$fl]"; }
  done
done
[ $dbad = 0 ] && echo "decode: $dok forms give the input back" \
  || { echo "decode: $dbad of $((dok+dbad)) DIFFER"; fail=1; }

# 4. the C must read what the Rust wrote.
for f in "$D"/ls.bed "$D"/db.rall.bed; do
  "$R" encode < "$f" > "$W/r.kst" 2>/dev/null
  "$C" decode "$W/r.kst" 2>/dev/null | cmp -s - "$f" \
    || { echo "  the C cannot decode the Rust's $(basename "$f")"; fail=1; }
done
echo "cross-decode: the C reads the Rust's files"

# 5. region: the same rows for every shape of query.
rok=0; rbad=0
for f in ls.bed s.bed db.rall.bed ord.bed; do
  for bl in "" "--block 4"; do
    "$C" $bl encode < "$D/$f" > "$W/t.kst" 2>/dev/null
    for seq in $(awk -F'\t' '{print $1}' "$D/$f" | uniq | head -3); do
      for q in "$seq" "$seq 50" "$seq 0 100" "$seq 0 100 2" "$seq 999999"; do
        "$C" region "$W/t.kst" $q > "$W/c.out" 2>/dev/null; c1=$?
        "$R" region "$W/t.kst" $q > "$W/r.out" 2>/dev/null; c2=$?
        if cmp -s "$W/c.out" "$W/r.out" && [ $c1 = $c2 ]; then rok=$((rok+1))
        else rbad=$((rbad+1)); echo "  region differs: $f [$bl] [$q]"; fi
      done
    done
  done
done
[ $rbad = 0 ] && echo "region: $rok queries identical" \
  || { echo "region: $rbad of $((rok+rbad)) DIFFER"; fail=1; }

# 6. index: the same trailer, and each reads the other's.
"$C" encode < "$D/db.rall.bed" > "$W/ci.kst" 2>/dev/null
cp "$W/ci.kst" "$W/ri.kst"
"$C" index "$W/ci.kst" >/dev/null 2>&1
"$R" index "$W/ri.kst" >/dev/null 2>&1
if cmp -s "$W/ci.kst" "$W/ri.kst"; then echo "index: byte-identical trailer"
else echo "index: TRAILERS DIFFER"; fail=1; fi
"$C" region "$W/ri.kst" seq0005 > "$W/a" 2>/dev/null
"$R" region "$W/ci.kst" seq0005 > "$W/b" 2>/dev/null
cmp -s "$W/a" "$W/b" || { echo "index: cross-read DIFFERS"; fail=1; }

# 6b. a sorted track: the Rust must resolve the ranks the same way.
sdok=0; sdbad=0
for f in ls.bed s.bed db.rall.bed ord.bed; do
  for fl in "--sort" "--sort --names keep"; do
    "$C" $fl encode < "$D/$f" > "$W/t.kst" 2>/dev/null
    "$C" decode "$W/t.kst" > "$W/c.out" 2>/dev/null; c1=$?
    "$R" decode "$W/t.kst" > "$W/r.out" 2>/dev/null; c2=$?
    if cmp -s "$W/c.out" "$W/r.out" && [ $c1 = $c2 ]; then sdok=$((sdok+1))
    else sdbad=$((sdbad+1)); echo "  sorted decode differs: $f [$fl]"; fi
  done
  "$C" --sort --names drop encode < "$D/$f" > "$W/t.kst" 2>/dev/null
  "$C" --names ordinal decode "$W/t.kst" > "$W/c.out" 2>/dev/null; c1=$?
  "$R" --names ordinal decode "$W/t.kst" > "$W/r.out" 2>/dev/null; c2=$?
  if cmp -s "$W/c.out" "$W/r.out" && [ $c1 = $c2 ]; then sdok=$((sdok+1))
  else sdbad=$((sdbad+1)); echo "  sorted+drop differs: $f"; fi
done
[ $sdbad = 0 ] && echo "sorted decode: $sdok forms identical" \
  || { echo "sorted decode: $sdbad DIFFER"; fail=1; }

# 7. the reader's flags: dropped names, a sidecar, --raw.
fok=0; fbad=0
for f in ls.bed s.bed db.rall.bed ord.bed; do
  "$C" encode --names drop < "$D/$f" > "$W/d.kst" 2>/dev/null
  for dfl in "--names ordinal" "--raw"; do
    "$C" $dfl decode "$W/d.kst" > "$W/c.out" 2>/dev/null; c1=$?
    "$R" $dfl decode "$W/d.kst" > "$W/r.out" 2>/dev/null; c2=$?
    if cmp -s "$W/c.out" "$W/r.out" && [ $c1 = $c2 ]; then fok=$((fok+1))
    else fbad=$((fbad+1)); echo "  decode $dfl differs on $f"; fi
  done
  "$C" encode --names "drop:$W/side" < "$D/$f" > "$W/s.kst" 2>/dev/null
  "$C" --names "$W/side" decode "$W/s.kst" > "$W/c.out" 2>/dev/null; c1=$?
  "$R" --names "$W/side" decode "$W/s.kst" > "$W/r.out" 2>/dev/null; c2=$?
  if cmp -s "$W/c.out" "$W/r.out" && [ $c1 = $c2 ]; then fok=$((fok+1))
  else fbad=$((fbad+1)); echo "  sidecar differs on $f"; fi
  "$C" encode < "$D/$f" > "$W/p.kst" 2>/dev/null
  "$C" --raw decode "$W/p.kst" > "$W/c.out" 2>/dev/null
  "$R" --raw decode "$W/p.kst" > "$W/r.out" 2>/dev/null
  cmp -s "$W/c.out" "$W/r.out" && fok=$((fok+1)) \
    || { fbad=$((fbad+1)); echo "  --raw differs on $f"; }
done
[ $fbad = 0 ] && echo "reader flags: $fok forms identical" \
  || { echo "reader flags: $fbad DIFFER"; fail=1; }

# 8. smooth: the same rows, over a real hierarchy shape.
printf 'feature_set\tchild\tparent\nfs\tA\tcategorized\nfs\tB\tcategorized\nfs\tC\tcategorized\nfs\tD\tcategorized\n' > "$W/h.tsv"
sok=0; sbad=0
for f in ls.bed s.bed db.rall.bed; do
  "$C" encode < "$D/$f" > "$W/t.kst" 2>/dev/null
  for q in "" "--v1" "--max-gap 5" "--max-gap 100000" "--halo 1"; do
    "$C" $q smooth "$W/t.kst" "$W/h.tsv" fs > "$W/c.out" 2>/dev/null; c1=$?
    "$R" $q smooth "$W/t.kst" "$W/h.tsv" fs > "$W/r.out" 2>/dev/null; c2=$?
    if cmp -s "$W/c.out" "$W/r.out" && [ $c1 = $c2 ]; then sok=$((sok+1))
    else sbad=$((sbad+1)); echo "  smooth differs: $f [$q]"; fi
  done
  for seq in $(awk -F'\t' '{print $1}' "$D/$f" | uniq | head -2); do
    for q in "$seq" "$seq 0 100"; do
      "$C" smooth "$W/t.kst" "$W/h.tsv" fs $q > "$W/c.out" 2>/dev/null; c1=$?
      "$R" smooth "$W/t.kst" "$W/h.tsv" fs $q > "$W/r.out" 2>/dev/null; c2=$?
      if cmp -s "$W/c.out" "$W/r.out" && [ $c1 = $c2 ]; then sok=$((sok+1))
      else sbad=$((sbad+1)); echo "  smooth differs: $f [$q]"; fi
    done
  done
done
[ $sbad = 0 ] && echo "smooth: $sok forms identical" \
  || { echo "smooth: $sbad DIFFER"; fail=1; }

# 9. --db, both directions.
DB="$D/db.acrocentric.bed $D/db.chromosome.bed $D/db.gene.bed $D/db.region.bed $D/db.repeat.bed $D/db.subtelomeric.bed"
"$C" encode --db "$D/db.features.tsv" $DB > "$W/cdb.kst" 2>/dev/null
"$R" encode --db "$D/db.features.tsv" $DB > "$W/rdb.kst" 2>/dev/null
if cmp -s "$W/cdb.kst" "$W/rdb.kst"; then echo "--db encode: identical"
else echo "--db encode: DIFFERS"; fail=1; fi
rm -rf "$W/co" "$W/ro"; mkdir -p "$W/co" "$W/ro"
"$C" decode --db "$D/db.features.tsv" --out "$W/co" "$W/cdb.kst" >/dev/null 2>&1
"$R" decode --db "$D/db.features.tsv" --out "$W/ro" "$W/cdb.kst" >/dev/null 2>&1
if diff -r "$W/co" "$W/ro" >/dev/null 2>&1; then echo "--db decode: identical"
else echo "--db decode: DIFFERS"; fail=1; fi
dn=0
for b in $DB; do
  n=$(basename "$b" .bed); n=${n#db.}
  cmp -s "$W/ro/$n.bed" "$b" || dn=$((dn+1))
done
[ $dn = 0 ] && echo "--db decode: restores all six inputs" \
  || { echo "--db decode: $dn inputs not restored"; fail=1; }

# 10. --aln, both directions, and aln-pack.
aok=0; abad=0
for pair in "s.bed p.paf ref.bed" "ls.bed lp.paf lref.bed"; do
  set -- $pair
  for fl in "" "--block 4" "--slow"; do
    "$C" $fl encode --aln "$D/$2" --ref "$D/$3" < "$D/$1" > "$W/c.kst" 2>/dev/null; c1=$?
    "$R" $fl encode --aln "$D/$2" --ref "$D/$3" < "$D/$1" > "$W/r.kst" 2>/dev/null; c2=$?
    if cmp -s "$W/c.kst" "$W/r.kst" && [ $c1 = $c2 ]; then aok=$((aok+1))
    else abad=$((abad+1)); echo "  --aln encode differs: $1 [$fl]"; fi
  done
  "$C" encode --aln "$D/$2" --ref "$D/$3" < "$D/$1" > "$W/a.kst" 2>/dev/null
  "$C" decode --aln "$D/$2" --ref "$D/$3" "$W/a.kst" > "$W/c.out" 2>/dev/null; c1=$?
  "$R" decode --aln "$D/$2" --ref "$D/$3" "$W/a.kst" > "$W/r.out" 2>/dev/null; c2=$?
  if cmp -s "$W/c.out" "$W/r.out" && [ $c1 = $c2 ] && cmp -s "$W/r.out" "$D/$1"; then aok=$((aok+1))
  else abad=$((abad+1)); echo "  --aln decode differs: $1"; fi
done
for pair in "p.paf ref.bed" "lp.paf lref.bed" "db.self.paf db.rall.bed"; do
  set -- $pair
  "$C" aln-pack "$D/$1" "$D/$2" "$W/c.pk" >/dev/null 2>&1; c1=$?
  "$R" aln-pack "$D/$1" "$D/$2" "$W/r.pk" >/dev/null 2>&1; c2=$?
  if cmp -s "$W/c.pk" "$W/r.pk" && [ $c1 = $c2 ]; then aok=$((aok+1))
  else abad=$((abad+1)); echo "  aln-pack differs: $1"; fi
done
[ $abad = 0 ] && echo "--aln and aln-pack: $aok forms identical" \
  || { echo "--aln: $abad DIFFER"; fail=1; }

# 11. the same refusals, with the same words.
for bad_in in 'r1\t0\t10' 'r1\t5\t10\tA' 'r1\t0\t0\tA' '\t0\t10\tA' 'r1\tx\t10\tA'; do
  cm=$(printf "$bad_in\n" | "$C" encode 2>&1 >/dev/null | tail -1)
  rm2=$(printf "$bad_in\n" | "$R" encode 2>&1 >/dev/null | tail -1)
  [ "$cm" = "$rm2" ] || { echo "  refusal differs on [$bad_in]: c=[$cm] rust=[$rm2]"; fail=1; }
done
echo "refusals: same message from both"

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
  co=$("$C" "$@" </dev/null 2>&1 | tr -d '\0'); cr=$?
  ro=$("$R" "$@" </dev/null 2>&1 | tr -d '\0'); rr=$?
  if [ "$co" = "$ro" ] && [ "$cr" = "$rr" ]; then cok=$((cok+1)); return; fi
  cbad=$((cbad+1)); echo "  cli differs: [$*]"; echo "    c=[$(echo "$co" | head -1)]"
  echo "    rust=[$(echo "$ro" | head -1)]"
}
cli; cli --help; cli -h; cli --version; cli --bogus; cli --bogus encode
cli decode; cli encode extra; cli --profile; cli --threads; cli --threads 4
cli badcmd; cli region; cli region f; cli index; cli smooth f; cli aln-pack
cli encode --names ordinal; cli decode --tail 1; cli index --slow
cli encode --block; cli encode --level 99; cli encode --block 0
cli -t; cli -t 0; cli -t 99999; cli --tmp-dir; cli decode --raw x --raw
cli index --block 0; cli aln-pack a b; cli smooth f h; cli region f s 1 2 3 4
cli encode --sort=0; cli encode --names foo; cli decode --names keep
cli aln-pack --names foo a b c; cli encode --db; cli --version encode
cli index f x; cli decode f x; cli smooth f h fs a b c d; cli aln-pack a b c d
cli smooth f h fs a b; cli smooth f h fs a b c; cli region f s 1
cli encode --tok-chunk 0; cli encode --id-codec zz; cli smooth --max-gap x f h fs
cli decode --out D; cli region --halo 3 f s; cli encode --tail 2
cli encode --id-codec 0x10; cli encode --id-codec 999999; cli encode --len-codec 0x10004
cli encode --id-codec 0x10 --len-codec 999999; cli encode --aln p; cli encode --ref r
cli decode --aln-names f; cli index --names ordinal; cli encode --sort --sort=1M
cli decode --aln p; cli decode --ref r; cli decode --aln p --ref r
cli encode --names drop --aln p --ref r; cli decode --db f; cli encode --db f
cli encode --names drop --names drop:x; cli encode --names drop:
cli encode --sort --names drop --aln p --ref r --aln-names q
[ $cbad = 0 ] && echo "cli: $cok forms give the same words and code" \
  || { echo "cli: $cbad of $((cok+cbad)) DIFFER"; fail=1; }

[ $fail = 0 ] && echo "DIFF_PASS" || echo "DIFF_FAIL"
exit $fail
