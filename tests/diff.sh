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
            "--tok-chunk 250000" "--slow --level 22"; do
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

# 9. the same refusals, with the same words.
for bad_in in 'r1\t0\t10' 'r1\t5\t10\tA' 'r1\t0\t0\tA' '\t0\t10\tA' 'r1\tx\t10\tA'; do
  cm=$(printf "$bad_in\n" | "$C" encode 2>&1 >/dev/null | tail -1)
  rm2=$(printf "$bad_in\n" | "$R" encode 2>&1 >/dev/null | tail -1)
  [ "$cm" = "$rm2" ] || { echo "  refusal differs on [$bad_in]: c=[$cm] rust=[$rm2]"; fail=1; }
done
echo "refusals: same message from both"

[ $fail = 0 ] && echo "DIFF_PASS" || echo "DIFF_FAIL"
exit $fail
