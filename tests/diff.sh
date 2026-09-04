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
  for fl in "" "--block 4" "--block all" "--level 3" "--level 22" "--names keep"; do
    "$C" $fl encode < "$f" > "$W/c.kst" 2>/dev/null
    "$R" $fl encode < "$f" > "$W/r.kst" 2>/dev/null
    if cmp -s "$W/c.kst" "$W/r.kst"; then ok=$((ok+1))
    else bad=$((bad+1)); echo "  encode differs: $(basename "$f") [$fl]"; fi
  done
done
[ $bad = 0 ] && echo "encode: $ok forms byte-identical" \
  || { echo "encode: $bad of $((ok+bad)) DIFFER"; fail=1; }

# 3. the C must read what the Rust wrote.
for f in "$D"/ls.bed "$D"/db.rall.bed; do
  "$R" encode < "$f" > "$W/r.kst" 2>/dev/null
  "$C" decode "$W/r.kst" 2>/dev/null | cmp -s - "$f" \
    || { echo "  the C cannot decode the Rust's $(basename "$f")"; fail=1; }
done
echo "cross-decode: the C reads the Rust's files"

# 4. the same refusals, with the same words.
for bad_in in 'r1\t0\t10' 'r1\t5\t10\tA' 'r1\t0\t0\tA' '\t0\t10\tA' 'r1\tx\t10\tA'; do
  cm=$(printf "$bad_in\n" | "$C" encode 2>&1 >/dev/null | tail -1)
  rm2=$(printf "$bad_in\n" | "$R" encode 2>&1 >/dev/null | tail -1)
  [ "$cm" = "$rm2" ] || { echo "  refusal differs on [$bad_in]: c=[$cm] rust=[$rm2]"; fail=1; }
done
echo "refusals: same message from both"

[ $fail = 0 ] && echo "DIFF_PASS" || echo "DIFF_FAIL"
exit $fail
