#!/bin/bash
# The real-data gate: every track under CORPUS, default and --fast, through
# both binaries. Fails on a round trip that differs, on C and Rust bytes that
# differ, on a default bigger than --fast, and, given a baseline build, on any
# track the change made bigger. Prints sizes and counts only, never a row.
#   tests/sweep.sh CORPUS [C] [RUST] [BASE_C]
set -u
D=$1; C=${2:-$(dirname "$0")/../c/build/karyozip}; R=${3:-$(dirname "$0")/../rust/target/release/karyozip}; B=${4:-}
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
n=0; bad=0; tn=0; tb=0; td=0
fail() { bad=$((bad+1)); echo "  $*"; }
while read -r f; do
  t=$(basename "$f" .bed.gz)
  for m in default --fast; do mf=""; [ $m = --fast ] && mf=--fast
    zcat "$f" | "$C" -t 8 encode $mf > "$W/c.kst" 2>/dev/null || { fail "encode refused: $t $m"; continue; }
    zcat "$f" | "$R" -t 8 encode $mf > "$W/r.kst" 2>/dev/null || { fail "rust encode refused: $t $m"; continue; }
    cmp -s "$W/c.kst" "$W/r.kst" || fail "C and Rust differ: $t $m"
    zcat "$f" | cmp -s - <("$R" decode "$W/c.kst" 2>/dev/null) || fail "round trip differs: $t $m"
    s=$(stat -c%s "$W/c.kst"); n=$((n+1)); tn=$((tn+s))
    [ $m = default ] && td=$s || { [ $td -gt $s ] && fail "default $td > --fast $s: $t"; }
    if [ -n "$B" ]; then
      zcat "$f" | "$B" -t 8 encode $mf > "$W/b.kst" 2>/dev/null && { o=$(stat -c%s "$W/b.kst"); tb=$((tb+o)); [ $s -gt $o ] && fail "bigger than baseline: $t $m $o -> $s"; }
    fi
  done
done < <(find "$D" -name '*.bed.gz' | sort)
[ -n "$B" ] && echo "sweep: $n encodes, $tn B against baseline $tb B ($(awk -v a=$tb -v b=$tn 'BEGIN{printf "%+.3f%%", (b-a)/a*100}'))" || echo "sweep: $n encodes, $tn B"
[ $bad = 0 ] && echo SWEEP_PASS || { echo "SWEEP_FAIL: $bad"; exit 1; }
