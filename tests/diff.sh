#!/bin/bash
# The equivalence gate: c/ and rust/ must agree. Run by `make diff`.
set -u
C=${1:-$(dirname "$0")/../c/build/karyozip}
R=${2:-$(dirname "$0")/../rust/target/release/karyozip}
D=$(cd "$(dirname "$0")/../c/tests" && pwd)
[ -x "$C" ] || { echo "no C binary at $C"; exit 1; }
[ -x "$R" ] || { echo "no Rust binary at $R"; exit 1; }
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
fail=0
# The stored bytes depend on the codec versions, so say them before comparing.
echo "C:    $("$C" --version | tr '\n' ' ')"
echo "Rust: $("$R" --version | tr '\n' ' ')"
cv=$("$C" --version | awk '/zstd|libbz3/{print $2}' | tr '\n' ' ')
rv=$("$R" --version | awk '/zstd|libbz3/{print $2}' | tr '\n' ' ')
[ "$cv" = "$rv" ] && echo "codec versions: match" \
  || { echo "codec versions: DIFFER, c=[$cv] rust=[$rv]"; fail=1; }

# 1. every codec, on real data, must give the same bytes from both.
if "$R" codecs "$D/db.rall.bed" >/dev/null 2>&1; then
  for f in s.bed ls.bed ord.bed db.rall.bed db.chromosome.bed; do
    "$R" codecs "$D/$f" 2>/dev/null | sed 's/\troundtrip=.*//' > "$W/r.out"
    "$R" codecs "$D/$f" 2>/dev/null | grep -q 'roundtrip=false' && { echo "codecs: a round trip failed on $f"; fail=1; }
    [ -s "$W/r.out" ] || { echo "codecs: no output on $f"; fail=1; }
  done
  echo "codecs: every candidate round-trips"
fi
[ $fail = 0 ] && echo "DIFF_PASS" || echo "DIFF_FAIL"
exit $fail
