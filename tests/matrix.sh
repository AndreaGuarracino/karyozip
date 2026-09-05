#!/bin/bash
# The wide check: one fixture cell (a real sample cut small), every feature
# set and smoothing, every strategy, every flag that changes the bytes, and
# every reader on every result, run through both binaries and compared.
# Prints only set, kind, flags and byte counts: never a row, never a name.
#   tests/matrix.sh CELLDIR [C] [RUST]
set -u
D=$1; C=${2:-$(dirname "$0")/../c/build/karyozip}; R=${3:-$(dirname "$0")/../rust/target/release/karyozip}
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
export TMPDIR=$W
ok=0; bad=0; refused=0; oracle_ok=0; oracle_miss=0
fail() { bad=$((bad+1)); echo "  DIFFER $*"; }

# both binaries, same stdin and args: stdout bytes and the last stderr line
# must match, and so must the exit code. A shared refusal counts as agreement.
same() {  # same LABEL [IN=file] -- args...
  local lbl=$1; shift; local in=/dev/null
  [ "$1" != -- ] && { in=$1; shift; }; shift
  "$C" "$@" < "$in" > "$W/c.out" 2> "$W/c.err"; local rc1=$?
  "$R" "$@" < "$in" > "$W/r.out" 2> "$W/r.err"; local rc2=$?
  if [ $rc1 != $rc2 ] || ! cmp -s "$W/c.out" "$W/r.out" \
     || [ "$(tail -1 "$W/c.err")" != "$(tail -1 "$W/r.err")" ]; then
    fail "$lbl: rc $rc1/$rc2, $(stat -c%s "$W/c.out")/$(stat -c%s "$W/r.out") B, c=[$(tail -1 "$W/c.err" | cut -c1-70)] r=[$(tail -1 "$W/r.err" | cut -c1-70)]"
    return 1
  fi
  [ $rc1 = 0 ] && ok=$((ok+1)) || refused=$((refused+1))
  return $rc1
}

# the fixture: tracks, references, PAF, tables, names
feat=$(ls "$D"/features.tsv "$D"/*features.txt 2>/dev/null | head -1)
hier=$(ls "$D"/hierarchy.tsv "$D"/*phylogeny.txt 2>/dev/null | head -1)
paf="$D/aln.paf"; names="$D/names.txt"
sets=$(head -1 "$feat" | cut -f2- | tr '\t' ' ')
ordinal=0; zcat "$(ls "$D"/*.bed.gz | grep -v '/ref\.' | head -1)" | head -1 | cut -f1 | grep -qE '^[0-9]+$' && ordinal=1
alnnames=""; [ $ordinal = 1 ] && alnnames="--aln-names $names"
echo "cell $(basename "$D"): sets [$sets], ordinal-named=$ordinal"

FLAGS=("" "--block all" "--block 4096" "--level 3" "--level 22" "--names keep" \
       "--names drop" "--sort" "--sort=1M" "--slow" "--tail 0" "--tail 1")

readers() {  # readers KST SEQ FS HIER "DECODE FLAGS"
  local k=$1 seq=$2 fs=$3 h=$4 dec=$5
  same "decode" -- $dec decode "$k"
  same "decode --raw" -- --raw decode "$k"
  same "region seq" -- $dec region "$k" "$seq" && cp "$W/c.out" "$W/region.out"
  same "region pos" -- $dec region "$k" "$seq" 100
  same "region range" -- $dec region "$k" "$seq" 0 5000
  same "region halo" -- $dec region "$k" "$seq" 0 5000 3
  cp "$k" "$W/ixc.kst"; cp "$k" "$W/ixr.kst"
  "$C" index "$W/ixc.kst" > /dev/null 2> "$W/c.err"; "$R" index "$W/ixr.kst" > /dev/null 2> "$W/r.err"
  if cmp -s "$W/ixc.kst" "$W/ixr.kst" && [ "$(tail -1 "$W/c.err")" = "$(tail -1 "$W/r.err")" ]; then ok=$((ok+1)); else fail "index: $(stat -c%s "$W/ixc.kst")/$(stat -c%s "$W/ixr.kst") B"; fi
  same "region after index" -- $dec region "$W/ixc.kst" "$seq" 0 5000
  same "smooth" -- $dec smooth "$k" "$h" "$fs" && cp "$W/c.out" "$W/smooth.out"
  same "smooth seq" -- $dec smooth "$k" "$h" "$fs" "$seq"
  same "smooth range" -- $dec smooth "$k" "$h" "$fs" "$seq" 0 5000
}
# the plain track's smooth and region are what every informed form must give
plainref() {  # plainref FS KIND
  cp "$W/smooth.out" "$W/plain.$1.$2.smooth"; cp "$W/region.out" "$W/plain.$1.$2.region"; cp "$W/c.out" "$W/plain.$1.$2.smrange"
}
mergecol() {  # mergecol COL: one label column of a --db region, runs merged
  cut -f1-3,"$1" | awk -F'\t' 'NR>1 && $1==p && $2==e && $4==l {e=$3; next} NR>1 {print p"\t"s"\t"e"\t"l} {p=$1; s=$2; e=$3; l=$4} END{if (NR) print p"\t"s"\t"e"\t"l}'
}

for kind in presmoothed smoothed; do
  for fs in $sets; do
    q=$(ls "$D"/*."$fs"."$kind".bed.gz 2>/dev/null | grep -v '/ref\.' | head -1); [ -n "$q" ] || continue
    ref=$(ls "$D"/ref.*."$fs"."$kind".bed.gz 2>/dev/null | head -1)
    zcat "$q" > "$W/in.bed"; seq=$(head -1 "$W/in.bed" | cut -f1)
    for strat in plain aln; do
      [ $strat = aln ] && { [ -s "$ref" ] && [ -s "$paf" ] || continue; }
      a=""; [ $strat = aln ] && a="--aln $paf --ref $ref $alnnames"
      for fl in "${FLAGS[@]}"; do
        lbl="$fs $kind $strat [$fl]"
        same "$lbl encode" "$q" -- $a $fl encode || continue
        cp "$W/c.out" "$W/t.kst"
        # the round trip, with the decoder given what the encoder was given
        dec=""; [ $strat = aln ] && dec="--aln $paf --ref $ref $alnnames"
        case "$fl" in *"--names drop"*) dec="$dec --names ordinal";; esac   # harmless when names were kept
        case "$fl" in *--sort*) sort "$W/in.bed" > "$W/want.bed"; "$C" $dec decode "$W/t.kst" 2>/dev/null | sort > "$W/got.bed";;
                      *) cp "$W/in.bed" "$W/want.bed"; "$C" $dec decode "$W/t.kst" > "$W/got.bed" 2>/dev/null;; esac
        case "$fl" in *"--names drop"*) cut -f2- "$W/want.bed" > "$W/w2"; cut -f2- "$W/got.bed" > "$W/g2"; mv "$W/w2" "$W/want.bed"; mv "$W/g2" "$W/got.bed";; esac
        cmp -s "$W/want.bed" "$W/got.bed" || fail "$lbl round trip"
        # every reader, on the default form only: the bytes are the same file
        if [ -z "$fl" ]; then
          readers "$W/t.kst" "$seq" "$fs" "$hier" "$dec"
          if [ $strat = plain ]; then plainref "$fs" "$kind"
          else cmp -s "$W/smooth.out" "$W/plain.$fs.$kind.smooth" || fail "$lbl smooth differs from the plain track's"
               cmp -s "$W/region.out" "$W/plain.$fs.$kind.region" || fail "$lbl region differs from the plain track's"
               cmp -s "$W/c.out" "$W/plain.$fs.$kind.smrange" || fail "$lbl smooth range differs from the plain track's"; fi
        fi
        # the oracle: smoothing the presmoothed track must give KaryoScope's
        # own smoothed track, from both. A shared miss is the smoother's
        # difference, not ours, and is reported as such.
        if [ -z "$fl" ] && [ $strat = plain ] && [ $kind = presmoothed ]; then
          smo=$(ls "$D"/*."$fs".smoothed.bed.gz 2>/dev/null | grep -v '/ref\.' | head -1)
          if [ -n "$smo" ]; then
            "$C" smooth "$W/t.kst" "$hier" "$fs" > "$W/sc.bed" 2>/dev/null
            "$R" smooth "$W/t.kst" "$hier" "$fs" > "$W/sr.bed" 2>/dev/null
            zcat "$smo" > "$W/ks.bed"
            if cmp -s "$W/sc.bed" "$W/ks.bed" && cmp -s "$W/sr.bed" "$W/ks.bed"; then oracle_ok=$((oracle_ok+1))
            elif cmp -s "$W/sc.bed" "$W/sr.bed"; then oracle_miss=$((oracle_miss+1)); echo "  oracle: $fs smooth differs from KaryoScope, identically in both"
            else fail "$fs smooth vs KaryoScope: C $(cmp -s "$W/sc.bed" "$W/ks.bed" && echo matches || echo differs), Rust $(cmp -s "$W/sr.bed" "$W/ks.bed" && echo matches || echo differs)"; fi
          fi
        fi
      done
      # the other input forms, default flags
      same "$fs $kind $strat gz stdin" "$q" -- $a encode
      same "$fs $kind $strat plain stdin" "$W/in.bed" -- $a encode
      same "$fs $kind $strat -t 1" "$q" -- $a -t 1 encode
    done
  done
  # --db: one BED per set, in header column order, and with --aln on top
  ins=""; for fs in $sets; do f=$(ls "$D"/*."$fs"."$kind".bed.gz 2>/dev/null | grep -v '/ref\.' | head -1); [ -n "$f" ] && ins="$ins $f"; done
  [ $(echo $ins | wc -w) = $(echo $sets | wc -w) ] || { echo "  (db: not every set present for $kind, skipped)"; continue; }
  for strat in db dbaln; do
    a="--db $feat"; [ $strat = dbaln ] && { ref=$(ls "$D"/ref.*.chromosome."$kind".bed.gz | head -1); a="$a --aln $paf --ref $ref $alnnames"; }
    for fl in "" "--block all" "--block 4096" "--sort" "--slow" "--names drop"; do
      lbl="db $kind $strat [$fl]"
      same "$lbl encode" -- $a $fl encode $ins || continue
      cp "$W/c.out" "$W/t.kst"
      mkdir -p "$W/oc" "$W/or"; rm -f "$W/oc"/* "$W/or"/*
      dec="--db $feat --out"; [ $strat = dbaln ] && dec="--aln $paf --ref $ref $alnnames $dec"
      case "$fl" in *"--names drop"*) dec="--names ordinal $dec";; esac
      "$C" $dec "$W/oc" decode "$W/t.kst" >/dev/null 2>&1; "$R" $dec "$W/or" decode "$W/t.kst" >/dev/null 2>&1
      for f in "$W/oc"/*.bed; do cmp -s "$f" "$W/or/$(basename "$f")" || fail "$lbl decode --db $(basename "$f")"; done
      [ -z "$fl" ] && for fs in $sets; do
        zcat "$(ls "$D"/*."$fs"."$kind".bed.gz | grep -v '/ref\.' | head -1)" | cmp -s - "$W/oc/$fs.bed" || fail "$lbl restores $fs"
      done
      # --db readers: the smooth of every set and one label column of region,
      # merged, must equal the plain track's
      [ -z "$fl" ] && { seq=$(zcat $(echo $ins | cut -d' ' -f1) | head -1 | cut -f1); rd=${dec% --out}
        same "$lbl region" -- $rd region "$W/t.kst" "$seq" && { j=0; for fs in $sets; do j=$((j+1))
          mergecol $((j+3)) < "$W/c.out" | cmp -s - <(mergecol 4 < "$W/plain.$fs.$kind.region") || fail "$lbl region column $fs differs from the plain track's"; done; }
        for fs in $sets; do same "$lbl smooth $fs" -- $rd smooth "$W/t.kst" "$hier" "$fs" && { cmp -s "$W/c.out" "$W/plain.$fs.$kind.smooth" || fail "$lbl smooth $fs differs from the plain track's"; }; done; }
    done
  done
done
echo "matrix: $ok identical, $refused refused alike, $bad differ; smooth oracle: $oracle_ok match KaryoScope, $oracle_miss differ in both"
[ $bad = 0 ] && echo MATRIX_PASS || echo MATRIX_FAIL
