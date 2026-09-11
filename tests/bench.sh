#!/bin/bash
# The benchmark table: one TSV row per cell, feature set, smoothing, strategy,
# seekability and implementation, on full-size real tracks. Runs on the
# cluster; prints numbers only, never a row or a sequence name.
#   tests/bench.sh CELLDIR PAF FEATURES HIER REFDIR [C] [RUST] [KIND:PART[:SEEK[:IMPL]]] > table.tsv
# CELLDIR holds <sample>.<db>.<set>.<kind>.bed.gz; REFDIR the reference
# tracks named *.<set>.<kind>.bed.gz; PAF may be a packed alignment.
# KIND:PART[:SEEK[:IMPL]] runs one node's share, PART a set or ALL for the
# --db rows; empty runs the lot.
set -u; BAD=0
D=$1; PAF=$2; FEAT=$3; HIER=$4; REFD=$5
C=${6:-$(dirname "$0")/../c/build/karyozip}; R=${7:-$(dirname "$0")/../rust/target/release/karyozip}
PART=${8:-}
T=8; NQ=50; NQ_ALN=50; REPS=1   # an aln-informed query parses the reference track
BLOCKS="65536 1048576"; LEVELS="9 22"; IMPLS="C Rust"   # swept beside the 262144 and 17 defaults
KINDS="presmoothed smoothed"; ONLY=""; SEEKS="seekable non-seekable"
if [ -n "$PART" ]; then IFS=: read -r KINDS ONLY sk im <<< "$PART"
  [ -n "$sk" ] && SEEKS=$sk
  [ -n "$im" ] && case " $IMPLS " in *" $im "*) IMPLS=$im;; *) IMPLS="";; esac; fi
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; export TMPDIR=$W
LASTIN=""; gz=0; plain=0; rows=0; seqs=0
cell=$(basename "$D"); smp=${cell%%.*}; dbv=$(echo "$cell" | cut -d. -f2); backend=$(echo "$cell" | cut -d. -f3 | tr a-z A-Z)
case $smp in asm_*) class=assembly;; ont_*) class=long_ONT;; hifi_*) class=long_HiFi;; short_*) class=short;; *) class=$smp;; esac
ordinal=0; zcat "$(ls "$D"/*.bed.gz | head -1)" | head -1 | cut -f1 | grep -qE '^[0-9]+$' && ordinal=1
alnnames=""; [ $ordinal = 1 ] && [ -s "$D/names.txt" ] && alnnames="--aln-names $D/names.txt"
allsets=$(head -1 "$FEAT" | cut -f2- | tr '\t' ' ')
sets=$(for s in $allsets; do ls "$D"/*."$s".presmoothed.bed.gz >/dev/null 2>&1 && echo $s; done | tr '\n' ' ')   # the ones this cell has

printf 'class\tdb\tbackend\tset\ttype\tinformed\tseekable\tmode\timpl\trows\tseqs\tbed_gz_B\tbed_plain_B\tkst_B\tindex_B\tcompanion_B\ttotal_B\tratio_gz\tratio_plain\tB_per_row\tencode_s\tencode_MB\tdecode_s\tdecode_MB\tsmooth_s\tsmooth_MB\troundtrip\tquery_mean_ms\tquery_sd_ms\tquery_p50_ms\tquery_p99_ms\tthreads\n'

best() { sort -n | head -1; }
refof() {  # refof SET KIND: a cut cell keeps its reference beside its tracks, named ref.*
  { ls "$REFD"/ref.*."$1"."$2".bed.gz; ls "$REFD"/*."$1"."$2".bed.gz; } 2>/dev/null | head -1
}
kstref() {  # kstref BED: the reference as a track, what --ref takes; stored once for the archive, so charged to no row
  zcat "$1" | "$C" -t $T encode > "$W/ref.kst" 2>/dev/null && echo "$W/ref.kst"
}
timed() {  # timed OUTFILE IN cmd... : fastest wall, largest RSS, or "refused refused"
  local out=$1 in=$2; shift 2; : > "$W/t"
  for i in $(seq $REPS); do
    /usr/bin/time -f '%e %M' -a -o "$W/t" "$@" < "$in" > "$out" 2>/dev/null || { echo "refused refused"; return 1; }
  done
  echo "$(cut -d' ' -f1 "$W/t" | best) $(( $(cut -d' ' -f2 "$W/t" | sort -n | tail -1) / 1024 ))"
}
ratio() { awk -v a=$1 -v b=$2 'BEGIN{ if (b > 0) printf "%.2f", a/b; else print "NA" }'; }
queries() {  # queries BIN "FLAGS" KST NAMEFILE : NQ random point queries, ms stats; NA unless every one answered
  local b=$1 fl=$2 k=$3 nf=$4; : > "$W/q"; : > "$W/qbad"
  shuf -r -n $NQ --random-source=<(yes) "$nf" | while read -r nm; do   # with replacement: a cut cell may hold 10 sequences
    s=$(date +%s%N); "$b" $fl region "$k" "$nm" 100 > "$W/qo" 2>/dev/null || echo x >> "$W/qbad"; e=$(date +%s%N)
    [ -s "$W/qo" ] || echo x >> "$W/qbad"
    echo $(( (e - s) / 1000 )) >> "$W/q"
  done
  if [ -s "$W/qbad" ]; then echo "queries: $(wc -l < "$W/qbad") of $NQ refused or empty for $b $fl $k" >&2; echo "NA NA NA NA"; return; fi
  awk '{x[NR]=$1/1000; s+=x[NR]; ss+=x[NR]*x[NR]} END{n=NR; m=s/n; sd=sqrt(ss/n-m*m); asort(x); printf "%.3f %.3f %.3f %.3f", m, sd, x[int(n*0.5)+1], x[int(n*0.99)]}' "$W/q"
}

row() {  # row SET KIND INFORMED SEEK MODE IMPL BIN INPUT(s) ENCFLAGS DECFLAGS COMPANION_B
  local set=$1 kind=$2 inf=$3 seek=$4 mode=$5 impl=$6 bin=$7 in=$8 ef=$9 df=${10} comp=${11}
  local kc="$W/kC.$set.$kind.$inf.$seek.$mode.kst"   # the C file a Rust row must match
  [ "$impl" = C ] && rm -f "$kc"
  # the input's plain text, row and sequence counts and name list: once per
  # input, not per row, since 16 rows share each one and Illumina is 3.3 GB
  if [ "$in" != "$LASTIN" ]; then LASTIN=$in
    gz=0; for f in $in; do gz=$((gz + $(stat -c%s $f))); done
    zcat $in > "$W/plain.bed"; plain=$(stat -c%s "$W/plain.bed"); rows=$(wc -l < "$W/plain.bed")
    cut -f1 "$W/plain.bed" | uniq > "$W/names"; seqs=$(wc -l < "$W/names")
  fi
  # ALL not-informed is the baseline --db is judged against: the same sets as
  # separate plain tracks, their bytes and times summed.
  if [ "$set" = ALL ] && [ "$inf" = not-informed ]; then
    local ts=0 tm=0 tb=0 n=0
    for f in $in; do n=$((n+1))
      set -- $(timed "$W/k.$n.kst" "$f" "$bin" -t $T $ef encode)
      [ "$1" = refused ] && { ts=refused; break; }
      ts=$(awk -v a=$ts -v b=$1 'BEGIN{printf "%.2f", a+b}'); [ $2 -gt $tm ] && tm=$2
      tb=$((tb + $(stat -c%s "$W/k.$n.kst")))
    done
    set -- $ts $tm
  elif [ "$set" = ALL ]; then set -- $(timed "$W/k.kst" /dev/null "$bin" -t $T $ef encode $in)
  else set -- $(timed "$W/k.kst" "$in" "$bin" -t $T $ef encode); fi
  local es=$1 em=$2
  if [ "$es" = refused ]; then
    echo "bench: encode refused on $set $kind $inf $seek $mode $impl" >&2; BAD=$((BAD+1))
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tNA\tNA\t%s\tNA\tNA\tNA\tNA\trefused\tNA\tNA\tNA\tNA\tNA\trefused\tNA\tNA\tNA\tNA\t%s\n' \
      "$class" "$dbv" "$backend" "$set" "$kind" "$inf" "$seek" "$mode" "$impl" "$rows" "$seqs" "$gz" "$plain" "$comp" "$T"
    return
  fi
  local kst=$(stat -c%s "$W/k.kst")
  [ "$comp" = side ] && comp=$(stat -c%s "$W/side")
  local ds dm rt
  if [ "$set" = ALL ] && [ "$inf" = not-informed ]; then
    kst=$tb; ds=0; dm=0; rt=exact; n=0
    for f in $in; do n=$((n+1))
      set -- $(timed "$W/back.bed" /dev/null "$bin" -t $T $df decode "$W/k.$n.kst")
      ds=$(awk -v a=$ds -v b=$1 'BEGIN{printf "%.2f", a+b}'); [ $2 -gt $dm ] && dm=$2
      zcat "$f" | cmp -s - "$W/back.bed" || rt=DIFFERS
    done
    rm -f "$W"/k.*.kst
  elif [ "$set" = ALL ]; then mkdir -p "$W/o"; rm -f "$W/o"/*
    set -- $(timed /dev/null /dev/null "$bin" -t $T $df --out "$W/o" decode "$W/k.kst"); ds=$1; dm=$2
    rt=exact; for f in $in; do s2=$(basename "$f" | sed -E 's/.*\.([a-z_]+)\.(pre)?smoothed\.bed\.gz/\1/'); zcat "$f" | cmp -s - "$W/o/$s2.bed" || rt=DIFFERS; done
  else set -- $(timed "$W/back.bed" /dev/null "$bin" -t $T $df decode "$W/k.kst"); ds=$1; dm=$2
    cmp -s "$W/plain.bed" "$W/back.bed" && rt=exact || rt=DIFFERS
  fi
  # sorting permutes whole records, so both sides are sorted before the
  # comparison, as the matrix does: the verdict stays exact or DIFFERS
  if [ "$inf" = sorted ]; then
    sort "$W/plain.bed" > "$W/w.bed"; sort "$W/back.bed" | cmp -s - "$W/w.bed" && rt=exact || rt=DIFFERS
  fi
  if [ "$rt" != exact ]; then echo "bench: round trip $rt on $set $kind $inf $seek $mode $impl" >&2; BAD=$((BAD+1))
  fi
  # one format, two writers: a Rust row must match the same row's C file
  if [ "$set" != ALL ] || [ "$inf" != not-informed ]; then
    if [ "$impl" = C ]; then cp "$W/k.kst" "$kc"
    elif [ -s "$kc" ] && ! cmp -s "$W/k.kst" "$kc"; then
      echo "bench: C and Rust bytes differ on $set $kind $inf $seek $mode" >&2; BAD=$((BAD+1)); fi
  fi
  # an unsmoothed track is also smoothed, timed: what the tool replaces
  local sm="NA NA"
  if [ "$kind" = presmoothed ] && [ "$set" != ALL ] && [ -s "$HIER" ]; then
    set -- $(timed "$W/sm.bed" /dev/null "$bin" -t $T $df smooth "$W/k.kst" "$HIER" "$set"); sm="$1 $2"
    # and must give KaryoScope's own smoothed track, unless a window clipped the ends
    local smo=$(ls "$D"/*."$set".smoothed.bed.gz 2>/dev/null | grep -v '/ref\.' | head -1)
    if [ "$inf" = not-informed ] && [ -s "$smo" ] && [ ! -e "$D/window" ] && ! zcat "$smo" | cmp -s - "$W/sm.bed"; then
      echo "bench: smooth differs from KaryoScope on $set $seek $mode $impl" >&2; BAD=$((BAD+1)); fi
  fi
  local ix=0 q="NA NA NA NA"
  if [ "$seek" = seekable ] && [ "$set" != ALL ]; then
    cp "$W/k.kst" "$W/ki.kst"
    # a track whose names are in a sidecar indexes from that sidecar, and
    # only that flag: index takes no --aln or --db
    local ixf=""; case "$df" in "--names "*) ixf="$df";; esac
    if "$bin" $ixf index "$W/ki.kst" >/dev/null 2>&1; then ix=$(( $(stat -c%s "$W/ki.kst") - kst ))
      q=$(NQ=$(case "$inf" in not-informed|db-informed|names-dropped|sorted|block-*|level-*) echo $NQ;; *) echo $NQ_ALN;; esac) queries "$bin" "$df" "$W/ki.kst" "$W/names")
    else echo "index refused on $set $kind $inf $mode $impl" >&2; ix=NA; fi
  fi
  local tot=$((kst + comp))
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$class" "$dbv" "$backend" "$set" "$kind" "$inf" "$seek" "$mode" "$impl" "$rows" "$seqs" "$gz" "$plain" "$kst" "$ix" "$comp" "$tot" \
    "$(ratio $gz $tot)" "$(ratio $plain $tot)" "$(awk -v a=$kst -v b=$rows 'BEGIN{ if (b>0) printf "%.3f", a/b; else print "NA" }')" \
    "$es" "$em" "$ds" "$dm" $sm "$rt" $q "$T"
}

pafB=0; [ -s "$PAF" ] && pafB=$(stat -c%s "$PAF")   # no alignment, no aln rows
for kind in $KINDS; do
  for set in $sets; do [ -n "$ONLY" ] && [ "$ONLY" != "$set" ] && continue
    q=$(ls "$D"/*."$set"."$kind".bed.gz 2>/dev/null | grep -v '/ref\.' | head -1); [ -n "$q" ] || { echo "bench: no $set $kind track in the cell" >&2; BAD=$((BAD+1)); continue; }
    ref=""; [ -s "$PAF" ] && { ref=$(refof "$set" "$kind"); [ -s "$ref" ] && ref=$(kstref "$ref"); [ -s "$ref" ] || { echo "bench: no $set $kind reference track" >&2; BAD=$((BAD+1)); }; }
    for seek in $SEEKS; do bf=""; [ $seek = seekable ] && bf="--queryable"
      for mode in default fast; do mf=""; [ $mode = fast ] && mf="--fast"
        for impl in $IMPLS; do bin=$C; [ $impl = Rust ] && bin=$R
          row "$set" "$kind" not-informed "$seek" "$mode" "$impl" "$bin" "$q" "$bf $mf" "" 0
          [ -s "$ref" ] && row "$set" "$kind" aln-informed "$seek" "$mode" "$impl" "$bin" "$q" "$bf $mf --aln $PAF --ref $ref $alnnames" "--aln $PAF --ref $ref $alnnames" "$pafB"
          # what the names cost and what reordering buys, on the
          # presmoothed track the archive stores
          [ "$kind" = presmoothed ] && [ "$seek" = seekable ] && {
            # ordinal-named tracks have no names to move to a sidecar
            [ $ordinal = 0 ] && row "$set" "$kind" names-dropped "$seek" "$mode" "$impl" "$bin" "$q" "$bf $mf --names drop:$W/side" "--names $W/side" side
            [ "$mode" = default ] && row "$set" "$kind" sorted "$seek" "$mode" "$impl" "$bin" "$q" "$bf $mf --sort" "" 0
            # a block is decoded whole, so its size trades bytes against the
            # query; the level trades bytes against the encode
            [ "$mode" = default ] && for v in $BLOCKS; do
              row "$set" "$kind" "block-$v" "$seek" "$mode" "$impl" "$bin" "$q" "--block $v" "" 0; done
            [ "$mode" = default ] && for v in $LEVELS; do
              row "$set" "$kind" "level-$v" "$seek" "$mode" "$impl" "$bin" "$q" "--level $v" "" 0; done; }
        done
      done
    done
  done
  # --db takes every set at once: one row per kind, set = ALL
  [ -n "$ONLY" ] && [ "$ONLY" != ALL ] && continue
  ins=""; for set in $sets; do f=$(ls "$D"/*."$set"."$kind".bed.gz 2>/dev/null | grep -v '/ref\.' | head -1); [ -n "$f" ] && ins="$ins $f"; done
  [ $(echo $ins | wc -w) = $(echo $allsets | wc -w) ] || { echo "bench: $kind --db rows skipped, the cell lacks a set its database lists" >&2; continue; }
  ref=""; [ -s "$PAF" ] && { ref=$(refof chromosome "$kind"); [ -s "$ref" ] && ref=$(kstref "$ref"); [ -s "$ref" ] || { echo "bench: no chromosome $kind reference track" >&2; BAD=$((BAD+1)); }; }
  for seek in $SEEKS; do bf=""; [ $seek = seekable ] && bf="--queryable"
    for mode in default fast; do mf=""; [ $mode = fast ] && mf="--fast"
      for impl in $IMPLS; do bin=$C; [ $impl = Rust ] && bin=$R
        row ALL "$kind" not-informed "$seek" "$mode" "$impl" "$bin" "$ins" "$bf $mf" "" 0
        row ALL "$kind" db-informed "$seek" "$mode" "$impl" "$bin" "$ins" "--db $FEAT --hier $HIER $bf $mf" "--db $FEAT --hier $HIER" 0
        [ -s "$ref" ] && row ALL "$kind" db+aln-informed "$seek" "$mode" "$impl" "$bin" "$ins" "--db $FEAT --hier $HIER $bf $mf --aln $PAF --ref $ref $alnnames" "--db $FEAT --hier $HIER --aln $PAF --ref $ref $alnnames" "$pafB"
      done
    done
  done
done
# every row must round-trip: a refusal or a difference fails the run
[ "${BAD:-0}" = 0 ] || { echo "bench: $BAD rows failed a check" >&2; exit 1; }
