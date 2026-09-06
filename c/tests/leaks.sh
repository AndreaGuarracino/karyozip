#!/bin/bash
# LeakSanitizer over every command path. fuzz.sh builds an ASan binary but
# runs it on damaged files, where an early die is expected and a leak on that
# path tells you nothing. This runs the clean paths, where nothing may leak.
#   bash leaks.sh [SANITIZER_BINARY]
set -u
D=$(cd "$(dirname "$0")" && pwd)
SRC=$D/..
SAN=${1:-${KARYOZIP_SAN:-}}; BUILT=""
W=$(mktemp -d)
# a Rust binary under test runs under valgrind, which sees the same leaks
if [ -z "$SAN" ] && [ -n "${KT:-}" ] && grep -q rust_begin_unwind "$KT" 2>/dev/null; then
  command -v valgrind >/dev/null || { echo "leaks: FAILED, no valgrind for the Rust"; exit 1; }
  # glibc's thread-local storage of pool threads still alive at exit is not a leak
  printf '{\n  tls\n  Memcheck:Leak\n  match-leak-kinds: definite\n  fun:calloc\n  fun:allocate_dtv\n}\n' > "$W/supp"
  printf '#!/bin/bash\nulimit -n 4096 2>/dev/null\nexec valgrind -q --leak-check=full --errors-for-leak-kinds=definite --error-exitcode=99 --suppressions="%s" "%s" "$@"\n' "$W/supp" "$KT" > "$W/san"; chmod +x "$W/san"; SAN=$W/san
fi
if [ -z "$SAN" ]; then BUILT=1
  command -v gcc >/dev/null || { echo "leaks: skipped, no gcc"; exit 0; }
  SAN=$(mktemp -u /tmp/kz-lsan.XXXXXX)
  gcc -O1 -g -fsanitize=address -std=gnu11 -fopenmp -D_GNU_SOURCE -DVERSION='"1.5.3"' \
      -o "$SAN" "$SRC"/karyozip.c "$SRC"/htscodecs/tokenise_name3.c \
      "$SRC"/htscodecs/rANS_static4x16pr.c "$SRC"/htscodecs/rANS_static32x16pr.c \
      "$SRC"/htscodecs/arith_dynamic.c "$SRC"/htscodecs/rle.c "$SRC"/htscodecs/pack.c \
      "$SRC"/htscodecs/utils.c "$SRC"/bzip3/libbz3.c -I"$SRC" -lzstd -lz -lm -lpthread \
      >/dev/null 2>&1 || { echo "leaks: sanitizer build FAILED"; exit 1; }
fi
trap 'rm -rf "$W" ${BUILT:+"$SAN"}' EXIT
fail=0
run() {   # description, then the command; a refusal is a failure unless the description says (refused)
  local what=$1; shift
  printf '  %-46s ' "$what"
  "$@" >/dev/null 2>"$W/err"; local rc=$?
  if grep -qE "LeakSanitizer|definitely lost: [1-9]" "$W/err" || [ $rc = 99 ]; then echo "LEAKED"; fail=1
  elif [ $rc = 0 ] && [[ "$what" != *"(refused)"* ]]; then echo clean
  elif [ $rc = 1 ] && [[ "$what" == *"(refused)"* ]]; then echo "clean (refused)"
  else echo "WRONG EXIT $rc"; fail=1; fi
}
printf 'feature_set\tchild\tparent\nfs\tA\tcategorized\nfs\tB\tcategorized\nfs\tC\tcategorized\nfs\tD\tcategorized\n' > "$W/h.tsv"
"$SAN" encode < "$D/ls.bed" > "$W/p.kt" 2>/dev/null
"$SAN" encode --names drop:"$W/side" < "$D/ls.bed" > "$W/s.kt" 2>/dev/null
"$SAN" encode --aln "$D/lp.paf" --ref "$D/lref.bed" < "$D/ls.bed" > "$W/a.kt" 2>/dev/null
DB="$D/db.acrocentric.bed $D/db.chromosome.bed $D/db.gene.bed $D/db.region.bed $D/db.repeat.bed $D/db.subtelomeric.bed"
"$SAN" encode --db "$D/db.features.tsv" $DB > "$W/f.kt" 2>/dev/null
mkdir -p "$W/out"
run "encode"                       "$SAN" encode < "$D/ls.bed"
run "encode --slow"                "$SAN" --slow encode < "$D/ls.bed"
run "encode --block all"           "$SAN" --block all encode < "$D/ls.bed"
run "encode --sort"                "$SAN" --sort encode < "$D/ls.bed"
run "encode --sort=1M"             "$SAN" --tmp-dir "$W" --sort=1M encode < "$D/ls.bed"
run "encode --names drop:F"        "$SAN" encode --names drop:"$W/side2" < "$D/ls.bed"
run "encode --aln --ref"           "$SAN" encode --aln "$D/lp.paf" --ref "$D/lref.bed" < "$D/ls.bed"
run "encode --db"                  "$SAN" encode --db "$D/db.features.tsv" $DB
run "decode"                       "$SAN" decode "$W/p.kt"
run "decode --raw"                 "$SAN" decode --raw "$W/p.kt"
run "decode --names F"             "$SAN" --names "$W/side" decode "$W/s.kt"
run "decode --aln --ref"           "$SAN" decode --aln "$D/lp.paf" --ref "$D/lref.bed" "$W/a.kt"
run "decode --db --out"            "$SAN" decode --db "$D/db.features.tsv" --out "$W/out" "$W/f.kt"
run "region, a name that is there" "$SAN" region "$W/p.kt" r1
run "region, a name that is not (refused)" "$SAN" region "$W/p.kt" nosuchname
run "region SEQ POS"               "$SAN" region "$W/p.kt" r1 150
run "region SEQ START END HALO"    "$SAN" region "$W/p.kt" r1 0 200 3
run "index"                        "$SAN" index "$W/p.kt"
run "region, after index"          "$SAN" region "$W/p.kt" r1 150
run "smooth"                       "$SAN" smooth "$W/p.kt" "$W/h.tsv" fs
run "smooth, one sequence"         "$SAN" smooth "$W/p.kt" "$W/h.tsv" fs r1
run "aln-pack"                     "$SAN" aln-pack "$D/lp.paf" "$D/lref.bed" "$W/pk"
run "--help"                       "$SAN" --help
run "--version"                    "$SAN" --version
[ $fail = 0 ] && echo "leaks: clean" || echo "leaks: FOUND SOMETHING"
exit $fail
