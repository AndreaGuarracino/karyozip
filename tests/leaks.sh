#!/bin/bash
# LeakSanitizer over every command path. fuzz.sh builds an ASan binary but
# runs it on damaged files, where an early die is expected and a leak on that
# path tells you nothing. This runs the clean paths, where nothing may leak.
#   bash leaks.sh [SANITIZER_BINARY]
set -u
D=$(cd "$(dirname "$0")" && pwd)
SAN=${1:-${KARYOZIP_SAN:-}}
W=$(mktemp -d)
# a Rust binary under test runs under valgrind, which sees the same leaks
if ! [ -x "${SAN%% *}" ] && [ -n "${KT:-}" ] && grep -q rust_begin_unwind "$KT" 2>/dev/null; then
  command -v valgrind >/dev/null || { echo "leaks: FAILED, no valgrind for the Rust"; exit 1; }
  # glibc's thread-local storage of pool threads still alive at exit is not a leak
  printf '{\n  tls\n  Memcheck:Leak\n  match-leak-kinds: definite\n  fun:calloc\n  fun:allocate_dtv\n}\n' > "$W/supp"
  printf '#!/bin/bash\nulimit -n 4096 2>/dev/null\nexec valgrind -q --leak-check=full --errors-for-leak-kinds=definite --error-exitcode=99 --suppressions="%s" "%s" "$@"\n' "$W/supp" "$KT" > "$W/san"; chmod +x "$W/san"; SAN=$W/san
fi
[ -x "${SAN%% *}" ] || { echo "leaks: FAILED, no sanitizer build: make -C c build/karyozip-san"; exit 1; }
trap 'rm -rf "$W"' EXIT
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
DB="-i $D/db.acrocentric.bed -i $D/db.chromosome.bed -i $D/db.gene.bed -i $D/db.region.bed -i $D/db.repeat.bed -i $D/db.subtelomeric.bed"
"$SAN" encode --db "$D/db.features.tsv" $DB > "$W/f.kt" 2>/dev/null
mkdir -p "$W/out"
run "encode"                       "$SAN" encode < "$D/ls.bed"
run "encode --effort low"                "$SAN" --effort low encode < "$D/ls.bed"
run "encode --block 262144"       "$SAN" --block 262144 encode < "$D/ls.bed"
run "encode --sort"                "$SAN" --sort encode < "$D/ls.bed"
run "encode --sort=1M"             "$SAN" --tmp-dir "$W" --sort=1M encode < "$D/ls.bed"
run "encode --names drop:F"        "$SAN" encode --names drop:"$W/side2" < "$D/ls.bed"
run "encode --aln --ref"           "$SAN" encode --aln "$D/lp.paf" --ref "$D/lref.bed" < "$D/ls.bed"
run "encode --db"                  "$SAN" encode --db "$D/db.features.tsv" $DB
run "decode"                       "$SAN" decode -i "$W/p.kt"
run "decode --raw"                 "$SAN" decode --raw -i "$W/p.kt"
run "decode --names F"             "$SAN" --names "$W/side" decode -i "$W/s.kt"
run "decode --aln --ref"           "$SAN" decode --aln "$D/lp.paf" --ref "$D/lref.bed" -i "$W/a.kt"
run "decode --db --out"            "$SAN" decode --db "$D/db.features.tsv" --out-dir "$W/out" -i "$W/f.kt"
run "region, a name that is there" "$SAN" region -i "$W/p.kt" -r r1
run "region, a name that is not (refused)" "$SAN" region -i "$W/p.kt" -r nosuchname
run "region SEQ POS"               "$SAN" region -i "$W/p.kt" -r r1:150-151
run "region SEQ START END HALO"    "$SAN" region -i "$W/p.kt" -r r1:0-200 --halo 3
run "index"                        "$SAN" index -i "$W/p.kt"
run "region, after index"          "$SAN" region -i "$W/p.kt" -r r1:150-151
run "smooth"                       "$SAN" smooth -i "$W/p.kt" --hier "$W/h.tsv" -f fs
run "smooth, one sequence"         "$SAN" smooth -i "$W/p.kt" --hier "$W/h.tsv" -f fs -r r1
run "aln-pack"                     "$SAN" aln-pack -a "$D/lp.paf" --ref "$D/lref.bed" -o "$W/pk"
run "--help"                       "$SAN" --help
run "--version"                    "$SAN" --version
[ $fail = 0 ] && echo "leaks: clean" || echo "leaks: FOUND SOMETHING"
exit $fail
