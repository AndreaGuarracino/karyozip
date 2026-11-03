# Shared by run.sh, diff.sh and fuzz.sh: prerequisites, counters, concurrency, timeouts
#   source "$(dirname "$0")/lib.sh"; need TOOL...; ...; summary NAME; exit $fail
# Knobs, all environment variables:
#   CASE=REGEX         only checks whose description matches are reported and counted (see each script)
#   TEST_JOBS=N        the total CPU budget, shared, never multiplied; default SLURM_CPUS_PER_TASK, else nproc.
#                      `make test` runs min(3, N) suites at once and gives each floor(N / suites) (at least 1)
#                      as its own TEST_JOBS; inside a suite it bounds the cases run at once (run.sh's readers,
#                      fuzz.sh's xargs -P) and cargo test's -j and --test-threads
#   KZ_TEST_THREADS=N  the -t every karyozip run gets unless the case tests a thread count; default the
#                      suite's TEST_JOBS. N cases at once share it: each gets KZ_TEST_THREADS / N (share N).
#                      The Rust integration tests pass their share as -t to the CLI they spawn (Makefile)
#   TEST_LIGHT=1       serial and light: skip the stress cases, each with one "skipped (TEST_LIGHT): CASE"
#                      line, counted apart in the summary. Stress cases run the thread counts they test
#                      (-t 2, 4, 8, 1024, the default count, 8 readers at once), an address limit, or size:
#                      every generated input above 20000 rows, -e high on more than the fixtures and the
#                      200-read tail inputs, featureID 2^26-1 (about 1 GB of --db tables for one row);
#                      fuzz.sh drops its forged-index stage (140000 rows) and runs 10 cases per generator,
#                      not 60. Everything else runs; the largest generated input kept is 20000 rows.
#                      Every `stress` call in the suites is one. In every mode two karyozip outputs are
#                      compared one run after the other (seqcmp), not through `cmd | cmp - <(cmd2)`.
#                      TEST_JOBS=1 TEST_LIGHT=1 make test in a 1-CPU Slurm job (two runs, measured): 154 to
#                      170 s, at most 3 runnable threads (a -t 1 compress, its reader and names threads), one
#                      karyozip at a time, at most 302 MB MaxRSS. In an 8-CPU job with KZ_TEST_THREADS=1: 357 s, one
#                      karyozip at a time, at most 4 threads in any (a -t 1 compress; sampled, and strace -f over
#                      diff.sh's cli forms). Run the stress cases in a bounded Slurm job
#   TIMEOUT_SCALE=X    multiplies every fixed timeout; default 1, 4 under a sanitizer binary (see tscale)
#   FUZZ_SEED=N        fuzz.sh's mutation offsets and bytes; printed at the start of every run
# Read by the karyozip binaries themselves (C and Rust, release builds included), set only by run.sh's signal cases:
#   KZ_TEST_SIG_PAUSE_MS=N   a SIGINT/SIGTERM/SIGHUP handler sleeps N ms right after it decides to exit or to defer
#   KZ_TEST_SWAP_PAUSE_MS=N  an output commit (--out-dir's BEDs, -o's rename) sleeps N ms before its window opens and
#                            after each rename, so a case can signal just before and inside the window
#   KZ_TEST_STAGE_PAUSE_MS=N --out-dir sleeps N ms after creating each staged BED, before recording its name
#   Sleeps at those fixed points only: timing, never the bytes. Unset or 0, nothing sleeps. For the test suite only

# need TOOL...: every missing tool named at once, then exit 1
need() {
  local t miss=""
  for t in "$@"; do command -v "$t" >/dev/null 2>&1 || miss="$miss $t"; done
  if [ -n "$miss" ]; then echo "$(basename "$0"): missing tools:$miss (install them, or put them on PATH)"; exit 1; fi
}
# GNU forms the suites use: realpath -m, grep -P, stat -c%s, head -c -N
need_gnu() {
  local bad=""
  realpath -m /nonexistent/x >/dev/null 2>&1 || bad="$bad realpath-m"
  echo a | grep -qP '^a$' 2>/dev/null || bad="$bad grep-P"
  stat -c%s / >/dev/null 2>&1 || bad="$bad stat-c"
  [ "$(printf ab | head -c -1 2>/dev/null)" = a ] || bad="$bad head-c-negative"
  if [ -n "$bad" ]; then echo "$(basename "$0"): needs GNU coreutils and grep, missing:$bad"; exit 1; fi
}

# one variable for every bit of test concurrency
TEST_JOBS=${TEST_JOBS:-${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || echo 1)}}
case $TEST_JOBS in ''|*[!0-9]*|0) TEST_JOBS=1 ;; esac
export TEST_JOBS
KZ_TEST_THREADS=${KZ_TEST_THREADS:-$TEST_JOBS}
case $KZ_TEST_THREADS in ''|*[!0-9]*|0) KZ_TEST_THREADS=1 ;; esac
[ "$KZ_TEST_THREADS" -gt 1024 ] && KZ_TEST_THREADS=1024
export KZ_TEST_THREADS
# share N: the -t for each of N binaries running at once, KZ_TEST_THREADS / N, at least 1
share() { local t=$((KZ_TEST_THREADS / $1)); [ $t -lt 1 ] && t=1; echo $t; }
# pin N: "taskset -c CPUS", the first N CPUs this process may use, for a program with no -t (the C library's
# test, which takes the default count, the allowed CPUs); empty without taskset
pin() {
  command -v taskset >/dev/null 2>&1 || return 0
  echo "taskset -c $(python3 -c 'import os, sys; print(",".join(map(str, sorted(os.sched_getaffinity(0))[:int(sys.argv[1])])))' "$1")"
}

# tscale BINARY: TIMEOUT_SCALE, 4 if unset and BINARY carries ASan; tmo N: N seconds scaled
tscale() {
  if [ -z "${TIMEOUT_SCALE:-}" ]; then
    if [ -n "${1:-}" ] && grep -aq __asan_init "$1" 2>/dev/null; then TIMEOUT_SCALE=4; else TIMEOUT_SCALE=1; fi
  fi
  export TIMEOUT_SCALE
}
tmo() { awk -v a="$1" -v s="${TIMEOUT_SCALE:-1}" 'BEGIN { printf "%g\n", a * s }'; }

# counters: tpass DESC prints DESC, tfail DESC prints "FAIL: DESC"; both count. Neither ever returns non-zero
NPASS=0; NFAIL=0; NSKIP=0; fail=0
want_case() { [ -z "${CASE:-}" ] || [[ $1 =~ $CASE ]]; }
tpass() {
  if [ $# -gt 0 ] && ! want_case "$1"; then return 0; fi
  NPASS=$((NPASS + 1)); if [ $# -gt 0 ]; then echo "$1"; fi; return 0
}
tfail() {
  if ! want_case "$1"; then return 0; fi
  NFAIL=$((NFAIL + 1)); fail=1; echo "FAIL: $1"; return 0
}
# a loop leaf: lb before it, le after it counts one pass if no tfail ran in between
lb() { _LB=$NFAIL; }
le() { [ "$NFAIL" = "$_LB" ] && NPASS=$((NPASS + 1)); return 0; }
# stress DESC: true to run a stress case; under TEST_LIGHT=1 prints "skipped (TEST_LIGHT): DESC", counts one skip, false
stress() {
  [ "${TEST_LIGHT:-0}" = 1 ] || return 0
  if want_case "$1"; then NSKIP=$((NSKIP + 1)); echo "skipped (TEST_LIGHT): $1"; fi; return 1
}
# seqcmp CMD1... -- CMD2...: true if the two commands print the same bytes. They run one after the other, each into
# a file in $W (never at once, as `CMD1 | cmp - <(CMD2)` would run them), their stderr and exit codes ignored
seqcmp() {
  local a=()
  while [ $# -gt 0 ] && [ "$1" != -- ]; do a+=("$1"); shift; done; shift
  "${a[@]}" > "$W/seqcmp.1" 2>/dev/null || :
  "$@" > "$W/seqcmp.2" 2>/dev/null || :
  cmp -s "$W/seqcmp.1" "$W/seqcmp.2"
}
summary() { echo "$1: $NPASS passed, $NFAIL failed, $NSKIP skipped${CASE:+ (CASE=$CASE)}"; }
