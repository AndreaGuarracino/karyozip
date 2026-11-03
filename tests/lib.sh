# Shared by run.sh, diff.sh and fuzz.sh: prerequisites, counters, concurrency, timeouts
#   source "$(dirname "$0")/lib.sh"; need TOOL...; ...; summary NAME; exit $fail
# Knobs, all environment variables:
#   CASE=REGEX         only checks whose description matches are reported and counted (see each script)
#   TEST_JOBS=N        parallel test processes; default SLURM_CPUS_PER_TASK, else nproc
#   TIMEOUT_SCALE=X    multiplies every fixed timeout; default 1, 4 under a sanitizer binary (see tscale)
#   FUZZ_SEED=N        fuzz.sh's mutation offsets and bytes; printed at the start of every run

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

# tscale BINARY: TIMEOUT_SCALE, 4 if unset and BINARY carries ASan; tmo N: N seconds scaled
tscale() {
  if [ -z "${TIMEOUT_SCALE:-}" ]; then
    if [ -n "${1:-}" ] && grep -aq __asan_init "$1" 2>/dev/null; then TIMEOUT_SCALE=4; else TIMEOUT_SCALE=1; fi
  fi
  export TIMEOUT_SCALE
}
tmo() { awk -v a="$1" -v s="${TIMEOUT_SCALE:-1}" 'BEGIN { printf "%g\n", a * s }'; }

# counters: tpass DESC prints DESC, tfail DESC prints "FAIL: DESC"; both count. Neither ever returns non-zero
NPASS=0; NFAIL=0; fail=0
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
summary() { echo "$1: $NPASS passed, $NFAIL failed${CASE:+ (CASE=$CASE)}"; }
