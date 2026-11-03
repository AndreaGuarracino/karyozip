# C reference, Rust twin: same bytes
all:
	$(MAKE) -C c

# Test budget, one knob: TEST_JOBS is the total, shared, never multiplied; default the Slurm allocation,
# else every CPU. The three suites run SUITE_JOBS = min(3, TEST_JOBS) at a time, and each gets SUITE_TJ =
# floor(TEST_JOBS / SUITE_JOBS), at least 1, as its own TEST_JOBS: its cases at once, cargo's -j and
# --test-threads, and (as KZ_TEST_THREADS, see tests/lib.sh) the -t of every karyozip it runs.
# TEST_JOBS=1: the suites one after another, each serial. TEST_LIGHT=1 also skips the stress cases.
TEST_JOBS ?= $(or $(SLURM_CPUS_PER_TASK),$(shell nproc 2>/dev/null || echo 1))
SUITE_JOBS := $(shell n=$(TEST_JOBS); [ "$$n" -ge 1 ] 2>/dev/null || n=1; [ "$$n" -gt 3 ] && n=3; echo $$n)
SUITE_TJ := $(shell n=$(TEST_JOBS); [ "$$n" -ge 1 ] 2>/dev/null || n=1; echo $$((n / $(SUITE_JOBS))))
# rayon N: the Rust's own tests run N at once (--test-threads), so each Writer's rayon pool, and each karyozip
# CLI the tests spawn (as -t, from KZ_TEST_THREADS, see rust/tests/lib.rs), gets KZ_TEST_THREADS (default N) / N
# threads, at least 1
rayon = $(shell t=$${KZ_TEST_THREADS:-$(1)}; t=$$((t / $(1))); [ $$t -ge 1 ] || t=1; echo $$t)
export TEST_JOBS

rust:
	cd rust && cargo build --release -j $(TEST_JOBS)

san:
	$(MAKE) -C c build/karyozip-san

# C suite, Rust suite, C against Rust; every binary built first, so the suites do not compile at once
test: all rust san
	@$(MAKE) --no-print-directory -j$(SUITE_JOBS) -O test-c test-rust test-diff
# env, not make: a TEST_JOBS given on make's command line would reach sub-makes over the suite's share
test-c:
	@env TEST_JOBS=$(SUITE_TJ) bash tests/run.sh c/build/karyozip c/build/libtest c/build/karyozip-san
test-rust:
	@cd rust && env TEST_JOBS=$(SUITE_TJ) RAYON_NUM_THREADS=$(call rayon,$(SUITE_TJ)) KZ_TEST_THREADS=$(call rayon,$(SUITE_TJ)) cargo test --release -q -j $(SUITE_TJ) --test lib -- --test-threads=$(SUITE_TJ)
	@env TEST_JOBS=$(SUITE_TJ) bash tests/run.sh rust/target/release/karyozip c/build/libtest
test-diff:
	@env TEST_JOBS=$(SUITE_TJ) bash tests/diff.sh
# the Rust unit and integration tests in a debug build: overflow checks and debug assertions on
test-rust-debug:
	@cd rust && env RAYON_NUM_THREADS=$(call rayon,$(TEST_JOBS)) KZ_TEST_THREADS=$(call rayon,$(TEST_JOBS)) cargo test -q -j $(TEST_JOBS) -- --test-threads=$(TEST_JOBS)

clean:
	$(MAKE) -C c clean
	cd rust && cargo clean 2>/dev/null || true

.PHONY: all rust san test test-c test-rust test-diff test-rust-debug clean
