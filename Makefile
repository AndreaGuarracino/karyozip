# C reference, Rust twin: same bytes
all:
	$(MAKE) -C c

rust:
	cd rust && cargo build --release

# Test concurrency, one knob: the Slurm allocation, else every CPU. The three suites run up to
# TEST_JOBS at a time (at most 3), and each suite reads TEST_JOBS for its own parallel work.
TEST_JOBS ?= $(or $(SLURM_CPUS_PER_TASK),$(shell nproc 2>/dev/null || echo 1))
SUITE_JOBS := $(shell n=$(TEST_JOBS); [ "$$n" -gt 3 ] 2>/dev/null && n=3; echo $$n)
export TEST_JOBS

# C suite, Rust suite, C against Rust
test: all rust
	@$(MAKE) --no-print-directory -j$(SUITE_JOBS) -O test-c test-rust test-diff
test-c:
	@$(MAKE) --no-print-directory -C c test
test-rust:
	@cd rust && cargo test --release -q --test lib
	@bash tests/run.sh rust/target/release/karyozip c/build/libtest
test-diff:
	@bash tests/diff.sh
# the Rust unit and integration tests in a debug build: overflow checks and debug assertions on
test-rust-debug:
	@cd rust && cargo test -q

clean:
	$(MAKE) -C c clean
	cd rust && cargo clean 2>/dev/null || true

.PHONY: all rust test test-c test-rust test-diff test-rust-debug clean
