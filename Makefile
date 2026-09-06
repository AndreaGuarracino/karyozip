# Two implementations of the same format. c/ is the reference; rust/ must
# produce the same bytes, which tests/diff.sh asserts.
all:
	$(MAKE) -C c

test:
	$(MAKE) -C c test

rust:
	cd rust && cargo build --release

# The same suite the C answers to. It found 13 divergences that diff.sh,
# which only compares outputs, could not see.
rust-test: all rust
	cd rust && cargo test --release
	KT=$(CURDIR)/rust/target/release/karyozip LT=$(CURDIR)/c/build/libtest bash c/tests/run.sh
	bash tests/diff.sh

diff: all rust
	bash tests/diff.sh

# The wide check: every real fixture cell under tests/fixtures, every set,
# smoothing, strategy, flag and reader, through both binaries.
matrix: all rust
	@for d in tests/fixtures/*/; do bash tests/matrix.sh "$$d" || exit 1; done

clean:
	$(MAKE) -C c clean
	cd rust && cargo clean 2>/dev/null || true

.PHONY: all test rust rust-test diff matrix clean
