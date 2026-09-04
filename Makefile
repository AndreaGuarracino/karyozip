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
	KT=$(CURDIR)/rust/target/release/karyozip bash c/tests/run.sh
	KT=$(CURDIR)/rust/target/release/karyozip bash c/tests/fuzz.sh \
	  $(CURDIR)/rust/target/release/karyozip

diff: all rust
	bash tests/diff.sh

clean:
	$(MAKE) -C c clean
	cd rust && cargo clean 2>/dev/null || true

.PHONY: all test rust rust-test diff clean
