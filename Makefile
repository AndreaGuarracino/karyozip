# Two implementations of the same format. c/ is the reference; rust/ must
# produce the same bytes, which tests/diff.sh asserts.
all:
	$(MAKE) -C c

rust:
	cd rust && cargo build --release

# Every local check, one target, 6 minutes: the hand-computed cases, leaks
# and fuzz on the C, the same on the Rust (they found 13 divergences diff.sh
# could not see), and C against Rust on generated inputs. The matrix takes 18
# minutes a cell and runs on gemini, tests/matrix.sbatch.
test: all rust
	$(MAKE) -C c test
	KT=$(CURDIR)/rust/target/release/karyozip LT=$(CURDIR)/c/build/libtest bash c/tests/run.sh
	bash tests/diff.sh

# The real-data gate: every track of a corpus, both binaries, both modes;
# BASE=old/karyozip also fails on any track the change made bigger.
CORPUS ?= $(HOME)/karyoscope-realdata/kstracks
sweep: all rust
	bash tests/sweep.sh $(CORPUS) c/build/karyozip rust/target/release/karyozip $(BASE)

clean:
	$(MAKE) -C c clean
	cd rust && cargo clean 2>/dev/null || true

.PHONY: all rust test sweep clean
