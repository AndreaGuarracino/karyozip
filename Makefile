# Two implementations of the same format. c/ is the reference; rust/ must
# produce the same bytes, which tests/diff.sh asserts.
all:
	$(MAKE) -C c

rust:
	cd rust && cargo build --release

# Every local check, one target: the hand-computed cases, leaks and fuzz on
# the C, the same on the Rust (they found 13 divergences diff.sh could not
# see), and C against Rust on generated inputs. The three are independent and
# run side by side, each one's output kept together.
test: all rust
	@$(MAKE) --no-print-directory -j3 -O test-c test-rust test-diff
test-c:
	@$(MAKE) --no-print-directory -C c test
test-rust:
	@cd rust && cargo test --release -q --test lib
	@KT=$(CURDIR)/rust/target/release/karyozip LT=$(CURDIR)/c/build/libtest bash c/tests/run.sh
test-diff:
	@bash tests/diff.sh

clean:
	$(MAKE) -C c clean
	cd rust && cargo clean 2>/dev/null || true

.PHONY: all rust test test-c test-rust test-diff clean
