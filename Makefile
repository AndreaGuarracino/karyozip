# C reference, Rust twin: same bytes
all:
	$(MAKE) -C c

rust:
	cd rust && cargo build --release

# C suite, Rust suite, C against Rust
test: all rust
	@$(MAKE) --no-print-directory -j3 -O test-c test-rust test-diff
test-c:
	@$(MAKE) --no-print-directory -C c test
test-rust:
	@cd rust && cargo test --release -q --test lib
	@bash tests/run.sh rust/target/release/karyozip c/build/libtest
test-diff:
	@bash tests/diff.sh

clean:
	$(MAKE) -C c clean
	cd rust && cargo clean 2>/dev/null || true

.PHONY: all rust test test-c test-rust test-diff clean
