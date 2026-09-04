# Two implementations of the same format. c/ is the reference; rust/ must
# produce the same bytes, which tests/diff.sh asserts.
all:
	$(MAKE) -C c

test:
	$(MAKE) -C c test

rust:
	cd rust && cargo build --release

diff: all rust
	bash tests/diff.sh

clean:
	$(MAKE) -C c clean
	cd rust && cargo clean 2>/dev/null || true

.PHONY: all test rust diff clean
