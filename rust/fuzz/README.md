# Fuzz targets

Coverage-guided targets for [cargo-fuzz](https://github.com/rust-fuzz/cargo-fuzz). They are not part of `make test`.

```
cargo install cargo-fuzz
cd rust && cargo +nightly fuzz run reader_open    # a .kst from bytes: Reader::open_with, then decode
cd rust && cargo +nightly fuzz run unzip          # karyozip::unzip: plain, BGZF or gzip text input
```

A refusal (`karyozip::bail`) counts as a pass; any other panic, a sanitizer report or a hang is a finding. Seed `reader_open` with `.kst` files made by `karyozip compress` from `tests/fixtures.txt`, and `unzip` with gzip and BGZF files (`python3 tests/kt.py bgzf IN OUT`).
