# The KST00001 file format

This document specifies the `.kst` track that karyozip 0.2.0 writes, plus its two companion files: the `.kzi` name index and the `KSTAL1` packed alignment.
It was derived from the code, mainly `rust/src/writer.rs` (`Writer::close`), `rust/src/reader.rs` (`Reader::open_with`, `block_dec`, `rmheads`, `ilv_split`, `grp_unpack`, the name decoders, `build_index`, `idxlookup`), `rust/src/kst.rs`, `rust/src/aln.rs` and `c/karyozip.c`.
Where the code and this document disagree, the code is what existing files contain. Please report the disagreement.
Statements marked **(uncertain)** could not be settled from the code with certainty.

Read the [Compatibility](#compatibility) section before you rely on the magic as a version number.

## Conventions

- **u64** is an unsigned 64-bit integer, little-endian. The Rust writes every fixed-width integer with `to_le_bytes`. The C `fwrite`s some of them in host byte order (the header, the tokeniser chunk header, the index trailer and the `KSTAL1` reference hash). So the C writes conforming files only on little-endian hosts. Every supported platform is little-endian.
- **varint** is unsigned LEB128. Each byte carries 7 bits, least significant group first, and the high bit (0x80) means another byte follows. Writers emit the shortest encoding. Readers accept up to 10 bytes, the tenth only as 0 or 1 with no continuation bit (it holds bit 63 alone), and refuse anything past 64 bits or a varint that runs off the end of its buffer with `corrupt varint`. A longer encoding than the shortest is read when it fits.
- **zigzag** maps a signed value to an unsigned one: `zz(v) = (v << 1) ^ (v >> 63)` (arithmetic shift), and the inverse is `unzz(u) = (u >> 1) ^ -(u & 1)`.
- **FNV-1a 64**: start at `h = 0xcbf29ce484222325`. For each byte `c`, do `h = (h ^ c) * 0x100000001b3 mod 2^64`. `fnv(h, bytes)` continues a hash, so `fnv(fnv(h0, a), b) == fnv(h0, a ‖ b)`.
- **zstd frame** means one zstd frame from zstd 1.5.7 (the C links it, and the Rust pins `zstd-sys =2.1.2`, which bundles it). Each frame records its content size: the metadata streams are decoded by that size. The level is not stored and is not needed to decode. The writer uses the compress level (default 17) for the labels and names, and `min(level, 19)` for the index chunks, the index table and the block table.
- A **row** is one BED line `name, start, end, label`. A **sequence** (a read or a contig) is a maximal run of rows with the same name. The **run length** of a row is `end - start`.

## The data model

A track is an ordered list of sequences. Each sequence has a unique name and one or more rows. A sequence's rows start at 0, are gapless (each row starts where the one before ended) and have positive length.
Names are 1 to 511 bytes, with no tab, newline or NUL. Labels contain no newline or NUL.
Decoding reproduces the input BED byte for byte as `name\tstart\tend\tlabel\n`, with decimal coordinates.

## File layout

```
offset 0     magic             8 bytes   "KST00001"
offset 8     header            16 x u64  (128 bytes)
offset 136   labels            h[6] bytes
             names             h[1] & NAMEMASK bytes
             index chunks      h[2] bytes
             block table       h[3] bytes
             index table       h[11] bytes
data0        block bodies      the sum of every block's two frame sizes
```

`data0 = 136 + h[6] + (h[1] & NAMEMASK) + h[2] + h[3] + h[11]`.
The streams follow each other with no padding and no alignment. Writers write nothing after the last block body. Current readers do not refuse trailing bytes, but they ignore them.

## The header

| word | meaning |
|---|---|
| `h[0]` | number of sequences |
| `h[1]` | flags in bits 49 to 63, and the names stream length in bits 0 to 48 (see below) |
| `h[2]` | length of the index chunks stream |
| `h[3]` | length of the block table stream |
| `h[4]` | number of rows |
| `h[5]` | bits 0 to 31: the lengths codec id. Bits 32 to 63: under `DBFLAG`, the low 32 bits of the features hash, otherwise 0 |
| `h[6]` | length of the labels stream |
| `h[7]` | number of labels in the label table |
| `h[8]` | block size in rows. `UINT64_MAX` means one block. Under `VARFLAG` it is the row threshold at which a block closes. A reader treats 0 as 512, but no writer writes 0 |
| `h[9]` | bits 0 to 31: the labels codec id. Bits 32 to 63: under `DBFLAG`, the high 32 bits of the features hash, otherwise 0 |
| `h[10]` | number of index chunks |
| `h[11]` | length of the index table stream |
| `h[12]` | the tail predictor: the median of the first 64 sequence totals. It is written even when `TAILFLAG` is off, and is 0 below 64 sequences |
| `h[13]` | bits 0 to 7: the tail skip, always 64. Bits 8 to 63: the low 56 bits of a binding hash, the row hash, or under `ALNFLAG` the alignment hash |
| `h[14]` | reserved. Always 0, covered by the metadata hash, and not checked by readers |
| `h[15]` | metadata hash |

### Flags in h[1]

| bit | name | meaning |
|---|---|---|
| 63 | `TOKFLAG` | names stream is the tokeniser encoding |
| 62 | `INTFLAG` | names stream is the integer-delta encoding |
| 61 | `RUNFLAG` | each sequence's run count heads its own lengths in the lengths frame (headed, grouped and interleaved layouts) |
| 60 | `DROPPED` | names were deliberately not stored. A reader needs `--names FILE` (or `--raw`, which prints ordinals) |
| 59 | `ALNFLAG` | rows are coded against a reference through an alignment (`--aln`) |
| 58 | `DBFLAG` | labels are featureIDs of a `features.tsv` (`--db`) |
| 57 | `VARFLAG` | blocks hold whole sequences, vary in row count, and each block table record carries a sequence count |
| 56 | `LENFLAG` | lengths the reference implies may be stored as 0 (always set together with `ALNFLAG`) |
| 55 | `TAILFLAG` | the last run length of every sequence from the 64th on is coded against `h[12]` |
| 54 | (free) | never set |
| 53 | (free) | never set |
| 52 | `CARTFLAG` | some featureIDs are `2^26 +` a vocabulary rank (needs `DBFLAG`) |
| 51 | `ILVFLAG` | labels interleaved with lengths in one frame. The labels frame is empty |
| 50 | `FLDFLAG` | names stream is the field-delta encoding |
| 49 | `GRPFLAG` | each sequence's lengths are grouped by label |

`NAMEMASK` is every bit not in the table above, which is bits 0 to 48 plus the free bits 53 and 54.
Both readers compute the names length as `h[1] & ~FLAGS`, where `FLAGS` is the 13 named bits. A set bit 53 or 54 therefore becomes part of the names length (see [Compatibility](#compatibility)).
At most one of `TOKFLAG`, `INTFLAG` and `FLDFLAG` is set. If none is set and the names stream is non-empty, the names stream is zstd text.

### The metadata hash

`h[15]` is FNV-1a 64 over the 128 header bytes with `h[15]` set to 0 (and `h[1]` with its flags), continued over the labels, names, index chunks, block table and index table streams, in that order.
Block bodies are not covered. Each block has its own hash in the block table.

### The row hash and the binding hash in h[13]

Without `ALNFLAG`, `h[13] >> 8` is the low 56 bits of the **row hash**. It binds a track that later serves as a `--ref`. The row hash starts at `0xcbf29ce484222325`. For every row in order it folds in four 64-bit words, each as `h = (h ^ word) * 0x100000001b3`: the FNV-1a of the name, then start, then end, then the FNV-1a of the label (`kst::rowmix`).
Under `ALNFLAG`, `h[13] >> 8` is the low 56 bits of the **alignment hash**, the XOR of three parts:
- the alignment's segment hash (`Aln::seghash`)
- the reference's row hash (computed over the reference BED, or taken from a `.kst` reference's `h[13] >> 8`)
- with `--aln-names`, the FNV-1a over each line followed by a 0 byte

The segment hash is the wrapping sum over query names of `fnv(name ‖ 0x00 ‖ for each segment: q u64, t u32, n u32, tp u64, rev i32)`, with segments sorted by `(q, t, tp, rev, n)`. The Rust hashes these as little-endian integers. The C hashes the raw struct bytes **(uncertain: whether struct padding can differ between the two was not checked; the test suite shows the two agree on x86-64)**.

### The features hash (under DBFLAG)

The features hash is FNV-1a 64 over the raw bytes of `features.tsv`. Under `CARTFLAG` the hash is continued over the raw bytes of `hierarchy.tsv`. A reader given `-d` refuses a track whose stored hash differs.

## The streams

### Labels (h[6] bytes)

The labels stream is a zstd frame of the label table: each label followed by `\n`, in first-seen order, with `h[7]` labels. A label's id is its 0-based position in the table.
The stream is empty only when there are no rows, and then `h[7] = 0`.
Under `DBFLAG` the labels are decimal featureIDs: 0 means an all-`novel` tuple, `1` to `2^26 - 1` are rows of `features.tsv`, and under `CARTFLAG` a value of `2^26` or more is `2^26 +` the tuple's mixed-radix rank over each set's sorted vocabulary from `--hier` (`db.rs`, `cart_id`).
Under `ALNFLAG` the reference's labels come first, in the reference's table order, and the label `=` means "the label the reference gives this row".

### Names (h[1] & NAMEMASK bytes)

The names give one name per sequence, in sequence order. The names stream takes one of five forms:

- **empty, `DROPPED` clear**: the names are the ordinals `0`, `1`, ..., `h[0] - 1`. The writer stores nothing when the names are exactly these.
- **empty, `DROPPED` set**: the names are absent, and the reader must be given them.
- **zstd text**, no name flag: a zstd frame of every name followed by `\n`.
- **integer deltas**, `INTFLAG`: a zstd frame of `h[0]` varints. Each varint is the zigzag of `name_i - name_{i-1}`, with `name_{-1} = 0`. It is used only when every name is a decimal integer up to `2^50` without a leading zero. A reader refuses a value past `2^50` and bytes left over.
- **tokeniser**, `TOKFLAG`: a chunked htscodecs `tok3` blob. All u32 are little-endian:

  ```
  u32 chunk          names per chunk at encode time (250000, 1000000, 3000000 or 9000000); informative
  u32 nchunk
  u32 raw[nchunk]    each chunk's decoded length in bytes
  u32 coded[nchunk]  each chunk's tok3 blob length
  blobs              tok3_encode_names(level 9, no arith) output, back to back
  ```

  Each chunk holds whole names. `tok3` returns names separated by NUL, and a reader turns every NUL into `\n`. The concatenated chunks are the names, each followed by `\n`.
- **field deltas**, `FLDFLAG`: `varint rawlen`, then a bzip3 frame (codec `0x20000 | 64`, see Codecs) that decodes to `rawlen` bytes `v`. `v` is laid out as follows:

  ```
  varint nf                1 to 16 fields
  varint sep               the separator byte, ':' or '-'
  nf times:
    varint type            0 constant, 1 decimal, 2 hex
    type 0: varint len, then len bytes, the field's text in every name
    type 2: varint width   even, 2 to 64: the field's length in hex digits
  then for each name, for each field in order:
    type 1: varint zigzag(value - previous name's value), previous starting at 0
    type 2: width/2 bytes, each two lowercase hex digits, high nibble first
  ```

  A name is its fields joined by `sep`. A decimal field has 1 to 18 digits and no leading zero. A reader refuses `rawlen` of 0 or past `2^40`, a type above 2, a bad width, a negative value, and bytes left over.

After decoding, the name count must equal `h[0]`. No two names may be equal: the writer refuses repeats, and a reader that builds a name table refuses them.

### Index chunks (h[2] bytes) and index table (h[11] bytes)

The **sequence index** is one varint per sequence, its run count (number of rows), in sequence order.
It is cut into chunks of 65536 sequences (`ICHUNK`, the last chunk shorter). Each chunk is its own zstd frame, and the frames are concatenated into the index chunks stream.
The **index table** is a zstd frame of `h[10]` pairs of varints, one pair per chunk: the chunk's frame length in bytes, and the number of rows before the chunk's first sequence.
Without `RUNFLAG`, `h[10] = ceil(h[0] / 65536)`.
Under `RUNFLAG` the run counts live in the lengths frame instead. Then the index chunks stream is empty, `h[2] = 0` and `h[10] = 0`, and the index table is a zstd frame of zero bytes.
A sequence's first row is the row count before its chunk plus the run counts before it in the chunk. Without `VARFLAG`, that row lies in block `row / h[8]` at offset `row % h[8]`.

### Block table (h[3] bytes)

The block table is a zstd frame of one record per block, in file order. Each record is 7 varints, or 8 under `VARFLAG`:

| # | field |
|---|---|
| 1 | lengths frame size in bytes |
| 2 | labels frame size in bytes (0 under `ILVFLAG`) |
| 3 | decoded size of the labels region |
| 4 | decoded size of the lengths region. Under `LENFLAG` this is `size << 2 \| spelling` (see the `--aln` section) |
| 5 | rows in the block |
| 6 | start coordinate of the block's first row within its sequence (0 when the block opens a sequence) |
| 7 | FNV-1a 64 of the block body (the lengths frame then the labels frame) |
| 8 | `VARFLAG` only: sequences in the block |

Block bodies follow `data0` in the same order. Each body is the lengths frame immediately followed by the labels frame.
Under `ILVFLAG` the single lengths frame decodes to field 4 plus field 3 bytes.

Block sizes:
- **fixed** (no `VARFLAG`): every block but the last has exactly `h[8]` rows, and the last has 1 to `h[8]`. A block may end mid-sequence.
- **whole sequences** (`VARFLAG`): a block closes at the first sequence start at or after `h[8]` rows, so blocks hold whole sequences. There are at least two blocks, and field 8 sums to `h[0]`.

`RUNFLAG` without `VARFLAG` requires exactly one block.

## Codecs

Each frame is decoded by the codec id in `h[5]` (lengths frames) or `h[9]` (labels frames). One id applies to every block.
The decoded length is always known from the block table, and the decoder must produce exactly that many bytes. A reader dispatches on the id as follows:

| id | codec | frame |
|---|---|---|
| `0` | zstd | one zstd frame |
| has `0x80000` | libzpaq 7.15, `0x80000 \| MB` | `u32 np`, then `np` pairs `(u32 raw, u32 coded)`, then the `np` zpaq blocks. Each piece compresses on its own, with method `5` or the karyozip ZPAQL model `ilv` (`c/zpaq/kzzpaq.cpp`). Each zpaq block names its own model, so a decoder does not need to know which. MB (64, `ZPQMB`) is the piece size the writer cut at |
| has `0x40000` | LZMA2, raw (no `.xz` container), liblzma 5.8.4 | `lzma_raw_buffer_encode` with the settings pinned in `c/xz/kzlzma.c`: 64 MiB dictionary, lc 3, lp 0, pb 0, normal mode, nice_len 273, BT4 |
| has `0x20000` | bzip3 1.5.4, `0x20000 \| MB` | a `bz3_compress` frame with block size `MB << 20` (MB is 1 to 256). The vendored `bz3_compress_cuts` may end blocks early, at sequence starts. The frame records each block's size, so an unmodified `bz3_decompress` reads it |
| has `0x10000` | htscodecs 1.6.7 adaptive arithmetic coder | `arith_compress_to(in, order = id & 0xffff)` |
| otherwise | htscodecs 1.6.7 rANS 4x16 | `rans_compress_to_4x16(in, order = id)` |

Readers test the bits in that order: zpaq, LZMA, bzip3, arith, then rANS. The writer never combines two of these bits.

The rANS and arith ids pass htscodecs flags through. These are `0x01` order 1, `0x04` 32-way (rANS), `0x08` stripe, `0x40` RLE and `0x80` pack, and for rANS also `0x100` (`X_SIMD_AUTO`) and a stripe count in bits 8 to 15 with `0x08`.
Plain rANS order 0 is written as `0x100`, so that it is not the zstd id 0. The writer's candidates are listed in `KCAND` in `writer.rs`.
Both htscodecs frames start with their own format byte, so a decoder only needs to know "rANS" or "arith".
A reader refuses an id with bit 31 set. It also refuses a rANS or arith id with `0x10` (`X_NOSZ`), because that frame does not record its size and so cannot be checked.

The tokeniser names use htscodecs `tok3`. The field-delta names use bzip3 with a 64 MB block. The `KSTAL1` streams use the same ids, offset by 3 (see below).

## Decoding a block

1. Read the body. Check its FNV-1a against field 7, then decode the lengths frame and, unless `ILVFLAG`, the labels frame. This gives a **lengths region** (field 4 bytes) and a **labels region** (field 3 bytes).
2. **`ILVFLAG`**: the one frame holds, for each sequence in the block, `varint k`, then `k` pairs `(varint flip(label), lengths value)`, where `flip(v) = v < 128 ? 127 - v : v`, which is its own inverse.
   Split it back into a headed lengths region (each `k` followed by its `k` lengths values) and a labels region (`unflip` of each label). The split must give exactly fields 4 and 3 bytes, and the `k` must sum to the block's rows.
3. **`RUNFLAG`**: the lengths region is, for each sequence in the block, `varint k` (its run count, 1 or more) followed by its `k` lengths values. Lift the `k` out into a per-block run table and keep the bare lengths.
   The number of heads must equal the block's sequence count: `h[0]` for a single block, field 8 under `VARFLAG`. Their sum must equal the block's rows.
4. **`GRPFLAG`** (needs `RUNFLAG`, excludes `ILVFLAG`): within each sequence the lengths are in label order. They are sorted stably by label id, ties by original position. The labels region is in row order.
   To undo it, for each sequence sort `(label id, position)` and give the j-th length to the j-th position of that order.
5. You now have, for every row of the block in order, one lengths value and one label id (a varint index into the label table).

**Plain layout** (no `RUNFLAG`): the lengths region is one varint per row and the labels region is one varint per row. Run counts come from the sequence index.

### Lengths values

A lengths value is normally the run length, which is 1 or more. Two exceptions apply:

- **Tail predictor** (`TAILFLAG`): for the last row of sequence `i` with `i >= h[13] & 0xff` (which must be 64), the value is `zz(total - h[12])`, where `total` is the sequence's end coordinate. The run length is `h[12] + unzz(value) - start`, and it must be positive.
  The writer turns this on when the first 64 totals lie close to their median: the sum of `|total - median| * 12` is at most `median * 64`. That threshold is writer policy, not format.
- **Zero** is a sentinel, allowed only under `LENFLAG` in a block whose spelling allows it (see the next section). Anywhere else a zero length is corrupt, unless the reader runs with `--raw`, which prints it as 0.

## Alignment-informed tracks (ALNFLAG, LENFLAG)

The track is coded against a reference track through a PAF alignment. Neither of these is stored in the file; both are bound by the alignment hash in `h[13]`.
From each PAF record with a `cg:Z:` CIGAR, only the `M`, `=` and `X` operations become **segments** `(q, t, n, tp, rev)`.

- A **label** `=` means the reference's label for the row. That label is the majority of the reference labels at query positions `start + len/4`, `start + len/2` and `start + 3*len/4`: if the first two agree, that label; otherwise the first if it matches the third, else the second. A row shorter than 2 bases uses the label at `start`.
- A query position maps to the reference through the last segment with `q <= pos < q + n`: `tpos = tp + (pos - q)`, or on the reverse strand `tp + (q + n - 1 - pos)`. The reference run covering `tpos` gives the label.
- A **length 0** means "up to where the reference run under `start` ends". On the forward strand that is `q + (run_end - tp)`. On the reverse strand it is `q + n - (run_start - tp)`. In both cases the run must lie within the segment.

Each block stores its lengths in one of three **spellings**, chosen per block as the smallest frame. The spelling is the low two bits of block table field 4:

| value | spelling |
|---|---|
| 0 | a sentinel wherever the reference implies the length |
| 1 | plain lengths, no zero allowed |
| 2 | sentinels only where the segment under the row is on the forward strand |

Value 3 is never written, and readers do not refuse it **(uncertain whether that is intended)**.
Exact segment semantics, including how the PAF is parsed, are in `rust/src/aln.rs` (`load_paf`, `segat`, `end`, `row`). Reimplementers should treat that file as normative.

## Validation a reader must perform

Both readers refuse a file, with a message, on each of the following. A reimplementation should do the same.

- The magic is not `KST00001`, or the file is shorter than 136 bytes.
- Any stream length (`h[1] & ~FLAGS`, `h[2]`, `h[3]`, `h[6]`, `h[11]`) exceeds the file size, or `136 +` their sum does.
- The metadata hash does not match.
- A labels stream with `h[7] = 0`, a label count without a stream, more labels than lines, or more lines than `h[7]`.
- More than one of `TOKFLAG`, `INTFLAG` and `FLDFLAG`.
- `GRPFLAG` with `ILVFLAG` or without `RUNFLAG`, or `CARTFLAG` without `DBFLAG`.
- `TAILFLAG` with `h[13] & 0xff != 64`.
- Codec id: bit 31 set, or an `X_NOSZ` rANS or arith id.
- Block table, for each record:
  - The frame sizes run past the file.
  - Rows exceed `h[4]`.
  - Rows are positive but a decoded size is 0.
  - A decoded size exceeds `16 * h[4] + 64`.
  - Under `VARFLAG`, a sequence count of 0 or above the block's rows.
- Over the block table as a whole:
  - Without `VARFLAG`, the blocks are not `h[8]` rows each with a last block of 1 to `h[8]`.
  - The rows do not sum to `h[4]`.
  - Under `VARFLAG`, the sequence counts do not sum to `h[0]`, or there are fewer than two blocks.
  - Under `RUNFLAG` without `VARFLAG`, more than one block.
- Index table: a chunk runs past the stream, the lengths do not sum to `h[2]`, the row starts are out of order or the first is not 0, or (without `RUNFLAG`) `h[10] != ceil(h[0] / 65536)`.
  Each decoded chunk must hold exactly its count of run counts, each below `2^32`, summing to the next chunk's start.
- Block body: the hash does not match, a frame does not decode to exactly its stated size, a run count of 0 or past the block's rows, the run counts do not match, the interleaved or grouped frames do not line up, a label id past the table (when the table is non-empty), a zero length outside a sentinel block, or a non-positive tail length.
- Names: the decoded count differs from `h[0]`, or any malformation listed under Names.
- Companions: reading an `ALNFLAG`, `DBFLAG` or `LENFLAG` track without its companion, or with a companion whose hash differs, unless `--raw`.
- After decoding, the rows emitted differ from `h[4]`.

Several 32-bit limits are part of the current readers: a sequence's run count and its offset within its first block, the block number of a fixed-size block, the run-count prefix sums under `RUNFLAG`, and at most `2^32` sequences for name lookup and the `.kzi`. The writer refuses a sequence whose first block, offset or run count does not fit 32 bits. It does not check the others **(uncertain whether every reader limit is matched by a writer check)**.

## The .kzi name index

`karyozip index -i TRACK` writes `TRACK.kzi` beside the track and never changes the track. The index maps a name to its ordinal without decoding the names stream.

```
counts      2^k bytes       entries in each bucket (0 to 255)
super       ceil(2^k/1024) u64   entries before bucket 1024*g
entries     ceil(n*(32+w)/8) bytes, a little-endian bitstream
pad         8 zero bytes
trailer     5 x u64: n, k, bind, pay, "KSTIDX01" (0x313058444954534b)
```

- `n = h[0]`, and `w` is the bit width of `n - 1` (at least 1, at most 32). `pay` is the byte count of everything before the trailer, which must equal `2^k + 8*ceil(2^k/1024) + ceil(n*(32+w)/8) + 8`.
- `h64b(name)`: start at `h = 0x9e3779b97f4a7c15`. For each byte `c`, do `h ^= c`, then `h *= 0xff51afd7ed558ccd`, then `h ^= h >> 32`.
- A name's bucket is the top `k` bits of `h64b` (bucket 0 when `k = 0`). Its 32-bit fingerprint is `(h64b ^ ordinal * 0x9e3779b97f4a7c15) & 0xffffffff`. Binding the fingerprint to the ordinal makes a damaged entry fail to match.
- Each entry is `32 + w` bits, least significant bit first: the fingerprint, then the ordinal. Entries are sorted by `(bucket, fingerprint, ordinal)`.
- `k` starts at the smallest value with `2^(k+2) >= n` (at most 30) and grows until no bucket holds more than 255 entries.
- `bind = h[15] XOR namesrc`. `namesrc` is 0, unless the index was built under `--names FILE`. Then it is FNV-1a 64 over the file's path, continued over its `st_size`, `st_ino`, `st_dev`, `st_mtime`, `st_mtime_nsec`, `st_ctime` and `st_ctime_nsec`, each as u64.

A reader uses the index only when all of these hold: the magic matches, `n == h[0]`, `k < 40`, `pay` fits the file, `bind` matches the track and names source, and `pay` matches the formula. Otherwise it silently falls back to decoding the names.
A lookup reads the name's bucket and keeps the entries whose fingerprint matches its ordinal.
**Lookup is probabilistic**: with exactly one match and the names not yet decoded, the ordinal is returned without comparing the name. So a name absent from the track can be answered as a present one, about 4 times in 2^32 lookups.
With several matches, the names are decoded and compared.

## The KSTAL1 packed alignment

`karyozip aln-pack` writes the segments of a PAF (only what `--aln` reads) so that many tracks can share one alignment.

```
magic     6 bytes  "KSTAL1"   (no newline is written; the C's ALNMAGIC string has one, but only 6 bytes are used)
refhash   u64      FNV-1a 64 over each reference contig name followed by '\n', in reference order
10 times: varint codec, varint rawlen, varint codedlen, codedlen bytes
```

Stream codec values:
- `0`: the stream is empty.
- `1`: zstd (level 19).
- `2`: the tokeniser blob of the Names section.
- `c + 3`: payload codec id `c` from the Codecs table. aln-pack tries only the rANS and arith candidates of `-e low`.

The decoded stream must be `rawlen` bytes. The ten streams are:

| # | stream | content |
|---|---|---|
| 0 | names | query names, each followed by `\n`, only those with segments |
| 1 | counts | varint segment count per query (1 or more) |
| 2 | deltas | a gap code per segment that does not start a target run |
| 3 | lengths | varint `n` per segment (1 to `2^32 - 1`) |
| 4 | targets | per target run: varint `t*2 + rev` (contig index into the reference, strand), varint run length in segments |
| 5 | qstarts | varint `q` of each target run's first segment |
| 6 to 9 | tpos0..3 | byte planes, least significant first, of the 32-bit `tp` of each target run's first segment |

A **target run** is a maximal stretch of a query's segments, in sorted order, with the same `t` and strand.
A segment that is not first in its run follows the segment before it (`pq`, `pt` that segment's end on the query and target, `ps` its `tp`). It is coded by a gap pair `(dq, dt)`:
- `q = pq + dq`
- `tp = pt + dt` on the forward strand, or `tp = ps - n - dt` on the reverse strand

The gap code is one varint `v`:
- `0` means `(0, 0)`.
- `3` is followed by `zz(dq)` and `zz(dt)`.
- `v % 3 == 1` means `(v/3 + 1, 0)`.
- `v % 3 == 2` means `(0, v/3 + 1)`.

Any other `v` is corrupt.

A reader refuses the following:
- a reference hash that differs from the given reference
- a stream that overruns the file
- a repeated query name
- a target index past the reference
- an empty target run
- planes of unequal length
- any byte left over in any stream

## Compatibility

Plainly stated:

- **The magic is not a version.** `KST00001` has not changed while the format gained `VARFLAG` (whole-sequence blocks, an 8th block table varint), `GRPFLAG`, `ILVFLAG`, `FLDFLAG`, `CARTFLAG`, the zpaq and LZMA codecs and bzip3 cut lists. A file from an older karyozip may therefore carry the same magic as one this version cannot write, and the reverse holds too. No release has been tagged, so no file is guaranteed to be readable by a later karyozip.
- **Unknown flags are not detected.** Bits 53 and 54 are free. Both readers compute the names length as `h[1]` minus the known flags, so a file with bit 53 or 54 set (say, from a newer writer) is refused as `corrupt header: the stream sizes do not fit the file`, not as "written by a newer karyozip". A new codec id fails later, as a frame decode error.
- **`h[14]` is unused** and is the only spare header word. It is covered by the metadata hash but not checked.
- C and Rust agreement is tested byte for byte. That cannot catch a mistake both share, and there are no golden files from earlier versions in `tests/`.

Open decisions for the maintainer:

1. **Bump the magic or add a version field.** One option is to bump the magic (`KST00002`) whenever the layout changes, and keep readers for older magics or refuse them by name. The other is to reserve a field now, for example `h[14]` as a format minor version or a "required features" bit mask, or bits 53 and 54. Readers would then refuse unknown bits and versions explicitly.
2. **Freeze at the first tag.** Decide whether the format is frozen when 0.2.0 (or 1.0) is tagged. If so, add golden `.kst`, `.kzi` and `KSTAL1` files from that tag to `tests/` and decode them in CI.
3. **Narrow the names length.** Decide whether readers should refuse any set bit outside the known flags and bits 0 to 48, instead of folding it into the names length.
4. **Byte order.** Decide whether to make the C write little-endian explicitly, so that the format does not depend on the host.
5. **Undefined values.** Decide what a reader does with spelling value 3, with `ILVFLAG` without `RUNFLAG` (the writer never produces it, and readers do not refuse it), and with trailing bytes after the last block.
