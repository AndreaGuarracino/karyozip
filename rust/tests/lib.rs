//! The library writes the same bytes as the CLI and refuses the same contract violations: the Rust half of tests/libtest.c
use karyozip::{ffi, writer::{self, Options, Writer}, Error, ErrorKind, KzResult};
use std::io::Write;

/// One BED through `Writer::row`, as an annotator would call it
fn build(bed: &[u8], opt: Options) -> Vec<u8> {
    let dir = std::env::temp_dir();
    build_in(bed, opt, &dir)
}
fn build_in(bed: &[u8], opt: Options, dir: &std::path::Path) -> Vec<u8> {
    let mut w = Writer::new(opt, dir).expect("temp file");
    for line in bed.split(|&c| c == b'\n') {
        if line.is_empty() { continue; }
        let f: Vec<&[u8]> = line.splitn(4, |&c| c == b'\t').collect();
        let st: u64 = std::str::from_utf8(f[1]).unwrap().parse().unwrap();
        let en: u64 = std::str::from_utf8(f[2]).unwrap().parse().unwrap();
        w.row(f[0], en - st, f[3], st).expect("row");
    }
    let mut out = Vec::new();
    w.close(&mut out).expect("close");
    out
}

/// The -t of every CLI child: KZ_TEST_THREADS (the Makefile passes each test's share), else 1; the bytes never depend on it
fn test_threads() -> String {
    match std::env::var("KZ_TEST_THREADS").ok().and_then(|v| v.parse::<u32>().ok()) {
        Some(n) if (1..=1024).contains(&n) => n.to_string(),
        _ => "1".into(),
    }
}

fn cli(args: &[&str], bed: &[u8]) -> Vec<u8> {
    let exe = env!("CARGO_BIN_EXE_karyozip");
    let mut c = std::process::Command::new(exe)
        .arg("-t").arg(test_threads())
        .args(args).arg("compress")
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::null())
        .spawn().expect("spawn");
    c.stdin.take().unwrap().write_all(bed).unwrap();
    c.wait_with_output().unwrap().stdout
}

fn fixture() -> Vec<u8> {
    let mut b = Vec::new();
    for s in 0..40 {
        let mut p = 0u64;
        for k in 0..25 {
            let l = 7 + (s * 13 + k * 5) % 400;
            b.extend_from_slice(format!("r{s}\t{p}\t{}\tL{}\n", p + l, k % 6).as_bytes());
            p += l;
        }
    }
    b
}

#[test]
fn same_bytes_as_the_cli() {
    let bed = fixture();
    for (args, opt) in [
        (vec![], Options::default()),
        (vec!["--block", "4"], Options { block: 4, ..Options::default() }),
        (vec!["--block", "all"], Options { block: u64::MAX, ..Options::default() }),
        (vec!["--names", "drop"], Options { drop_names: true, ..Options::default() }),
    ] {
        assert_eq!(build(&bed, opt), cli(&args, &bed), "differs for {args:?}");
    }
}

/// Every contract violation is refused, not written out, in the CLI's words
#[test]
fn refuses_what_the_header_promises() {
    let big = vec![b'x'; 4096];
    type Case = (Box<dyn Fn(&mut Writer) -> KzResult<()>>, &'static str);
    let bad: Vec<Case> = vec![
        (Box::new(|w| w.row(b"r1", 10, b"A", 5)), "a sequence does not start at 0"),
        (Box::new(|w| { w.row(b"r1", 10, b"A", 0)?; w.row(b"r1", 10, b"A", 50) }), "gap or overlap: start is not the previous end"),
        (Box::new(|w| { w.row(b"r1", 10, b"A", 0)?; w.row(b"r1", 10, b"A", 5) }), "gap or overlap: start is not the previous end"),
        (Box::new(|w| { w.row(b"r1", 10, b"A", 0)?; w.row(b"r2", 10, b"A", 0)?;
                        w.row(b"r1", 10, b"A", 10) }), "a sequence does not start at 0"),   // a sequence returns
        (Box::new(|w| w.row(b"r1", 0, b"A", 0)), "end equals start: a run must cover at least one base"),
        (Box::new(move |w| w.row(&big, 10, b"A", 0)), "name too long"),
        (Box::new(|w| w.row(b"", 10, b"A", 0)), "empty sequence name"),
        (Box::new(|w| w.row(b"r1", 10, b"A\nB", 0)), "a newline inside a name or label"),
        (Box::new(|w| w.row(b"r\0x", 10, b"A", 0)), "a NUL byte inside a row"),
        (Box::new(|w| w.row(b"a\tb", 10, b"A", 0)), "a tab inside a sequence name"),
        (Box::new(|w| { w.row(b"r1", u64::MAX, b"A", 0)?;
                        w.row(b"r1", 10, b"A", u64::MAX) }), "a coordinate does not fit in 64 bits"),
    ];
    let dir = std::env::temp_dir();
    for (i, (f, want)) in bad.iter().enumerate() {
        let mut w = Writer::new(Options::default(), &dir).expect("temp file");
        let e = f(&mut w).expect_err(&format!("case {i} was accepted"));
        assert_eq!((e.message(), e.kind()), (*want, ErrorKind::Input), "case {i}");
    }
    // a repeated name is found at the close
    let mut w = Writer::new(Options::default(), &dir).expect("temp file");
    for n in [&b"r1"[..], b"r2", b"r1"] { w.row(n, 10, b"A", 0).expect("row"); }
    let e = w.close(&mut Vec::new()).err().expect("a repeated name was accepted");
    assert_eq!(e.message(), "a sequence name repeats");
}

/// A directory of its own under cargo's scratch space, so a test can count what is left in it
fn scratch(name: &str) -> std::path::PathBuf {
    let d = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("kz-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).expect("scratch dir");
    d
}
/// karyozip.* temporaries left under `d`
fn leftovers(d: &std::path::Path) -> Vec<String> {
    std::fs::read_dir(d).map(|r| r.filter_map(|e| e.ok()).map(|e| e.file_name().to_string_lossy().into_owned())
                                  .filter(|n| n.starts_with("karyozip")).collect()).unwrap_or_default()
}
fn tmp_msg(d: &std::path::Path) -> String { format!("cannot create a temp file under {} (see --tmp-dir)", d.display()) }

fn feed(w: &mut Writer, bed: &[u8]) {
    for line in bed.split(|&c| c == b'\n').filter(|l| !l.is_empty()) {
        let f: Vec<&[u8]> = line.splitn(4, |&c| c == b'\t').collect();
        let st: u64 = std::str::from_utf8(f[1]).unwrap().parse().unwrap();
        let en: u64 = std::str::from_utf8(f[2]).unwrap().parse().unwrap();
        w.row(f[0], en - st, f[3], st).expect("row");
    }
}

/// A temp directory the writer cannot use is an Err, in the CLI's words, and the process goes on
#[test]
fn a_bad_tmp_dir_is_an_error() {
    let bed = fixture();
    let d = scratch("tmpdir");
    // missing
    let gone = d.join("nosuch");
    let e = Writer::new(Options::default(), &gone).err().expect("a missing tmp dir was accepted");
    assert_eq!((e.message(), e.kind()), (tmp_msg(&gone).as_str(), ErrorKind::Io));
    assert_eq!(e.to_string(), tmp_msg(&gone));
    // moved away between new and close (a rename, since NFS keeps an unlinked open file as .nfs*): the close's body file has nowhere to go
    let sub = d.join("moved"); std::fs::create_dir(&sub).unwrap();
    let mut w = Writer::new(Options { block: 4, ..Options::default() }, &sub).expect("temp file");
    feed(&mut w, &bed);
    std::fs::rename(&sub, d.join("moved.away")).unwrap();
    let e = w.close(&mut Vec::new()).err().expect("close without a tmp dir was accepted");
    assert_eq!((e.message(), e.kind()), (tmp_msg(&sub).as_str(), ErrorKind::Io));
    assert!(leftovers(&d.join("moved.away")).is_empty(), "left: {:?}", leftovers(&d.join("moved.away")));
    // read-only (root writes anyway: then only that it works)
    use std::os::unix::fs::PermissionsExt;
    let ro = d.join("ro"); std::fs::create_dir(&ro).unwrap();
    std::fs::set_permissions(&ro, std::fs::Permissions::from_mode(0o555)).unwrap();
    let root = std::fs::File::create(ro.join("probe")).is_ok();
    if !root {
        let e = Writer::new(Options::default(), &ro).err().expect("a read-only tmp dir was accepted");
        assert_eq!(e.message(), tmp_msg(&ro));
        // writable at new, read-only at close
        std::fs::set_permissions(&ro, std::fs::Permissions::from_mode(0o755)).unwrap();
        let mut w = Writer::new(Options::default(), &ro).expect("temp file");
        feed(&mut w, &bed);
        std::fs::set_permissions(&ro, std::fs::Permissions::from_mode(0o555)).unwrap();
        let e = w.close(&mut Vec::new()).err().expect("close in a read-only dir was accepted");
        assert_eq!(e.message(), tmp_msg(&ro));
        assert!(leftovers(&ro).is_empty(), "left in the tmp dir: {:?}", leftovers(&ro));
    }
    std::fs::set_permissions(&ro, std::fs::Permissions::from_mode(0o755)).unwrap();
    // still alive, and the next writer in the same directory gives the CLI's bytes
    assert_eq!(build_in(&bed, Options::default(), &d), cli(&[], &bed));
    assert!(leftovers(&d).is_empty(), "left in the tmp dir: {:?}", leftovers(&d));
    let _ = std::fs::remove_dir_all(&d);
}

/// An output that fails after `room` bytes, as a full disk does
struct Full { room: usize, got: Vec<u8> }
impl Write for Full {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if self.got.len() + b.len() > self.room {
            return Err(std::io::Error::new(std::io::ErrorKind::StorageFull, "no space left"));
        }
        self.got.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> { Ok(()) }
}

/// A failing output is "write failed", as the CLI says it; no temp file stays, and the process goes on
#[test]
fn a_failing_output_is_an_error() {
    let bed = fixture();
    let d = scratch("full");
    let whole = build_in(&bed, Options { block: 4, ..Options::default() }, &d);
    // nothing, inside the header, inside the tables, inside the body, one byte short
    for room in [0, 8, 100, whole.len() / 2, whole.len() - 1] {
        let mut w = Writer::new(Options { block: 4, ..Options::default() }, &d).expect("temp file");
        feed(&mut w, &bed);
        let mut out = Full { room, got: Vec::new() };
        let e = w.close(&mut out).err().unwrap_or_else(|| panic!("an output full at {room} B was accepted"));
        assert_eq!((e.message(), e.kind()), ("write failed", ErrorKind::Io), "full at {room} B");
        assert!(out.got.len() <= room);
        assert_eq!(out.got[..], whole[..out.got.len()], "what was written is the file's start");
        assert!(leftovers(&d).is_empty(), "left in the tmp dir: {:?}", leftovers(&d));
    }
    let mut out = Full { room: whole.len(), got: Vec::new() };
    let mut w = Writer::new(Options { block: 4, ..Options::default() }, &d).expect("temp file");
    feed(&mut w, &bed);
    w.close(&mut out).expect("an output with room");
    assert_eq!(out.got, whole);
    // a Writer dropped unclosed, rows and all, leaves nothing either
    let mut w = Writer::new(Options::default(), &d).expect("temp file");
    feed(&mut w, &bed);
    drop(w);
    assert!(leftovers(&d).is_empty(), "left in the tmp dir: {:?}", leftovers(&d));
    let _ = std::fs::remove_dir_all(&d);
}

/// Writers one after another in one process, each the CLI's bytes
#[test]
fn writers_in_sequence() {
    let bed = fixture();
    let want = cli(&["--block", "4"], &bed);
    for _ in 0..3 { assert_eq!(build(&bed, Options { block: 4, ..Options::default() }), want, "a later writer in the same process"); }
}

/// TEST_LIGHT=1: serial and light (tests/lib.sh); a test that starts threads of its own returns early, saying so
fn light(what: &str) -> bool {
    let on = std::env::var("TEST_LIGHT").is_ok_and(|v| v == "1");
    if on { eprintln!("skipped (TEST_LIGHT): {what}"); }
    on
}

/// Two writers at once on two threads, each the CLI's bytes; skipped under TEST_LIGHT=1, as it runs two at once
#[test]
fn writers_at_once() {
    if light("two writers at once") { return; }
    let bed = fixture();
    let want = cli(&["--block", "4"], &bed);
    let opt = || Options { block: 4, ..Options::default() };
    let (a, b) = std::thread::scope(|s| {
        let a = s.spawn(|| build(&bed, opt()));
        let b = s.spawn(|| build(&bed, opt()));
        (a.join().unwrap(), b.join().unwrap())
    });
    assert_eq!((a == want, b == want), (true, true), "two writers at once");
}

/// The error type: a message, a kind, and std's Error
#[test]
fn the_error_type() {
    let e = Error::input("name too long");
    assert_eq!((e.message(), e.kind(), e.to_string()), ("name too long", ErrorKind::Input, "name too long".to_string()));
    let b: Box<dyn std::error::Error + Send + Sync> = Box::new(e.clone());
    assert_eq!(b.to_string(), "name too long");
}

/// One file of tests/fixtures.txt, where "@@ NAME" opens each
fn fixfile(name: &str) -> Vec<u8> {
    let (all, mut out, mut on) = (std::fs::read_to_string(format!("{}/../tests/fixtures.txt", env!("CARGO_MANIFEST_DIR"))).unwrap(), String::new(), false);
    for l in all.lines() { if let Some(n) = l.strip_prefix("@@ ") { on = n == name } else if on { out += l; out.push('\n') } }
    out.into_bytes()
}

/// Every codec candidate, and the name tokeniser, round-trips its input
#[test]
fn every_codec_round_trips() {
    for f in ["db.rall.bed", "db.chromosome.bed", "ls.bed"] {
        let mut src = fixfile(f);
        let mut n = 0;
        for c in writer::candlist(true) {
            if let Some(v) = ffi::paycomp(&src, c) {
                n += 1;
                assert_eq!(ffi::paydecomp(&v, c, src.len()).as_deref(), Some(&src[..]), "codec 0x{c:x} on {f}");
            }
        }
        assert!(n > 0, "no codec took {f}");
        if let Some(v) = ffi::tok_encode(&mut src) {
            assert_eq!(ffi::tok_decode(&v).as_deref(), Some(&src[..]), "tokeniser on {f}");
        }
    }
}

/// db::encode over per-set .kst inputs the library wrote: the CLI's bytes from the BEDs; a set missing is an Err
#[test]
fn db_encode_from_per_set_kst() {
    let d = scratch("dbkst");
    let feats_path = d.join("features.tsv");
    std::fs::write(&feats_path, fixfile("db.features.tsv")).unwrap();
    let sets = ["acrocentric", "chromosome", "gene", "region", "repeat", "subtelomeric"];
    let (mut beds, mut ksts) = (Vec::new(), Vec::new());
    for (i, s) in sets.iter().enumerate() {
        let bed = fixfile(&format!("db.{s}.bed"));
        let (b, k) = (d.join(format!("{s}.bed")), d.join(format!("{s}.kst")));
        std::fs::write(&b, &bed).unwrap();
        std::fs::write(&k, build_in(&bed, Options { block: [u64::MAX, 4, 262144][i % 3], ..Options::default() }, &d)).unwrap();
        beds.push(b.display().to_string()); ksts.push(k.display().to_string());
    }
    let mut c = std::process::Command::new(env!("CARGO_BIN_EXE_karyozip"));
    c.arg("-t").arg(test_threads()).arg("compress").arg("--db").arg(&feats_path);
    for b in &beds { c.arg("-i").arg(b); }
    let want = c.stderr(std::process::Stdio::null()).output().expect("spawn").stdout;
    assert!(!want.is_empty());
    let feats = karyozip::db::Feats::load(&feats_path.display().to_string()).expect("features");
    let mut w = Writer::new(Options::default(), &d).expect("temp file");
    let used = karyozip::db::encode(&ksts, &feats, &mut w).expect("encode");
    w.set_feathash(if used { feats.hierhash } else { feats.hash }, used);
    let mut out = Vec::new();
    w.close(&mut out).expect("close");
    assert_eq!(out, want, "per-set .kst through the library");
    // a set missing: the CLI's refusal, returned
    let mut w = Writer::new(Options::default(), &d).expect("temp file");
    let e = karyozip::db::encode(&ksts[..5], &feats, &mut w).expect_err("five sets of six were accepted");
    assert!(e.message().starts_with("--db: give one BED per feature set"), "{}", e.message());
    assert!(leftovers(&d).is_empty());
    let _ = std::fs::remove_dir_all(&d);
}

/// A packed alignment whose tokenised names section is malformed: an Err in the CLI's words from the loader, and the
/// caller goes on (the sentinel after each call), where the token decoder once exited the process
#[test]
fn a_bad_packed_names_section_is_an_error() {
    use karyozip::aln;
    let d = scratch("pkname");
    let vput = |o: &mut Vec<u8>, mut v: u64| { while v >= 0x80 { o.push(v as u8 | 0x80); v >>= 7; } o.push(v as u8); };
    let mut alive = 0;
    let cases: [(&[u8], &str); 4] = [
        (b"", "truncated tokenised names"),                                     // codec 2, raw 0, coded 0
        (b"\0\0\0\0\x05\0\0\0", "corrupt tokenised names header"),           // five chunks, no chunk table
        (b"\0\0\0\0\x01\0\0\0\x04\0\0\0\x64\0\0\0", "corrupt tokenised names: chunks overrun the stream"),
        (b"\0\0\0\0\x01\0\0\0\x04\0\0\0\x03\0\0\0\x01\x02\x03", "name tokeniser decode failed"),
    ];
    for (i, (names, want)) in cases.iter().enumerate() {
        let mut a = aln::Aln::default();
        let mut f = b"KSTAL1".to_vec();
        f.extend_from_slice(&aln::refhash(&a).to_le_bytes());
        vput(&mut f, 2); vput(&mut f, 0); vput(&mut f, names.len() as u64); f.extend_from_slice(names);
        for _ in 1..10 { vput(&mut f, 0); vput(&mut f, 0); vput(&mut f, 0); }   // the other nine streams empty
        let p = d.join(format!("t{i}.pk"));
        std::fs::write(&p, &f).unwrap();
        let e = a.load_aln(&p.display().to_string()).expect_err(&format!("case {i} was accepted"));
        assert_eq!((e.message(), e.kind()), (*want, ErrorKind::Input), "case {i}");
        alive += 1;
    }
    assert_eq!(alive, cases.len(), "the caller is still running");
    let _ = std::fs::remove_dir_all(&d);
}

/// Alignment coordinates that leave 64 bits: each an Err in the CLI's words, and the caller goes on (the sentinel
/// after each call), where an overflow-checking build once panicked and a release build wrapped. Then the boundary
/// cases, a block ending at 2^64 - 1, still accepted, from a PAF and from a packed alignment.
#[test]
fn alignment_coordinates_past_64_bits_are_errors() {
    use karyozip::aln;
    let d = scratch("alnovf");
    let refp = d.join("ref.bed");
    std::fs::write(&refp, "s\t0\t10\tA\n").unwrap();
    let refs = refp.display().to_string();
    let fresh = || { let mut a = aln::Aln::default(); a.load_ref(&refs, aln::Lazy::No).expect("reference"); a };
    const M: u64 = u64::MAX;
    let mut alive = 0;
    // PAF: forward query end, target start, a deletion's target, a reverse CIGAR running before 0
    let pafs: [(String, &str); 5] = [
        (format!("q\t1\t{M}\t{M}\t+\ts\t10\t0\t1\t1\t1\t60\tcg:Z:1M"), "PAF coordinates overflow 64 bits"),
        (format!("q\t9\t{}\t{M}\t+\ts\t10\t0\t2\t1\t1\t60\tcg:Z:1M1I", M - 1), "PAF coordinates overflow 64 bits"),
        (format!("q\t1\t0\t1\t+\ts\t10\t{M}\t{M}\t1\t1\t60\tcg:Z:1M"), "PAF coordinates overflow 64 bits"),
        (format!("q\t1\t0\t1\t+\ts\t10\t{}\t{M}\t1\t1\t60\tcg:Z:1M1D", M - 1), "PAF coordinates overflow 64 bits"),
        ("q\t1\t0\t0\t-\ts\t10\t0\t1\t1\t1\t60\tcg:Z:1M".to_string(), "PAF CIGAR runs before the query start"),
    ];
    for (i, (line, want)) in pafs.iter().enumerate() {
        let p = d.join(format!("o{i}.paf"));
        std::fs::write(&p, format!("{line}\n")).unwrap();
        let mut a = fresh();
        let e = a.load_aln(&p.display().to_string()).expect_err(&format!("PAF case {i} was accepted"));
        assert_eq!((e.message(), e.kind()), (*want, ErrorKind::Input), "PAF case {i}");
        alive += 1;
    }
    // the boundary: blocks ending at 2^64 - 1 on both sides
    let p = d.join("edge.paf");
    std::fs::write(&p, format!("q\t9\t{}\t{M}\t+\ts\t10\t{}\t{M}\t1\t1\t60\tcg:Z:1M\n", M - 1, M - 1)).unwrap();
    let mut a = fresh();
    a.load_aln(&p.display().to_string()).expect("a block ending at 2^64 - 1");
    assert_eq!(a.qs.len(), 1);
    // packed: names, counts, deltas, lengths, targets, qstarts, four tpos planes; each stream zstd (codec 1)
    let vput = |o: &mut Vec<u8>, mut v: u64| { while v >= 0x80 { o.push(v as u8 | 0x80); v >>= 7; } o.push(v as u8); };
    let zz = |v: i64| ((v << 1) ^ (v >> 63)) as u64;
    let packed = |cnt: u64, deltas: &[u8], tgt: u64, qs: u64| -> Vec<u8> {
        let a = fresh();
        let mut f = b"KSTAL1".to_vec();
        f.extend_from_slice(&aln::refhash(&a).to_le_bytes());
        let (mut c, mut l, mut t, mut q) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        vput(&mut c, cnt);
        for _ in 0..cnt { vput(&mut l, 1); }
        vput(&mut t, tgt); vput(&mut t, cnt);
        vput(&mut q, qs);
        for st in [&b"q\n"[..], &c, deltas, &l, &t, &q, &[0], &[0], &[0], &[0]] {
            if st.is_empty() { vput(&mut f, 0); vput(&mut f, 0); vput(&mut f, 0); continue; }
            let z = ffi::zstd(st, 3);
            vput(&mut f, 1); vput(&mut f, st.len() as u64); vput(&mut f, z.len() as u64); f.extend_from_slice(&z);
        }
        f
    };
    let (mut fwd_back, mut zero) = (Vec::new(), Vec::new());
    vput(&mut fwd_back, 3); vput(&mut fwd_back, zz(0)); vput(&mut fwd_back, zz(-2));   // tp = 1 - 2: wraps to 2^64 - 1
    vput(&mut zero, 0);                                                                  // reverse: tp = 0 - 1
    let pks: [(Vec<u8>, &str); 3] = [
        (packed(1, &[], 0, M), "one block at query position 2^64 - 1"),
        (packed(2, &fwd_back, 0, 0), "a forward target position wrapping below 0"),
        (packed(2, &zero, 1, 5), "a reverse target position running below 0"),
    ];
    for (i, (f, what)) in pks.iter().enumerate() {
        let p = d.join(format!("o{i}.pk"));
        std::fs::write(&p, f).unwrap();
        let mut a = fresh();
        let e = a.load_aln(&p.display().to_string()).expect_err(&format!("packed case {i}, {what}, was accepted"));
        assert_eq!((e.message(), e.kind()), ("corrupt packed alignment: a coordinate overflows", ErrorKind::Input), "packed case {i}: {what}");
        alive += 1;
    }
    assert_eq!(alive, pafs.len() + pks.len(), "the caller is still running");
    // the boundary, packed: one block at 2^64 - 2, ending at 2^64 - 1
    let p = d.join("edge.pk");
    std::fs::write(&p, packed(1, &[], 0, M - 1)).unwrap();
    let mut a = fresh();
    a.load_aln(&p.display().to_string()).expect("a packed block ending at 2^64 - 1");
    assert_eq!((a.qs.len(), a.qs[0].1[0].q), (1, M - 1));
    let _ = std::fs::remove_dir_all(&d);
}

/// A features table with its header alone: a novel row encodes as the CLI's bytes, any other label is the CLI's
/// refusal, returned; neither panics
#[test]
fn db_encode_with_a_header_only_table() {
    let d = scratch("dbhdr");
    let tsv = d.join("features.tsv");
    std::fs::write(&tsv, "featureID\tsetA\n").unwrap();
    let feats = karyozip::db::Feats::load(&tsv.display().to_string()).expect("features");
    let (nb, xb) = (d.join("novel.bed"), d.join("x.bed"));
    std::fs::write(&nb, "chr1\t0\t5\tnovel\n").unwrap();
    std::fs::write(&xb, "chr1\t0\t5\tx\n").unwrap();
    let want = std::process::Command::new(env!("CARGO_BIN_EXE_karyozip")).arg("-t").arg(test_threads())
        .arg("compress").arg("--db").arg(&tsv).arg("-i").arg(&nb)
        .stderr(std::process::Stdio::null()).output().expect("spawn").stdout;
    assert!(!want.is_empty());
    let mut w = Writer::new(Options::default(), &d).expect("temp file");
    let used = karyozip::db::encode(&[nb.display().to_string()], &feats, &mut w).expect("a novel row");
    w.set_feathash(if used { feats.hierhash } else { feats.hash }, used);
    let mut out = Vec::new();
    w.close(&mut out).expect("close");
    assert_eq!(out, want, "a header-only table through the library");
    let mut w = Writer::new(Options::default(), &d).expect("temp file");
    let e = karyozip::db::encode(&[xb.display().to_string()], &feats, &mut w).expect_err("a label not in the table was accepted");
    assert_eq!((e.message(), e.kind()), ("a label tuple is not in the features table", ErrorKind::Input));
    drop(w);
    assert!(leftovers(&d).is_empty());
    let _ = std::fs::remove_dir_all(&d);
}
