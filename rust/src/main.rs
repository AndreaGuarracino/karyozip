mod ffi;
mod kst;
mod reader;
mod writer;

use std::io::{BufWriter, Read, Write};

fn die(msg: &str) -> ! {
    eprintln!("karyozip: {msg}");
    std::process::exit(1);
}
pub fn bail(msg: &str) -> ! { die(msg) }

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.iter().any(|a| a == "--version") {
        println!("karyozip {} (rust), format KST00001", env!("CARGO_PKG_VERSION"));
        println!("  zstd       {}", ffi::zstd_version());
        println!("  libbz3     {}", ffi::bzip3_version());
        return;
    }
    // The codec probe tests/diff.sh uses.
    if args.len() == 3 && args[1] == "codecs" { return probe(&args[2]); }

    let mut opt = writer::Options::default();
    let mut tmpdir = std::env::var("TMPDIR").unwrap_or_else(|_| "/tmp".into());
    let mut cmd = None;
    let mut i = 1;
    while i < args.len() {
        let a = args[i].as_str();
        let mut val = || -> String {
            i += 1;
            args.get(i).cloned().unwrap_or_else(|| die(&format!("{a} needs a value")))
        };
        match a {
            "--block" => {
                let v = val();
                opt.block = if v == "all" { u64::MAX }
                            else { v.parse().unwrap_or_else(|_| die("--block is a number of runs, 1 or more, or \"all\"")) };
                if opt.block == 0 { die("--block is a number of runs, 1 or more, or \"all\""); }
            }
            "--level" => {
                let v = val();
                opt.level = v.parse().unwrap_or_else(|_| die("--level is a zstd level, 1 to 22"));
                if !(1..=22).contains(&opt.level) { die("--level is a zstd level, 1 to 22"); }
            }
            "--names" => { if val() == "keep" { opt.keep_names = true; } }
            "--tmp-dir" => tmpdir = val(),
            "--id-codec" => { opt.id_codec = Some(parse_codec(&val())); }
            "--len-codec" => { opt.len_codec = Some(parse_codec(&val())); }
            "-t" | "--threads" => {
                let v: usize = val().parse().unwrap_or_else(|_| die("--threads is a number of 1 to 1024"));
                if !(1..=1024).contains(&v) { die("--threads is a number of 1 to 1024"); }
                // Explicit, so RAYON_NUM_THREADS cannot reach it.
                rayon::ThreadPoolBuilder::new().num_threads(v).build_global().ok();
            }
            _ if a.starts_with('-') && a.len() > 1 => die(&format!("unknown flag {a}")),
            _ => { if cmd.is_none() { cmd = Some(a.to_string()); } }
        }
        i += 1;
    }
    match cmd.as_deref() {
        Some("encode") => encode(opt, std::path::Path::new(&tmpdir)),
        Some("decode") => {
            let f = args.iter().skip(1).find(|a| !a.starts_with('-') && *a != "decode")
                        .unwrap_or_else(|| die("decode needs a track"));
            let mut r = reader::Reader::open(f);
            let stdout = std::io::stdout();
            let mut o = BufWriter::new(stdout.lock());
            r.decode(&mut o);
            o.flush().unwrap_or_else(|e| die(&e.to_string()));
        }
        Some(c) => die(&format!("bad command {c}")),
        None => {
            eprintln!("usage: karyozip encode [--block N|all] [--level N] [--names keep] < in.bed > out.kst");
            std::process::exit(1);
        }
    }
}

fn parse_codec(s: &str) -> i32 {
    let v = if let Some(h) = s.strip_prefix("0x") {
        i64::from_str_radix(h, 16)
    } else { s.parse::<i64>() };
    match v {
        Ok(v) if (0..=0x7fff_ffff).contains(&v) => v as i32,
        _ => die("--id-codec and --len-codec take a codec id, not any number"),
    }
}

fn encode(opt: writer::Options, tmpdir: &std::path::Path) {
    let mut w = writer::Writer::new(opt, tmpdir)
        .unwrap_or_else(|e| die(&format!("cannot create a temp file under {}: {e}", tmpdir.display())));
    // A slab at a time, like the C. Reading the whole input first would cost
    // the input's size in memory, 2.3 GB on the benchmark track.
    const SLAB: usize = 1 << 24;
    let mut slab = vec![0u8; SLAB];
    let mut have = 0usize;
    let mut stdin = std::io::stdin().lock();
    let mut eof = false;
    loop {
        while !eof && have < SLAB {
            match stdin.read(&mut slab[have..]) {
                Ok(0) => eof = true,
                Ok(n) => have += n,
                Err(e) => die(&e.to_string()),
            }
        }
        let end = match slab[..have].iter().rposition(|&c| c == b'\n') {
            Some(i) => i + 1,
            None => {
                if eof && have > 0 { die("input does not end with a newline"); }
                if !eof && have == SLAB { die("a BED line is longer than the read buffer"); }
                break;
            }
        };
        for line in slab[..end].split(|&c| c == b'\n') {
            if line.is_empty() { continue; }
            row(&mut w, line);
        }
        slab.copy_within(end..have, 0);
        have -= end;
        if eof && have == 0 { break; }
    }

    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    let st = w.close(&mut out).unwrap_or_else(|m| die(&m));
    out.flush().unwrap_or_else(|e| die(&e.to_string()));
    eprintln!("codec: len=0x{:x} id=0x{:x} from {} sample blocks, {}",
              st.len_codec, st.id_codec,
              if st.blocks == 0 { 0 } else if st.blocks > 2 { 3 } else if st.blocks > 1 { 2 } else { 1 },
              if st.seekable { "seekable" } else { "one block" });
    eprintln!("encoded rows={} reads={} blocks={} names={}B{} index={}B offs={}B lens={}B labs={}B raw={}B",
              st.rows, st.reads, st.blocks, st.names,
              if st.ordinal { " (ordinal, not stored)" }
              else if st.tokused { " (tokenised)" }
              else if st.intused { " (int deltas)" } else { "" },
              st.index, st.offs, st.lens, st.labs, st.raw);
}

fn row(w: &mut writer::Writer, line: &[u8]) {
    if line.contains(&0) { die("a NUL byte inside a row"); }
    let mut it = line.splitn(4, |&c| c == b'\t');
    let (name, s, e, label) = match (it.next(), it.next(), it.next(), it.next()) {
        (Some(a), Some(b), Some(c), Some(d)) => (a, b, c, d),
        _ => die("bad line"),
    };
    if s.is_empty() || e.is_empty() { die("empty start or end column"); }
    let st = bedcoord(s, "start is not a number");
    let en = bedcoord(e, "end is not a number");
    if en < st { die("end is before start"); }
    if en == st { die("end equals start: a run must cover at least one base"); }
    if name.is_empty() { die("empty sequence name"); }
    w.row(name, en - st, label, st).unwrap_or_else(|m| die(&m));
}

fn bedcoord(b: &[u8], what: &str) -> u64 {
    let mut v: u64 = 0;
    for &c in b {
        if !c.is_ascii_digit() { die(what); }
        if v > (u64::MAX - (c - b'0') as u64) / 10 { die("a coordinate does not fit in 64 bits"); }
        v = v * 10 + (c - b'0') as u64;
    }
    v
}

fn probe(path: &str) {
    let mut src = Vec::new();
    std::fs::File::open(path).unwrap().read_to_end(&mut src).unwrap();
    let mut out = std::io::stdout().lock();
    for lvl in [1, 3, 17, 19, 22] {
        let z = ffi::zstd(&src, lvl);
        writeln!(out, "zstd:{lvl}\t{}\t{:016x}", z.len(), kst::fnv(kst::FNV0, &z)).unwrap();
    }
    for c in writer::CAND_FAST {
        match ffi::paycomp(&src, c) {
            Some(v) => {
                let ok = ffi::paydecomp(&v, c, src.len()).as_deref() == Some(&src[..]);
                writeln!(out, "codec:0x{c:x}\t{}\t{:016x}\troundtrip={}",
                         v.len(), kst::fnv(kst::FNV0, &v), ok).unwrap();
            }
            None => writeln!(out, "codec:0x{c:x}\trefused").unwrap(),
        }
    }
    match ffi::tok_encode(&src) {
        Some(v) => writeln!(out, "tok\t{}\t{:016x}\troundtrip={}", v.len(),
                            kst::fnv(kst::FNV0, &v),
                            ffi::tok_decode(&v).as_deref() == Some(&src[..])).unwrap(),
        None => writeln!(out, "tok\trefused").unwrap(),
    }
}
