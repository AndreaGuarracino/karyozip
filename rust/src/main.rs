mod aln;
mod db;
mod ffi;
mod kst;
mod reader;
mod smooth;
mod writer;

use std::io::{BufWriter, Read, Write};

fn die(msg: &str) -> ! {
    eprintln!("karyozip: {msg}");
    std::process::exit(1);
}
pub fn bail(msg: &str) -> ! { die(msg) }

const USAGE: &str = "usage:
  karyozip encode [--block N|all] [--level N] [--slow: smaller, slower]
                  [--names keep|drop|drop:FILE]
  karyozip decode F   --names F|ordinal   --raw
  karyozip region F SEQ | SEQ POS | SEQ START END [HALO]
  karyozip index  F                     appends an index to the track itself
  any command: --threads N, -t N  --tmp-dir DIR  --help, -h  --version";

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.iter().any(|a| a == "--help" || a == "-h") {
        println!("{USAGE}");
        return;
    }
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
    let mut ro = reader::ReadOpts::default();
    let (mut maxgap, mut halo, mut v1) = (1000u64, 200u64, false);
    let (mut dbfeats, mut outdir): (Option<String>, Option<String>) = (None, None);
    let (mut alnpaf, mut alnref): (Option<String>, Option<String>) = (None, None);
    let mut alnnames: Option<String> = None;
    let mut cmd = None;
    let mut pos: Vec<String> = Vec::new();
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
                opt.levelset = true;
            }
            "--names" => {
                let v = val();
                if v == "keep" { opt.keep_names = true; }
                else if v == "ordinal" { ro.ordnames = true; }
                else if v == "drop" { opt.drop_names = true; }
                else if let Some(f) = v.strip_prefix("drop:") {
                    if f.is_empty() { die("--names drop:FILE needs a file"); }
                    opt.names_out = Some(f.to_string());
                } else { ro.names_file = Some(v); }
            }
            "--slow" => opt.slow = true,
            "--sort" => opt.sort = true,
            "--db" => { dbfeats = Some(val()); ro.db = true; }
            "--out" => outdir = Some(val()),
            "--aln" => alnpaf = Some(val()),
            "--ref" => alnref = Some(val()),
            "--aln-names" => alnnames = Some(val()),
            "--raw" => ro.raw = true,
            "--max-gap" => maxgap = val().parse().unwrap_or_else(|_| die("--max-gap needs a number")),
            "--halo" => halo = val().parse().unwrap_or_else(|_| die("--halo needs a number")),
            "--v1" => v1 = true,
            "--tail" => {
                let v = val();
                opt.tail = Some(match v.as_str() { "0" => false, "1" => true,
                                                   _ => die("--tail is 0 or 1") });
            }
            "--tok-chunk" => {
                let v: usize = val().parse().unwrap_or_else(|_| die("--tok-chunk needs a number"));
                if v == 0 { die("--tok-chunk needs a number"); }
                opt.tok_chunk = Some(v.min(9_000_000));
            }
            "--profile" => {}
            "--tmp-dir" => tmpdir = val(),
            "--id-codec" => { opt.id_codec = Some(parse_codec(&val())); }
            "--len-codec" => { opt.len_codec = Some(parse_codec(&val())); }
            "-t" | "--threads" => {
                let v: usize = val().parse().unwrap_or_else(|_| die("--threads is a number of 1 to 1024"));
                if !(1..=1024).contains(&v) { die("--threads is a number of 1 to 1024"); }
                // Explicit, so RAYON_NUM_THREADS cannot reach it.
                rayon::ThreadPoolBuilder::new().num_threads(v).build_global().ok();
            }
            _ if a.starts_with("--sort=") => {
                opt.sort = true;
                mem_arg(&a[7..]);
            }
            _ if a.starts_with('-') && a.len() > 1 => die(&format!("unknown flag {a}")),
            _ => { if cmd.is_none() { cmd = Some(a.to_string()); } else { pos.push(a.to_string()); } }
        }
        i += 1;
    }
    match cmd.as_deref() {
        Some("encode") => {
            match dbfeats {
                Some(fp) => {
                    let feats = db::Feats::load(&fp);
                    eprintln!("features: {} ids, {} sets", feats.id2lab.len(), feats.sets.len());
                    let mut w = writer::Writer::new(opt, std::path::Path::new(&tmpdir))
                        .unwrap_or_else(|e| die(&e.to_string()));
                    w.set_feathash(feats.hash);
                    db::encode(&pos, &feats, &mut w);
                    let stdout = std::io::stdout();
                    let mut o = BufWriter::new(stdout.lock());
                    let st = w.close(&mut o).unwrap_or_else(|m| die(&m));
                    o.flush().unwrap_or_else(|e| die(&e.to_string()));
                    report(&st);
                }
                None => {
                    let a = load_aln(&alnpaf, &alnref, &alnnames);
                    encode(opt, std::path::Path::new(&tmpdir), a)
                }
            }
        }
        Some("smooth") => {
            if pos.len() < 3 { die("bad command or missing arguments"); }
            let n = |s: &str, w: &str| -> u64 { s.parse().unwrap_or_else(|_| die(w)) };
            let (seq, lo, hi) = match pos.len() {
                3 => (None, 0, 0),
                4 => (Some(pos[3].as_bytes()), 0, 0),
                6 => (Some(pos[3].as_bytes()),
                      n(&pos[4], "smooth: start is not a number"),
                      n(&pos[5], "smooth: end is not a number")),
                _ => die("bad command or missing arguments"),
            };
            let stdout = std::io::stdout();
            let mut o = BufWriter::new(stdout.lock());
            smooth::run(&pos[0], &pos[1], pos[2].as_bytes(), seq, lo, hi, v1, maxgap, halo, &ro, &mut o);
            o.flush().unwrap_or_else(|e| die(&e.to_string()));
        }
        Some("aln-pack") => {
            if pos.len() != 3 { die("usage: karyozip aln-pack P.paf REF.bed OUT"); }
            let mut a = aln::Aln::default();
            a.load_ref(&pos[1]);
            let rh = aln::refhash(&a);
            a.load_paf(&pos[0]);
            aln::pack(&a, &pos[2], rh);
        }
        Some("index") => {
            let f = pos.first().unwrap_or_else(|| die("bad command or missing arguments"));
            reader::Reader::open_with(f, &ro).build_index(f);
        }
        Some("region") => {
            if pos.len() < 2 { die("bad command or missing arguments"); }
            let n = |s: &str, w: &str| -> u64 {
                s.parse().unwrap_or_else(|_| die(w)) };
            let (lo, hi, halo) = match pos.len() {
                2 => (0, 0, 0),
                3 => { let p = n(&pos[2], "region: start is not a number"); (p, p + 1, 0) }
                4 | 5 => {
                    let a = n(&pos[2], "region: start is not a number");
                    let b = n(&pos[3], "region: end is not a number");
                    if a >= b { die("region: give START END with END past START"); }
                    (a, b, if pos.len() == 5 { n(&pos[4], "region: halo is not a number") } else { 0 })
                }
                _ => die("bad command or missing arguments"),
            };
            let mut r = reader::Reader::open_with(&pos[0], &ro);
            let stdout = std::io::stdout();
            let mut o = BufWriter::new(stdout.lock());
            r.region(pos[1].as_bytes(), lo, hi, halo, &mut o);
            o.flush().unwrap_or_else(|e| die(&e.to_string()));
        }
        Some("decode") => {
            if pos.len() > 1 { die("decode takes one track"); }
            let f = pos.first().unwrap_or_else(|| die("bad command or missing arguments"));
            ro.aln = load_aln(&alnpaf, &alnref, &alnnames);
            let mut r = reader::Reader::open_with(f, &ro);
            match (&dbfeats, &outdir) {
                (Some(fp), Some(d)) => {
                    let feats = db::Feats::load(fp);
                    eprintln!("features: {} ids, {} sets", feats.id2lab.len(), feats.sets.len());
                    if r.feathash() != feats.hash {
                        die("this features.tsv is not the one the track was encoded with");
                    }
                    db::decode(&mut r, &feats, d);
                }
                (Some(_), None) => die("decode --db needs --out DIR for the per-set BEDs"),
                _ => {
                    let stdout = std::io::stdout();
                    let mut o = BufWriter::new(stdout.lock());
                    r.decode(&mut o);
                    o.flush().unwrap_or_else(|e| die(&e.to_string()));
                }
            }
        }
        Some(c) => die(&format!("bad command {c}")),
        None => {
            eprintln!("usage: karyozip encode [--block N|all] [--level N] [--names keep] < in.bed > out.kst");
            std::process::exit(1);
        }
    }
}

/// --sort=MEM: digits and an optional K, M, G or T. The C spells the
/// refusals this way, so these do too.
fn mem_arg(p: &str) -> u64 {
    if p.is_empty() { die("--sort=MEM must be positive"); }
    let b = p.as_bytes();
    let mut val: u64 = 0;
    let mut i = 0;
    while i < b.len() && b[i].is_ascii_digit() {
        if val > (u64::MAX - (b[i] - b'0') as u64) / 10 {
            die("--sort=MEM is larger than the address space");
        }
        val = val * 10 + (b[i] - b'0') as u64;
        i += 1;
    }
    if i == 0 || val == 0 { die("--sort=MEM must be positive"); }
    let mul: u64 = if i == b.len() { 1 } else {
        if i + 1 != b.len() { die("--sort=MEM: nothing may follow the unit"); }
        match b[i] {
            b'k' | b'K' => 1 << 10, b'm' | b'M' => 1 << 20,
            b'g' | b'G' => 1 << 30, b't' | b'T' => 1 << 40,
            _ => die("--sort=MEM: use a plain number of bytes or a K, M, G or T suffix"),
        }
    };
    if val > u64::MAX / mul { die("--sort=MEM is larger than the address space"); }
    let m = val * mul;
    if m < (1 << 20) { die("--sort=MEM must be at least 1M"); }
    m
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

/// --aln needs --ref, and the two are read once, before any row.
fn load_aln(paf: &Option<String>, refbed: &Option<String>,
            ordfile: &Option<String>) -> Option<aln::Aln> {
    if paf.is_none() && ordfile.is_some() { die("--aln-names does nothing without --aln"); }
    match (paf, refbed) {
        (None, None) => None,
        (Some(_), None) => die("encode --aln needs --ref"),
        (None, Some(_)) => die("--ref does nothing without --aln"),
        (Some(p), Some(r)) => {
            let mut a = aln::Aln::default();
            a.load_ref(r);
            if let Some(o) = ordfile { a.load_ordinals(o); }
            a.load_paf(p);
            Some(a)
        }
    }
}

fn encode(opt: writer::Options, tmpdir: &std::path::Path, a: Option<aln::Aln>) {
    let mut w = writer::Writer::new(opt, tmpdir)
        .unwrap_or_else(|e| die(&format!("cannot create a temp file under {}: {e}", tmpdir.display())));
    let aln_on = a.is_some();
    if let Some(x) = a { w.set_aln(x); }
    // A slab at a time, like the C. Reading the whole input first would cost
    // the input's size in memory, 2.3 GB on the benchmark track.
    const SLAB: usize = 1 << 24;
    let mut slab = vec![0u8; SLAB];
    let mut have = 0usize;
    let nthr = rayon::current_num_threads().min(32);
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
        // Split the lines on threads, append in order. The append is
        // serial because the streams come out in row order.
        let rows = parse_slab(&slab[..end], nthr);
        for r in &rows {
            w.row(&slab[r.nm.0..r.nm.1], r.en - r.st,
                  &slab[r.lab.0..r.lab.1], r.st)
                .unwrap_or_else(|m| die(&m));
        }
        slab.copy_within(end..have, 0);
        have -= end;
        if eof && have == 0 { break; }
    }

    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    let hits = w.aln_hit();
    let rows_before = 0u64;
    let _ = rows_before;
    let st = w.close(&mut out).unwrap_or_else(|m| die(&m));
    out.flush().unwrap_or_else(|e| die(&e.to_string()));
    if aln_on {
        eprintln!("aln: {hits} of {} rows matched the reference ({:.1}%)",
                  st.rows, if st.rows > 0 { 100.0 * hits as f64 / st.rows as f64 } else { 0.0 });
    }
    report(&st);
}

fn report(st: &writer::Stats) {
    eprintln!("encoded rows={} reads={} blocks={} names={}B{} index={}B offs={}B lens={}B labs={}B raw={}B",
              st.rows, st.reads, st.blocks, st.names,
              if st.ordinal { " (ordinal, not stored)" }
              else if st.tokused { " (tokenised)" }
              else if st.intused { " (int deltas)" } else { "" },
              st.index, st.offs, st.lens, st.labs, st.raw);
}

/// One parsed row: byte ranges into the slab, and the coordinates.
struct Row { nm: (usize, usize), lab: (usize, usize), st: u64, en: u64 }

/// Splits `s` at line boundaries into `t` chunks, parses each on its own
/// thread, and returns the rows in order. A worker records the first thing
/// it refuses; main reports the earliest chunk's.
fn parse_slab(s: &[u8], t: usize) -> Vec<Row> {
    use rayon::prelude::*;
    let mut cut = vec![0usize];
    for k in 1..t {
        let mut g = s.len() * k / t;
        while g < s.len() && s[g] != b'\n' { g += 1; }
        cut.push(if g < s.len() { g + 1 } else { s.len() });
    }
    cut.push(s.len());
    for k in 1..cut.len() { if cut[k] < cut[k - 1] { cut[k] = cut[k - 1]; } }
    let parts: Vec<Result<Vec<Row>, &'static str>> = (0..cut.len() - 1)
        .into_par_iter()
        .map(|k| parse_chunk(s, cut[k], cut[k + 1]))
        .collect();
    let mut out = Vec::new();
    for p in parts {
        match p { Ok(mut v) => out.append(&mut v), Err(m) => die(m) }
    }
    out
}

fn parse_chunk(s: &[u8], from: usize, to: usize) -> Result<Vec<Row>, &'static str> {
    let mut out = Vec::new();
    let mut p = from;
    while p < to {
        let nl = match s[p..to].iter().position(|&c| c == b'\n') {
            Some(i) => p + i,
            None => break,
        };
        let line = &s[p..nl];
        if line.contains(&0) { return Err("a NUL byte inside a row"); }
        let t1 = line.iter().position(|&c| c == b'\t').ok_or("bad line")?;
        let t2 = line[t1 + 1..].iter().position(|&c| c == b'\t').ok_or("bad line")? + t1 + 1;
        let t3 = line[t2 + 1..].iter().position(|&c| c == b'\t').ok_or("bad line")? + t2 + 1;
        if t1 + 1 == t2 || t2 + 1 == t3 { return Err("empty start or end column"); }
        let st = bedcoord_r(&line[t1 + 1..t2], "start is not a number")?;
        let en = bedcoord_r(&line[t2 + 1..t3], "end is not a number")?;
        if en < st { return Err("end is before start"); }
        if en == st { return Err("end equals start: a run must cover at least one base"); }
        if t1 == 0 { return Err("empty sequence name"); }
        out.push(Row { nm: (p, p + t1), lab: (p + t3 + 1, nl), st, en });
        p = nl + 1;
    }
    Ok(out)
}

fn bedcoord_r(b: &[u8], what: &'static str) -> Result<u64, &'static str> {
    let mut v: u64 = 0;
    for &c in b {
        if !c.is_ascii_digit() { return Err(what); }
        if v > (u64::MAX - (c - b'0') as u64) / 10 {
            return Err("a coordinate does not fit in 64 bits");
        }
        v = v * 10 + (c - b'0') as u64;
    }
    Ok(v)
}

pub fn bedcoord(b: &[u8], what: &str) -> u64 {
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
