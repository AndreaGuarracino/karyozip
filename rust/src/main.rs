mod cli;

use karyozip::{aln, bail as die, bedcoord_r, db, ffi, kst, prof, reader, smooth, writer};

use cli::{Cli, Cmd};
use std::io::{BufWriter, Read, Write};


/// The C's wording for a number that is not one.
fn num(v: &Option<String>, lo: u64, hi: u64, what: &'static str) -> Option<u64> {
    let s = v.as_ref()?;
    let n: u64 = s.parse().unwrap_or_else(|_| die(what));
    if n < lo || n > hi { die(what); }
    Some(n)
}

/// A positional the command cannot run without.
fn need(v: &Option<String>) -> &str {
    v.as_deref().unwrap_or_else(|| die("bad command or missing arguments"))
}

fn main() {
    ffi::ignore_sigxfsz();
    for fd in 0..3 {
        if std::fs::read_link(format!("/proc/self/fd/{fd}")).is_err() {
            die("stdin, stdout or stderr is closed");
        }
    }
    let c = Cli::get();
    // --threads, else OMP_NUM_THREADS, else one per core capped at 8, which
    // is what the C does. The pool is always explicit, so RAYON_NUM_THREADS
    // cannot reach it and the two binaries read the one same variable.
    const THRCAP: usize = 8;
    let n = match num(&c.threads, 1, 1024, "--threads is a number of 1 to 1024") {
        Some(t) => t as usize,
        None => match std::env::var("OMP_NUM_THREADS").ok().and_then(|s| s.parse().ok()) {
            Some(m @ 1..) => usize::min(m, 1024),
            _ => std::thread::available_parallelism().map_or(1, |p| p.get()).min(THRCAP),
        },
    };
    rayon::ThreadPoolBuilder::new().num_threads(n).build_global().ok();
    let tmpdir = c.tmp_dir.clone()
        .or_else(|| std::env::var("TMPDIR").ok())
        .unwrap_or_else(|| "/tmp".into());

    let mut opt = writer::Options::default();
    if let Some(b) = &c.block {
        opt.block = if b == "all" { u64::MAX } else {
            b.parse().unwrap_or_else(|_| die("--block is a number of runs, 1 or more, or \"all\""))
        };
        if opt.block == 0 { die("--block is a number of runs, 1 or more, or \"all\""); }
    }
    if let Some(l) = num(&c.level, 1, 22, "--level is a zstd level, 1 to 22") {
        opt.level = l as i32;
        opt.levelset = true;
    }
    opt.slow = c.slow;
    if let Some(m) = &c.sort {
        opt.sort = true;
        // clap cannot tell "--sort" from "--sort=", and the C refuses the
        // second, so the empty size is looked for in argv itself.
        if m.is_empty() && std::env::args().any(|a| a == "--sort=") { mem_arg(""); }
        if !m.is_empty() { opt.sortmem = mem_arg(m); }
    }
    if let Some(t) = &c.tail {
        opt.tail = Some(match t.as_str() { "0" => false, "1" => true, _ => die("--tail is 0 or 1") });
    }
    if let Some(k) = num(&c.tok_chunk, 1, 0x7fff_ffff, "--tok-chunk needs a number") {
        opt.tok_chunk = Some((k as usize).min(9_000_000));
    }
    opt.id_codec = c.id_codec.as_deref().map(|v| parse_codec(v, "--id-codec"));
    opt.len_codec = c.len_codec.as_deref().map(|v| parse_codec(v, "--len-codec"));

    let Some(cmd) = c.cmd.as_ref() else { die(cli::USAGE) };
    prof::start(c.profile, n);
    // The C tests each predicate over both codecs before the next one.
    let ks = [opt.id_codec, opt.len_codec];
    let any = |f: fn(i32) -> bool| ks.iter().flatten().any(|&k| f(k));
    if any(ffi::badcodec) { die("--id-codec and --len-codec take a codec id, not any number"); }
    if any(ffi::nosz) { die("a NOSZ codec (0x10) stores no size and is not accepted"); }
    if any(ffi::noext) { die("arith X_EXT (0x04) needs libbz2 and is not accepted"); }
    c.check_scope().unwrap_or_else(|m| die(&m));
    let mut ro = reader::ReadOpts { raw: c.raw, db: c.db.is_some(), ..Default::default() };
    // The C keeps a mode and a file separately, so a repeat of one does not
    // clear the other, and drop plus drop:FILE is a refusal.
    for n in &c.names {
        match n.as_str() {
            "keep" => {
                if opt.drop_names { die("--names keep and --names drop are exclusive"); }
                opt.keep_names = true;
            }
            "drop" => {
                if opt.keep_names { die("--names keep and --names drop are exclusive"); }
                opt.drop_names = true;
            }
            _ if n.starts_with("drop:") => {
                if n.len() == 5 { die("--names drop:FILE needs a file"); }
                opt.names_out = Some(n[5..].to_string());
            }
            "ordinal" => ro.ordnames = true,
            _ => ro.names_file = Some(n.clone()),
        }
    }

    match cmd {
        Cmd::Encode { files } => {
            if c.db.is_none() && !files.is_empty() {
                die("without --db, encode reads one BED on stdin");
            }
            if opt.drop_names && opt.names_out.is_some() {
                die("--names drop and --names drop:FILE are exclusive");
            }
            if c.aln.is_some() && c.refbed.is_none() { die("--aln needs --ref"); }
            if c.aln.is_some() && opt.drop_names && c.aln_names.is_none() {
                die("--aln with --names drop needs --aln-names: without the names \
                     decode cannot reach the alignment");
            }
            if c.aln.is_some() && opt.sort && opt.drop_names {
                die("--aln with --sort and --names drop cannot be decoded: sorting \
                     reorders the records, so the ordinal name map is invalid");
            }
            if c.db.is_some() && files.is_empty() { die("--db: give one BED per feature set"); }
            let a = load_aln(&c.aln, &c.refbed, &c.aln_names, "");
            match &c.db {
                Some(fp) => {
                    let feats = db::Feats::load(fp);
                    eprintln!("features: {} ids, {} sets", feats.id2lab.len(), feats.sets.len());
                    let mut w = writer::Writer::new(opt, std::path::Path::new(&tmpdir))
                        .unwrap_or_else(|e| die(&e.to_string()));
                    w.set_feathash(feats.hash);
                    let aln_on = a.is_some();
                    if let Some(x) = a { w.set_aln(x); }
                    let files: Vec<String> = if files.is_empty() { vec!["-".into()] } else { files.clone() };
                    db::encode(&files, &feats, &mut w);
                    finish(w, aln_on);
                }
                None => {
                    if !files.is_empty() { die("without --db, encode reads one BED on stdin"); }
                    encode(opt, std::path::Path::new(&tmpdir), a)
                }
            }
        }
        Cmd::Decode { track, rest } => {
            if !rest.is_empty() { die("decode takes one track"); }
            let track = need(track);
            ro.aln = load_aln(&c.aln, &c.refbed, &c.aln_names, "decode ");
            // The C refuses the flag pair before it opens the track.
            if c.db.is_some() && c.out.is_none() { die("decode --db needs --out DIR for the per-set BEDs"); }
            let mut r = reader::Reader::open_with(track, &ro);
            match (&c.db, &c.out) {
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
                    let so = std::io::stdout();
                    let mut o = BufWriter::new(so.lock());
                    r.decode(&mut o);
                    o.flush().unwrap_or_else(|e| die(&e.to_string()));
                }
            }
        }
        Cmd::Region { track, seq, a, b, hl, rest: _ } => {
            let (track, seq) = (need(track), need(seq));
            let p = |v: &Option<String>, w: &'static str| -> u64 {
                v.as_ref().unwrap().parse().unwrap_or_else(|_| die(w))
            };
            let (lo, hi, hlo) = match (a, b) {
                (None, _) => (0, 0, 0),
                (Some(_), None) => { let x = p(a, "region: start is not a number"); (x, x + 1, 0) }
                (Some(_), Some(_)) => {
                    let x = p(a, "region: start is not a number");
                    let y = p(b, "region: end is not a number");
                    if x >= y { die("region: give START END with END past START"); }
                    (x, y, hl.as_ref().map_or(0, |_| p(hl, "region: halo is not a number")))
                }
            };
            ro.aln = load_aln(&c.aln, &c.refbed, &c.aln_names, "");
            let mut r = reader::Reader::open_with(track, &ro);
            let so = std::io::stdout();
            let mut o = BufWriter::new(so.lock());
            r.region(seq.as_bytes(), lo, hi, hlo, &mut o);
            o.flush().unwrap_or_else(|e| die(&e.to_string()));
        }
        Cmd::Smooth { track, hier, fs, seq, start, end, rest } => {
            let (track, hier, fs) = (need(track), need(hier), need(fs));
            if !rest.is_empty() { die("smooth: give SEQ, or SEQ START END"); }
            let maxgap = num(&c.max_gap, 0, u64::MAX, "--max-gap needs a number").unwrap_or(1000);
            let halo = num(&c.halo, 0, u64::MAX, "--halo needs a number").unwrap_or(200);
            let (lo, hi) = match (start, end) {
                (Some(a), Some(b)) => (
                    a.parse().unwrap_or_else(|_| die("smooth: start is not a number")),
                    b.parse().unwrap_or_else(|_| die("smooth: end is not a number"))),
                (None, None) => (0, 0),
                _ => die("smooth: give SEQ, or SEQ START END"),
            };
            if start.is_some() && lo >= hi { die("smooth: end must be past start"); }
            let so = std::io::stdout();
            let mut o = BufWriter::new(so.lock());
            smooth::run(track, hier, fs.as_bytes(), seq.as_ref().map(|s| s.as_bytes()),
                        lo, hi, c.v1, maxgap, halo, &ro, &mut o);
            o.flush().unwrap_or_else(|e| die(&e.to_string()));
        }
        Cmd::Index { track, rest: _ } => {
            let track = need(track);
            reader::Reader::open_with(track, &ro).build_index(track)
        }
        Cmd::Codecs { file } => probe(need(file)),
        Cmd::AlnPack { paf, rbed, dest, rest } => {
            let (paf, rbed, dest) = match (paf, rbed, dest) {
                (Some(p), Some(r), Some(o)) if rest.is_empty() => (p, r, o),
                _ => die("usage: karyozip aln-pack P.paf REF.bed OUT"),
            };
            let mut a = aln::Aln::default();
            a.load_ref(rbed);
            let rh = aln::refhash(&a);
            a.load_paf(paf);
            aln::pack(&a, dest, rh);
        }
    }
    prof::dump();
}

/// Writes the file and reports, for both encode paths.
fn finish(w: writer::Writer, aln_on: bool) {
    let hits = w.aln_hit();
    let so = std::io::stdout();
    let mut o = BufWriter::new(so.lock());
    let st = w.close(&mut o).unwrap_or_else(|m| die(&m));
    o.flush().unwrap_or_else(|e| die(&e.to_string()));
    if aln_on {
        eprintln!("aln: {hits} of {} rows matched the reference ({:.1}%)",
                  st.rows, if st.rows > 0 { 100.0 * hits as f64 / st.rows as f64 } else { 0.0 });
    }
    report(&st);
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

fn parse_codec(s: &str, flag: &str) -> i32 {
    let v = if let Some(h) = s.strip_prefix("0x") {
        i64::from_str_radix(h, 16)
    } else { s.parse::<i64>() };
    match v {
        Ok(v) if (0..=0x7fff_ffff).contains(&v) => v as i32,
        _ => die(&format!("{flag} needs a codec id")),
    }
}

/// --aln needs --ref, and the two are read once, before any row.
fn load_aln(paf: &Option<String>, refbed: &Option<String>,
            ordfile: &Option<String>, who: &str) -> Option<aln::Aln> {
    if paf.is_none() && ordfile.is_some() { die("--aln-names does nothing without --aln"); }
    match (paf, refbed) {
        (None, None) => None,
        (Some(_), None) => die(&format!("{who}--aln needs --ref")),
        (None, Some(_)) => die("--ref does nothing without --aln"),
        (Some(p), Some(r)) => {
            // The C loads the ordinal map first, and the diagnostics are
            // compared line for line.
            let mut a = aln::Aln::default();
            if let Some(o) = ordfile { a.load_ordinals(o); }
            a.load_ref(r);
            a.load_aln(p);
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
    let mut stdin = karyozip::unzip(std::io::stdin());   // plain or .gz
    let mut eof = false;
    let _p = prof::span(0);
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
        let _s = prof::span(5);
        let rows = parse_slab(&slab[..end], nthr);
        drop(_s);
        let _s = prof::span(6);
        for r in rows.iter().flatten() {
            w.push(&slab[r.nm.0 as usize..r.nm.1 as usize], r.en - r.st,
                  &slab[r.lab.0 as usize..r.lab.1 as usize], r.st)
                .unwrap_or_else(|m| die(&m));
        }
        drop(_s);
        slab.copy_within(end..have, 0);
        have -= end;
        if eof && have == 0 { break; }
    }

    drop(_p);
    finish(w, aln_on);
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
/// Offsets into the slab, which is 16 MB, so u32 is enough and the row is
/// 32 bytes instead of 48.
struct Row { nm: (u32, u32), lab: (u32, u32), st: u64, en: u64 }

/// Splits `s` at line boundaries into `t` chunks, parses each on its own
/// thread, and returns the rows in order. A worker records the first thing
/// it refuses; main reports the earliest chunk's.
fn parse_slab(s: &[u8], t: usize) -> Vec<Vec<Row>> {
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
    parts.into_iter().map(|p| p.unwrap_or_else(|m| die(m))).collect()
}

fn parse_chunk(s: &[u8], from: usize, to: usize) -> Result<Vec<Row>, &'static str> {
    let mut out: Vec<Row> = Vec::new();
    let mut p = from;
    while p < to {
        let Some(nl) = memchr::memchr(b'\n', &s[p..to]).map(|i| p + i) else { break };
        let line = &s[p..nl];
        if memchr::memchr(0, line).is_some() { return Err("a NUL byte inside a row"); }
        let t1 = memchr::memchr(b'\t', line).ok_or("bad line")?;
        let t2 = memchr::memchr(b'\t', &line[t1 + 1..]).ok_or("bad line")? + t1 + 1;
        let t3 = memchr::memchr(b'\t', &line[t2 + 1..]).ok_or("bad line")? + t2 + 1;
        if t1 + 1 == t2 || t2 + 1 == t3 { return Err("empty start or end column"); }
        let st = bedcoord_r(&line[t1 + 1..t2], "start is not a number")?;
        let en = bedcoord_r(&line[t2 + 1..t3], "end is not a number")?;
        if en < st { return Err("end is before start"); }
        if en == st { return Err("end equals start: a run must cover at least one base"); }
        if t1 == 0 { return Err("empty sequence name"); }
        out.push(Row { nm: (p as u32, (p + t1) as u32),
                       lab: ((p + t3 + 1) as u32, nl as u32), st, en });
        p = nl + 1;
    }
    Ok(out)
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
