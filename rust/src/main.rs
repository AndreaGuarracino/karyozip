mod cli;

use karyozip::{aln, bail as die, bedcoord_r, db, ffi, kst, prof, reader, smooth, writer};

use cli::{Cli, Cmd};
use std::io::{BufWriter, Read, Write};


/// The C's wording for a number that is not one.
/// A count: digits only, at most 18 of them, in range, as the C's int_arg.
fn count(s: &str, lo: u64, hi: u64, what: &'static str) -> u64 {
    if s.len() > 18 { die(what); }
    let n = karyozip::bedcoord(s.as_bytes(), what);
    if n < lo || n > hi { die(what); }
    n
}
fn num(v: &Option<String>, lo: u64, hi: u64, what: &'static str) -> Option<u64> {   // a coordinate-sized number, as num_arg
    v.as_ref().map(|s| { let n = karyozip::bedcoord(s.as_bytes(), what); if n < lo || n > hi { die(what); } n })
}
fn cnt(v: &Option<String>, lo: u64, hi: u64, what: &'static str) -> Option<u64> { v.as_ref().map(|s| count(s, lo, hi, what)) }

/// A positional the command cannot run without.
fn need(v: &Option<String>) -> &str {
    v.as_deref().unwrap_or_else(|| die("bad command or missing arguments"))
}

fn main() {
    ffi::ignore_sigxfsz();
    if !ffi::STDFDS_OPEN.load(std::sync::atomic::Ordering::Relaxed) { die("stdin, stdout or stderr is closed"); }
    let c = Cli::get();
    // --threads, else OMP_NUM_THREADS, else one per core capped at 8, which
    // is what the C does. The pool is always explicit, so RAYON_NUM_THREADS
    // cannot reach it and the two binaries read the one same variable.
    const THRCAP: usize = 8;
    let n = match cnt(&c.threads, 1, 1024, "--threads is a number of 1 to 1024") {
        Some(t) => t as usize,
        None => match std::env::var("OMP_NUM_THREADS").ok().and_then(|s| s.parse().ok()) {
            Some(m @ 1..) => usize::min(m, 1024),
            _ => std::thread::available_parallelism().map_or(1, |p| p.get()).min(THRCAP),
        },
    };
    ffi::pin_mmap_threshold();
    rayon::ThreadPoolBuilder::new().num_threads(n).build_global().ok();
    let tmpdir = c.tmp_dir.clone()
        .or_else(|| std::env::var("TMPDIR").ok())
        .unwrap_or_else(|| "/tmp".into());

    let mut opt = writer::Options::default();
    if let Some(b) = &c.block {
        opt.block = if b == "all" { u64::MAX } else {
            count(b, 1, i64::MAX as u64, "--block is a number of runs, 1 or more, or \"all\"")
        };
    }
    let effort = c.effort.as_deref().unwrap_or("medium");
    if !["low", "medium", "high"].contains(&effort) { die("--effort is low, medium or high"); }
    opt.slow = effort != "low";
    opt.whole = effort == "high";
    if effort == "low" && c.block.is_none() { opt.block = 262144; }   // -e low: blocks and an index, for queries
    if let Some(m) = &c.sort {
        opt.sort = true;
        if !m.is_empty() { opt.sortmem = mem_arg(m); }
    }

    let Some(cmd) = c.cmd.as_ref() else { die(cli::USAGE) };
    prof::start(c.profile, n);
    let mut ro = reader::ReadOpts { raw: c.raw, db: c.db.is_some(), ..Default::default() };
    // The C keeps a mode and a file separately, so a repeat of one does not
    // clear the other, and drop plus drop:FILE is a refusal.
    for n in &c.names {
        match n.as_str() {
            "drop" => opt.drop_names = true,
            _ if n.starts_with("drop:") => {
                if n.len() == 5 { die("--names drop:FILE needs a file"); }
                opt.names_out = Some(n[5..].to_string());
            }
            "ordinal" => ro.ordnames = true,
            _ => ro.names_file = Some(n.clone()),
        }
    }

    let name = c.name();
    let pack = name == "aln-pack";
    if c.input.len() > 1 && !(name == "encode" && c.db.is_some()) { die("give one -i"); }
    if let (Some(o), false) = (&c.output, pack) { ffi::redirect(o, 1); }
    let track = || c.input.first().map(|s| s.as_str())
        .unwrap_or_else(|| die(&format!("{name} needs -i TRACK")));
    match cmd {
        Cmd::Encode => {
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
            if c.db.is_some() && c.input.is_empty() { die("--db: give one -i BED per feature set"); }
            let a = load_aln(&c.aln, &c.refbed, &c.aln_names, "", aln::Lazy::No);
            match &c.db {
                Some(fp) => {
                    let feats = load_feats(fp, c.hier.as_deref());
                    let mut w = writer::Writer::new(opt, std::path::Path::new(&tmpdir))
                        .unwrap_or_else(|e| die(&e.to_string()));
                    let aln_on = a.is_some();
                    if let Some(x) = a { w.set_aln(x); }
                    let used = db::encode(&c.input, &feats, &mut w);
                    w.set_feathash(if used { feats.hierhash } else { feats.hash }, used);
                    finish(w, aln_on);
                }
                None => {
                    if c.hier.is_some() { die("--hier does nothing without --db"); }
                    if let Some(i) = c.input.first() { ffi::redirect(i, 0); }
                    encode(opt, std::path::Path::new(&tmpdir), a)
                }
            }
        }
        Cmd::Decode if c.region.is_none() => {
            let track = track();
            if c.halo.is_some() { die("--halo needs -r NAME:START-END"); }
            let feats = informed(&c, &mut ro, "decode ", None, aln::Lazy::Contig);
            let mut r = reader::Reader::open_with(track, &ro);
            if let Some(f) = &feats { f.check_hash(&r); }
            match &feats {
                Some(f) => db::decode(&mut r, f, c.out.as_deref(), &f.setcols(c.set.as_deref())),
                None => {
                    let so = std::io::stdout();
                    let mut o = BufWriter::new(so.lock());
                    r.decode(&mut o);
                    o.flush().unwrap_or_else(|e| die(&e.to_string()));
                }
            }
        }
        Cmd::Decode => {   // one sequence, or the runs overlapping a range of it
            let track = track();
            if c.out.is_some() { die("-r and --out-dir: give one"); }
            if c.raw { die("-r and --raw: give one"); }
            let (seq, lo, hi) = regparse(c.region.as_deref().unwrap());
            if c.halo.is_some() && hi == 0 { die("--halo needs -r NAME:START-END"); }
            let hlo = c.halo.as_ref().map_or(0, |h| karyozip::bedcoord(h.as_bytes(), "--halo is not a number"));
            let feats = informed(&c, &mut ro, "decode ", None, aln::Lazy::OneSeq);
            let mut r = reader::Reader::open_with(track, &ro);
            if let Some(f) = &feats { f.check_hash(&r); }
            let cols = feats.as_ref().map_or(Vec::new(), |f| f.setcols(c.set.as_deref()));
            let so = std::io::stdout();
            let mut o = BufWriter::new(so.lock());
            r.region(seq.as_bytes(), lo, hi, hlo, feats.as_ref(), &cols, &mut o);
            o.flush().unwrap_or_else(|e| die(&e.to_string()));
        }
        Cmd::Smooth => {
            let track = track();
            let (Some(hier), Some(fs)) = (c.hier.as_deref(), c.feature_set.as_deref())
                else { die("smooth needs --hier HIER and -f FEATURE_SET") };
            let maxgap = num(&c.max_gap, 0, u64::MAX, "--max-gap needs a number").unwrap_or(1000);
            let (seq, lo, hi) = c.region.as_deref().map_or(("", 0, 0), regparse);
            let feats = informed(&c, &mut ro, "smooth ", Some(hier), aln::Lazy::Contig);
            // not stdout().lock(): a lock guard cannot cross to the smooth
            // workers, and they are serialized by turn anyway
            let mut o = BufWriter::new(std::io::stdout());
            smooth::run(track, hier, fs.as_bytes(), c.region.as_ref().map(|_| seq.as_bytes()),
                        lo, hi, c.v1, maxgap, &ro, feats.as_ref(), &mut o);
            o.flush().unwrap_or_else(|e| die(&e.to_string()));
        }
        Cmd::Index => {
            let track = track();
            ro.index = true;
            reader::Reader::open_with(track, &ro).build_index(track)
        }
        Cmd::Codecs { file } => probe(need(file)),
        Cmd::AlnPack => {
            let (Some(paf), Some(rbed), Some(dest)) = (&c.aln, &c.refbed, &c.output)
                else { die("usage: karyozip aln-pack -a P.paf --ref REF.bed|REF.kst -o OUT") };
            let mut a = aln::Aln::default();
            a.load_ref(rbed, aln::Lazy::No);
            let rh = aln::refhash(&a);
            a.load_paf(paf);
            aln::pack(&a, dest, rh);
        }
    }
    prof::dump();
}

/// Writes the file and reports, for both encode paths.
fn finish(w: writer::Writer, aln_on: bool) {
    let so = std::io::stdout();
    let mut o = BufWriter::new(so.lock());
    let st = w.close(&mut o).unwrap_or_else(|m| die(&m));
    o.flush().unwrap_or_else(|e| die(&e.to_string()));
    if aln_on {
        eprintln!("aln: {} of {} rows matched the reference ({:.1}%)", st.alnhit,
                  st.rows, if st.rows > 0 { 100.0 * st.alnhit as f64 / st.rows as f64 } else { 0.0 });
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

/// -r NAME or NAME:START-END, the range BED's: 0-based, END excluded. A
/// name may hold ':' itself, so only a digits-digits tail is a range.
fn regparse(r: &str) -> (&str, u64, u64) {
    let digits = |x: &str| !x.is_empty() && x.bytes().all(|b| b.is_ascii_digit());
    let Some((nm, tail)) = r.rsplit_once(':') else { return (r, 0, 0) };
    let Some((a, b)) = tail.split_once('-') else { return (r, 0, 0) };
    if !digits(a) || !digits(b) { return (r, 0, 0); }
    let lo = karyozip::bedcoord(a.as_bytes(), "-r: START is not a number");
    if lo == u64::MAX { die("-r: START is not a number"); }
    let hi = karyozip::bedcoord(b.as_bytes(), "-r: END is not a number");
    if lo >= hi { die("-r: END must be past START"); }
    (nm, lo, hi)
}

/// --aln needs --ref, and the two are read once, before any row.
/// --db and --aln companions for a reader, in the C's order: features first.
/// smooth's positional hierarchy stands in for --hier under --db.
fn informed(c: &cli::Cli, ro: &mut reader::ReadOpts, who: &str, hier: Option<&str>,
            lazy: aln::Lazy) -> Option<db::Feats> {
    if c.set.is_some() && c.db.is_none() { die("--set needs --db"); }
    if c.hier.is_some() && c.db.is_none() && hier.is_none() { die("--hier does nothing without --db"); }
    let feats = c.db.as_deref().map(|fp| load_feats(fp, c.hier.as_deref()));
    ro.aln = load_aln(&c.aln, &c.refbed, &c.aln_names, who, lazy).map(std::sync::Arc::new);
    feats
}

fn load_feats(fp: &str, hier: Option<&str>) -> db::Feats {
    let mut f = db::Feats::load(fp);
    eprintln!("features: {} ids, {} sets", f.nids, f.sets.len());
    if f.sets.len() <= 32 { eprintln!("features: integer-tuple index on"); }
    if let Some(h) = hier { f.load_hier(h); }
    f
}

fn load_aln(paf: &Option<String>, refbed: &Option<String>,
            ordfile: &Option<String>, who: &str, lazy: aln::Lazy) -> Option<aln::Aln> {
    match (paf, refbed) {
        (Some(_), None) => die(&format!("{who}--aln needs --ref")),
        (None, Some(_)) => die("--ref does nothing without --aln"),
        (None, None) if ordfile.is_some() => die("--aln-names does nothing without --aln"),
        (None, None) => None,
        (Some(p), Some(r)) => {
            // The C loads the ordinal map first, and the diagnostics are
            // compared line for line.
            let mut a = aln::Aln::default();
            if let Some(o) = ordfile { a.load_ordinals(o); }
            a.load_ref(r, lazy);
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
              if st.dropped { " (dropped, needs --names FILE)" }
              else if st.ordinal { " (ordinal, not stored)" }
              else if st.tokused { " (tokenised)" }
              else if st.intused { " (int deltas)" }
              else if st.fldused { " (field deltas)" } else { "" },
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
    // a row is at least 8 bytes: sized once, so no doubling copy and only
    // the rows used are touched, as the C's per-thread arrays
    let mut out: Vec<Row> = Vec::with_capacity((to - from) / 8 + 1);
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
    std::fs::File::open(path).and_then(|mut f| f.read_to_end(&mut src)).unwrap_or_else(|e| die(&format!("{path}: {e}")));
    let mut out = std::io::stdout().lock();
    for lvl in [1, 3, 17, 19, 22] {
        let z = ffi::zstd(&src, lvl);
        writeln!(out, "zstd:{lvl}\t{}\t{:016x}", z.len(), kst::fnv(kst::FNV0, &z)).unwrap();
    }
    for c in writer::candlist(true) {
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
