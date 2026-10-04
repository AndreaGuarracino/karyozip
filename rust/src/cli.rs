//! The command line, parsed as the C parses it: the same flag table, flags
//! on either side of the command, short forms rewritten to long ones, and
//! the C's own refusal for each mistake, in the same order.

use karyozip::bail as die;

/// The C's own usage text, lifted from karyozip.c by build.rs.
pub const USAGE: &str = include_str!(concat!(env!("OUT_DIR"), "/usage.txt"));

pub enum Cmd { Encode, Decode, Smooth, Index, AlnPack, Codecs { file: Option<String> } }

#[derive(Default)]
pub struct Cli {
    pub cmd: Option<Cmd>,
    pub input: Vec<String>,
    pub output: Option<String>,
    pub region: Option<String>,
    pub halo: Option<String>,
    pub feature_set: Option<String>,
    pub block: Option<String>,
    pub effort: Option<String>,
    pub sort: Option<String>,        // "" for a bare --sort
    pub names: Vec<String>,
    pub db: Option<String>,
    pub hier: Option<String>,
    pub out: Option<String>,
    pub set: Option<String>,
    pub aln: Option<String>,
    pub refbed: Option<String>,
    pub aln_names: Option<String>,
    pub raw: bool,
    pub max_gap: Option<String>,
    pub v1: bool,
    pub threads: Option<String>,
    pub tmp_dir: Option<String>,
    pub profile: bool,
}

/// Every flag and whether it takes a value, as the C's FLAGS.
const FLAGS: &[(&str, bool)] = &[
    ("--version", false), ("--threads", true), ("-t", true), ("--profile", false), ("--tmp-dir", true),
    ("--names", true), ("--max-gap", true), ("--v1", false), ("--aln-names", true), ("--sort", false),
    ("--effort", true), ("--block", true), ("--db", true), ("--aln", true), ("--ref", true),
    ("--out-dir", true), ("--set", true), ("--raw", false), ("--hier", true), ("--input", true),
    ("--output", true), ("--region", true), ("--halo", true), ("--feature-set", true),
];
const SHORT: &[(&str, &str)] = &[
    ("-i", "--input"), ("-o", "--output"), ("-r", "--region"), ("-f", "--feature-set"),
    ("-b", "--block"), ("-e", "--effort"), ("-n", "--names"), ("-d", "--db"), ("-a", "--aln"),
];

fn long(a: String) -> String {
    SHORT.iter().find(|(s, _)| *s == a).map_or(a, |(_, l)| l.to_string())
}

impl Cli {
    pub fn get() -> Cli {
        let mut av: Vec<String> = std::env::args().skip(1).collect();
        for a in &av {
            if a == "--help" || a == "-h" { println!("{USAGE}"); std::process::exit(0); }
            if a == "--version" {
                println!("karyozip {}, format KST00001\n  zstd       {}\n  libbz3     {}\n  htscodecs  {}\n  libdeflate {}",
                         env!("CARGO_PKG_VERSION"), karyozip::ffi::zstd_version(),
                         karyozip::ffi::bzip3_version(), karyozip::ffi::HTSCODECS_VERSION, karyozip::ffi::LIBDEFLATE_VERSION);
                std::process::exit(0);
            }
        }
        if av.is_empty() { die(USAGE); }
        // the command: the first word that is not a flag or its value
        let mut cmd = String::new();
        let mut i = 0;
        while i < av.len() {
            let a = long(av[i].clone());
            if !a.starts_with('-') || a.len() == 1 { cmd = a; break; }
            let nv = if a.starts_with("--sort=") { Some(false) }
                     else { FLAGS.iter().find(|(f, _)| *f == a).map(|f| f.1) };
            match nv { Some(v) => i += 1 + v as usize, None => die(&format!("unknown flag {a}")) }
        }
        av = av.into_iter().map(long).collect();
        let c = cmd.as_str();
        let (enc, dec, sm) = (c == "encode", c == "decode", c == "smooth");
        let rd = sm || dec || c == "index";
        let inf = enc || sm || dec;
        let pack = c == "aln-pack";
        let mut o = Cli::default();
        let mut i = 0;
        while i < av.len() {
            let a = av[i].clone();
            let v = av.get(i + 1).cloned();
            if FLAGS.iter().any(|(f, nv)| *f == a && *nv) && v.is_none() { die(&format!("{a} needs a value")); }
            let val = || v.clone().unwrap();
            let ok = match a.as_str() {
                "--threads" | "-t" => { crate::count(&val(), 1, 1024, "--threads is a number of 1 to 1024"); o.threads = Some(val()); true }
                "--profile" => { o.profile = true; true }
                "--tmp-dir" => { o.tmp_dir = Some(val()); true }
                "--names" => {
                    let n = val();
                    let ok = if n == "drop" || n.starts_with("drop:") {
                        if n == "drop:" { die("--names drop:FILE needs a file"); }
                        enc
                    } else if enc {
                        die(if n == "ordinal" { "encode takes --names drop or drop:FILE" }
                            else { "encode takes --names drop or drop:FILE; a name file is --aln-names" })
                    } else { rd };
                    o.names.push(n); ok
                }
                "--max-gap" => { crate::num(&v, 0, u64::MAX, "--max-gap needs a number"); o.max_gap = Some(val()); sm }
                "--v1" => { o.v1 = true; sm }
                "--aln-names" => { o.aln_names = Some(val()); inf }
                "--sort" => { o.sort = Some(String::new()); enc }
                s if s.starts_with("--sort=") => {
                    crate::mem_arg(&s[7..]);
                    o.sort = Some(s[7..].to_string()); enc
                }
                "--effort" => {
                    if !["low", "medium", "high"].contains(&val().as_str()) { die("--effort is low, medium or high"); }
                    o.effort = Some(val()); enc
                }
                "--block" => {
                    if val() != "all" { crate::count(&val(), 1, i64::MAX as u64, "--block is a number of runs, 1 or more, or \"all\""); }
                    o.block = Some(val()); enc
                }
                "--db" => { o.db = Some(val()); inf }
                "--hier" => { o.hier = Some(val()); inf }
                "--aln" => { o.aln = Some(val()); inf || pack }
                "--ref" => { o.refbed = Some(val()); inf || pack }
                "--out-dir" => { o.out = Some(val()); dec }
                "--input" => { o.input.push(val()); !pack }
                "--output" => { o.output = Some(val()); c != "index" }
                "--region" => { o.region = Some(val()); sm || dec }
                "--halo" => { o.halo = Some(val()); dec }
                "--feature-set" => { o.feature_set = Some(val()); sm }
                "--set" => { o.set = Some(val()); dec }
                "--raw" => { o.raw = true; dec }
                _ => { i += 1; continue; }
            };
            if !ok {
                if cmd.is_empty() { die(USAGE); }
                die(&format!("{a} does not apply to {cmd}"));
            }
            av.drain(i..i + if matches!(a.as_str(), "--profile" | "--v1" | "--raw") || a.starts_with("--sort") { 1 } else { 2 });
        }
        if cmd.is_empty() { die(USAGE); }
        o.cmd = Some(match c {
            "encode" => Cmd::Encode, "decode" => Cmd::Decode,
            "smooth" => Cmd::Smooth, "index" => Cmd::Index, "aln-pack" => Cmd::AlnPack,
            "codecs" => Cmd::Codecs { file: if av.len() > 1 { Some(av.remove(1)) } else { None } },
            _ => die("bad command or missing arguments"),
        });
        if let Some(a) = av.get(1) {   // every argument is a flag now
            die(&if a.starts_with('-') && a.len() > 1 { format!("unknown flag {a}") }
                 else { format!("unexpected argument {a}") });
        }
        o
    }

    pub fn name(&self) -> &'static str {
        match self.cmd.as_ref().unwrap() {
            Cmd::Encode => "encode", Cmd::Decode => "decode",
            Cmd::Smooth => "smooth", Cmd::Index => "index", Cmd::AlnPack => "aln-pack",
            Cmd::Codecs { .. } => "codecs",
        }
    }
}
