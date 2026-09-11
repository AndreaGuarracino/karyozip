//! The command line, in clap. Every flag is global, so it works on either
//! side of the command name as the C does, and the scope check afterwards
//! gives the C's own refusal for a flag on the wrong command.

use clap::{Parser, Subcommand};

/// The C's own usage text, lifted from karyozip.c by build.rs.
pub const USAGE: &str = include_str!(concat!(env!("OUT_DIR"), "/usage.txt"));

/// What --version prints: the tool, the format, and the linked codecs.
fn version_text() -> &'static str {
    static V: std::sync::OnceLock<String> = std::sync::OnceLock::new();
    V.get_or_init(|| format!("{}, format KST00001\n  zstd       {}\n  libbz3     {}\n  htscodecs  {}",
            env!("CARGO_PKG_VERSION"), karyozip::ffi::zstd_version(),
            karyozip::ffi::bzip3_version(), karyozip::ffi::HTSCODECS_VERSION))
}

#[derive(Parser)]
#[command(name = "karyozip", version = version_text(), disable_help_flag = true,
          disable_help_subcommand = true, args_override_self = true,
          about = "Reads and writes the KaryoScope Track (KST) format")]
pub struct Cli {
    #[command(subcommand)]
    pub cmd: Option<Cmd>,

    /// encode: the BED, repeated under --db; readers: the track
    #[arg(long, short = 'i', global = true, value_name = "FILE")]
    pub input: Vec<String>,
    /// where the output goes, else stdout
    #[arg(long, short = 'o', global = true, value_name = "FILE")]
    pub output: Option<String>,
    /// region, smooth: the sequence
    #[arg(long, short = 's', global = true, value_name = "SEQ")]
    pub seq: Option<String>,
    /// region: the run covering this position
    #[arg(long, short = 'p', global = true, value_name = "POS")]
    pub pos: Option<String>,
    #[arg(long, global = true, value_name = "N")]
    pub start: Option<String>,
    #[arg(long, global = true, value_name = "N")]
    pub end: Option<String>,
    /// region: runs each side of --start --end
    #[arg(long, global = true, value_name = "N")]
    pub halo: Option<String>,
    /// smooth: the feature set
    #[arg(long = "feature-set", short = 'f', global = true, value_name = "FS")]
    pub feature_set: Option<String>,
    /// rows per block, or "all" for one block
    #[arg(long, short = 'b', global = true, value_name = "N|all")]
    pub block: Option<String>,
    /// low: quick, with blocks and an index; medium; high: the smallest
    #[arg(long, short = 'e', global = true, value_name = "low|medium|high")]
    pub effort: Option<String>,
    /// order by label tuple, then lengths
    #[arg(long, global = true, num_args = 0..=1, require_equals = true,
          default_missing_value = "", value_name = "MEM")]
    pub sort: Option<String>,
    /// keep|drop|drop:FILE on encode, a file or "ordinal" on the readers
    #[arg(long, short = 'n', global = true, value_name = "WHAT")]
    pub names: Vec<String>,
    /// one BED per feature set, or one merged featureID track
    #[arg(long, short = 'd', global = true, value_name = "features.tsv")]
    pub db: Option<String>,
    /// with --db: label tuples outside features.tsv get an id from the hierarchy
    #[arg(long, global = true, value_name = "hierarchy.tsv")]
    pub hier: Option<String>,
    /// per-set BEDs go here
    #[arg(long = "out-dir", global = true, value_name = "DIR")]
    pub out: Option<String>,
    /// decode, region: one set of a --db track
    #[arg(long, global = true, value_name = "NAME")]
    pub set: Option<String>,
    /// do not store what the reference already gives
    #[arg(long, short = 'a', global = true, value_name = "P.paf")]
    pub aln: Option<String>,
    /// the reference track --aln deltas against
    #[arg(long = "ref", short = 'r', global = true, value_name = "R")]
    pub refbed: Option<String>,
    /// ordinal N is line N+1 of this file
    #[arg(long = "aln-names", global = true, value_name = "F")]
    pub aln_names: Option<String>,
    /// show what is stored, reversing no transform
    #[arg(long, global = true)]
    pub raw: bool,
    /// smooth: the largest gap two runs may span
    #[arg(long = "max-gap", global = true, value_name = "N")]
    pub max_gap: Option<String>,
    /// smooth: v1 label names, dropping the _specific suffix
    #[arg(long, global = true)]
    pub v1: bool,
    /// threads, 1 to 1024
    #[arg(long, short = 't', global = true, value_name = "N")]
    pub threads: Option<String>,
    /// where temp files go, else TMPDIR, else /tmp
    #[arg(long = "tmp-dir", global = true, value_name = "DIR")]
    pub tmp_dir: Option<String>,
    /// phase times and peak memory on stderr
    #[arg(long, global = true)]
    pub profile: bool,
}

#[derive(Subcommand)]
pub enum Cmd {
    /// compress a BED-like track
    Encode { rest: Vec<String> },
    /// give the track back
    Decode { rest: Vec<String> },
    /// the whole sequence, the run at a position, or an overlap
    Region { rest: Vec<String> },
    /// promote noisy runs through a hierarchy
    Smooth { rest: Vec<String> },
    /// append a name index to the track
    Index { rest: Vec<String> },
    /// pack a PAF to only what --aln reads
    #[command(name = "aln-pack")]
    AlnPack { rest: Vec<String> },
    /// every candidate codec over one file, for tests/diff.sh
    #[command(hide = true)]
    Codecs { file: Option<String> },
}

impl Cli {
    /// clap, then the C's wording for what clap refuses. --help and -h win
    /// from either side of the command name, as they do in the C.
    pub fn get() -> Cli {
        use clap::error::{ContextKind, ContextValue, ErrorKind};
        if std::env::args().skip(1).any(|a| a == "--help" || a == "-h") {
            println!("{USAGE}");
            std::process::exit(0);
        }
        Cli::try_parse().unwrap_or_else(|e| {
            if e.kind() == ErrorKind::DisplayVersion { print!("{e}"); std::process::exit(0); }
            let arg = match e.get(ContextKind::InvalidArg) {
                Some(ContextValue::String(f)) => f.split_whitespace().next().unwrap_or(""),
                _ => "",
            };
            // The C echoes the spelling that was used, and -t is the one short.
            let arg = if arg == "--threads" && std::env::args().any(|a| a == "-t")
                { "-t" } else { arg };
            karyozip::bail(&match e.kind() {
                ErrorKind::UnknownArgument => format!("unknown flag {arg}"),
                ErrorKind::InvalidValue if !arg.is_empty() => format!("{arg} needs a value"),
                _ => "bad command or missing arguments".to_string(),
            })
        })
    }

    pub fn name(&self) -> &'static str {
        match self.cmd.as_ref().unwrap() {
            Cmd::Encode { .. } => "encode",
            Cmd::Decode { .. } => "decode",
            Cmd::Region { .. } => "region",
            Cmd::Smooth { .. } => "smooth",
            Cmd::Index { .. } => "index",
            Cmd::AlnPack { .. } => "aln-pack",
            Cmd::Codecs { .. } => "codecs",
        }
    }

    /// A flag on the wrong command is refused by name, as the C does.
    pub fn check_scope(&self) -> Result<(), String> {
        let c = self.name();
        let enc = c == "encode";
        let dec = c == "decode";
        let sm = c == "smooth";
        let rd = sm || dec || c == "region" || c == "index";
        let bad = |f: &str| -> Result<(), String> {
            Err(format!("{f} does not apply to {c}"))
        };
        if self.block.is_some() && !enc { bad("--block")?; }
        if self.effort.is_some() && !enc { bad("--effort")?; }
        if self.sort.is_some() && !enc { bad("--sort")?; }
        if self.out.is_some() && !dec { bad("--out-dir")?; }
        let (pack, reg) = (c == "aln-pack", c == "region");
        if !self.input.is_empty() && pack { bad("--input")?; }
        if self.output.is_some() && c == "index" { bad("--output")?; }
        if self.seq.is_some() && !(sm || reg) { bad("--seq")?; }
        if self.pos.is_some() && !reg { bad("--pos")?; }
        if self.start.is_some() && !(sm || reg) { bad("--start")?; }
        if self.end.is_some() && !(sm || reg) { bad("--end")?; }
        if self.halo.is_some() && !reg { bad("--halo")?; }
        if self.feature_set.is_some() && !sm { bad("--feature-set")?; }
        if self.set.is_some() && !(dec || c == "region") { bad("--set")?; }
        if self.raw && !dec { bad("--raw")?; }
        let inf = enc || dec || sm || c == "region";
        if self.db.is_some() && !inf { bad("--db")?; }
        if self.hier.is_some() && !inf { bad("--hier")?; }
        if self.aln.is_some() && !(inf || pack) { bad("--aln")?; }
        if self.refbed.is_some() && !(inf || pack) { bad("--ref")?; }
        if self.aln_names.is_some() && !inf { bad("--aln-names")?; }
        if self.max_gap.is_some() && !sm { bad("--max-gap")?; }
        if self.v1 && !sm { bad("--v1")?; }
        for n in &self.names {
            let enc_form = n == "keep" || n == "drop" || n.starts_with("drop:");
            if enc && !enc_form {
                return Err(if n == "ordinal" {
                    "encode takes --names keep, drop or drop:FILE".into()
                } else {
                    "encode takes --names keep, drop or drop:FILE; a name file is --aln-names".into()
                });
            }
            if !enc && enc_form { bad("--names")?; }
            if !enc && !rd { bad("--names")?; }
        }
        Ok(())
    }
}
