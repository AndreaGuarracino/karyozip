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

    /// rows per block, or "all" for one block
    #[arg(long, global = true, value_name = "N|all")]
    pub block: Option<String>,
    /// zstd level, 1 to 22
    #[arg(long, global = true, value_name = "N")]
    pub level: Option<String>,
    /// every codec over every block: smaller, slower
    #[arg(long, global = true)]
    pub slow: bool,
    /// order by label tuple, then lengths
    #[arg(long, global = true, num_args = 0..=1, require_equals = true,
          default_missing_value = "", value_name = "MEM")]
    pub sort: Option<String>,
    /// keep|drop|drop:FILE on encode, a file or "ordinal" on the readers
    #[arg(long, global = true, value_name = "WHAT")]
    pub names: Vec<String>,
    /// one BED per feature set, or one merged featureID track
    #[arg(long, global = true, value_name = "features.tsv")]
    pub db: Option<String>,
    /// with --db: label tuples outside features.tsv get an id from the hierarchy
    #[arg(long, global = true, value_name = "hierarchy.tsv")]
    pub hier: Option<String>,
    /// per-set BEDs go here
    #[arg(long, global = true, value_name = "DIR")]
    pub out: Option<String>,
    /// decode, region: one set of a --db track
    #[arg(long, global = true, value_name = "NAME")]
    pub set: Option<String>,
    /// do not store what the reference already gives
    #[arg(long, global = true, value_name = "P.paf")]
    pub aln: Option<String>,
    /// the reference track --aln deltas against
    #[arg(long = "ref", global = true, value_name = "R")]
    pub refbed: Option<String>,
    /// ordinal N is line N+1 of this file
    #[arg(long = "aln-names", global = true, value_name = "F")]
    pub aln_names: Option<String>,
    /// show what is stored, reversing no transform
    #[arg(long, global = true)]
    pub raw: bool,
    /// force the tail predictor off or on
    #[arg(long, global = true, value_name = "0|1")]
    pub tail: Option<String>,
    /// names per tokeniser chunk
    #[arg(long = "tok-chunk", global = true, value_name = "N")]
    pub tok_chunk: Option<String>,
    #[arg(long = "id-codec", global = true, value_name = "N")]
    pub id_codec: Option<String>,
    #[arg(long = "len-codec", global = true, value_name = "N")]
    pub len_codec: Option<String>,
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
    Encode { files: Vec<String> },
    /// give the track back
    Decode { track: Option<String>, rest: Vec<String> },
    /// the whole sequence, the run at POS, or an overlap
    Region { track: Option<String>, seq: Option<String>, a: Option<String>,
             b: Option<String>, hl: Option<String>, rest: Vec<String> },
    /// promote noisy runs through a hierarchy
    /// (clap shares one id between a global flag and a positional of the same name)
    Smooth { track: Option<String>, #[arg(value_name = "HIER")] hpath: Option<String>, fs: Option<String>,
             seq: Option<String>, start: Option<String>, end: Option<String>,
             rest: Vec<String> },
    /// append a name index to the track
    Index { track: Option<String>, rest: Vec<String> },
    /// pack a PAF to only what --aln reads
    #[command(name = "aln-pack")]
    AlnPack { paf: Option<String>, rbed: Option<String>, dest: Option<String>,
              rest: Vec<String> },
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
        if self.level.is_some() && !enc { bad("--level")?; }
        if self.slow && !enc { bad("--slow")?; }
        if self.sort.is_some() && !enc { bad("--sort")?; }
        if self.tail.is_some() && !enc { bad("--tail")?; }
        if self.tok_chunk.is_some() && !enc { bad("--tok-chunk")?; }
        if self.id_codec.is_some() && !enc { bad("--id-codec")?; }
        if self.len_codec.is_some() && !enc { bad("--len-codec")?; }
        if self.out.is_some() && !dec { bad("--out")?; }
        if self.set.is_some() && !(dec || c == "region") { bad("--set")?; }
        if self.raw && !dec { bad("--raw")?; }
        let inf = enc || dec || sm || c == "region";
        if self.db.is_some() && !inf { bad("--db")?; }
        if self.hier.is_some() && !inf { bad("--hier")?; }
        if self.aln.is_some() && !inf { bad("--aln")?; }
        if self.refbed.is_some() && !inf { bad("--ref")?; }
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
