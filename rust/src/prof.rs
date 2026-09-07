//! --profile: the same phases, and the same closing line, as the C prints.
//! Times and memory are the machine's, so they are never compared byte for
//! byte with the C; only the shape of the output is the same.
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering::Relaxed};
use std::sync::OnceLock;
use std::time::Instant;

pub static ON: AtomicBool = AtomicBool::new(false);
static NS: [AtomicU64; 10] = [const { AtomicU64::new(0) }; 10];
static HW: [AtomicU64; 10] = [const { AtomicU64::new(0) }; 10];
static THR: AtomicU64 = AtomicU64::new(0);
static T0: OnceLock<Instant> = OnceLock::new();
static NAME: [&str; 9] = ["parse", "compress", "names", "index", "smooth",
                          "parse:lines", "parse:append", "write", "codec search"];

pub fn start(on: bool, threads: usize) {
    ON.store(on, Relaxed);
    THR.store(threads as u64, Relaxed);
    let _ = T0.set(Instant::now());
}

/// Charges its phase for as long as it lives.
pub struct Span(usize, Instant);
pub fn span(i: usize) -> Option<Span> {
    (ON.load(Relaxed) && i < NAME.len()).then(|| Span(i, Instant::now()))   // 9 is "not this phase"
}
/// The high-water mark is read from /proc at most every 10 ms, as the C.
static LAST: AtomicU64 = AtomicU64::new(0);
impl Drop for Span {
    fn drop(&mut self) {
        NS[self.0].fetch_add(self.1.elapsed().as_nanos() as u64, Relaxed);
        let now = T0.get().map_or(0, |t| t.elapsed().as_nanos() as u64);
        if now - LAST.load(Relaxed) > 10_000_000 {
            LAST.store(now, Relaxed);
            HW[self.0].store(field("/proc/self/status", "VmHWM:"), Relaxed);
        }
    }
}

/// One kB number from a "Key: N kB" line, or 0 when /proc is not there.
fn field(path: &str, key: &str) -> u64 {
    std::fs::read_to_string(path).ok()
        .and_then(|t| t.lines().find(|l| l.starts_with(key))
                       .and_then(|l| l.split_whitespace().nth(1)?.parse().ok()))
        .unwrap_or(0)
}

/// utime plus stime, in seconds. /proc/self/stat counts in 100 Hz ticks,
/// and the comm field can hold spaces, so the split starts after its ')'.
fn cpu() -> f64 {
    let Ok(t) = std::fs::read_to_string("/proc/self/stat") else { return 0.0 };
    let Some(i) = t.rfind(')') else { return 0.0 };
    let f: Vec<&str> = t[i + 1..].split_whitespace().collect();
    let get = |k: usize| f.get(k).and_then(|v| v.parse::<f64>().ok()).unwrap_or(0.0);
    (get(11) + get(12)) / 100.0
}

pub fn dump() {
    if !ON.load(Relaxed) { return; }
    let tot: f64 = NS.iter().map(|n| n.load(Relaxed) as f64).sum::<f64>() / 1e9;
    for i in 0..9 {
        let s = NS[i].load(Relaxed) as f64 / 1e9;
        if s > 0.0 {
            eprintln!("  phase {:<12} {:8.3} s  {:5.1}%   hwm {:5} MB",
                      NAME[i], s, 100.0 * s / if tot > 0.0 { tot } else { 1.0 },
                      HW[i].load(Relaxed) / 1024);
        }
    }
    let wall = T0.get().map_or(0.0, |t| t.elapsed().as_secs_f64());
    let (c, hwm) = (cpu(), field("/proc/self/status", "VmHWM:"));
    eprintln!("  cpu {c:.3} s over {wall:.3} s wall = {:.2} cores busy \
               ({} allowed), peak rss {:.2} GB",
              if wall > 0.0 { c / wall } else { 0.0 },
              THR.load(Relaxed), hwm as f64 / 1048576.0);
}
