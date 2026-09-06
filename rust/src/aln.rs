//! --aln: labels and run ends the reference already implies are not stored.
//! The alignment is read for one thing only, the collinear M/=/X segments.

use crate::kst::{self, unzz, vget, FNV0};
use std::collections::HashMap;

/// One collinear block: where it starts on the query, which reference
/// contig, how long, where it starts there, and the strand.
#[derive(Clone, Copy)]
pub struct Seg { pub q: u64, pub t: u32, pub n: u32, pub tp: u64, pub rev: bool }

/// A reference contig: run starts, its end, and a label id per run.
#[derive(Clone)]
pub struct Rcon { pub name: Vec<u8>, pub st: Vec<u64>, pub end: u64, pub lab: Vec<u32> }

#[derive(Default, Clone)]
pub struct Aln {
    pub rc: Vec<Rcon>,
    rcix: HashMap<Vec<u8>, u32>,
    pub qs: Vec<(Vec<u8>, Vec<Seg>)>,
    qsix: HashMap<Vec<u8>, usize>,
    pub labs: Vec<Vec<u8>>,
    labix: HashMap<Vec<u8>, u32>,
    /// --aln-names: ordinal N is line N+1, for tracks named by ordinal.
    pub ord: Vec<Vec<u8>>,
    pub hash: u64,
    /// the last segment answered, so a walk in order does not re-search
    cur: std::sync::Arc<std::sync::atomic::AtomicU64>,
}

impl Aln {
    fn lab(&mut self, s: &[u8]) -> u32 {
        if let Some(&i) = self.labix.get(s) { return i; }
        let i = self.labs.len() as u32;
        self.labs.push(s.to_vec());
        self.labix.insert(s.to_vec(), i);
        i
    }

    /// The reference track: one contig's rows together, sorted, gapless.
    pub fn load_ref(&mut self, path: &str) {
        let mut lines = crate::Lines::open(path, "cannot open the reference track");
        let mut rh = FNV0;
        let (mut badln, mut cur) = (0u64, usize::MAX);
        while let Some(raw) = lines.next_line() {
            let mut line = raw;
            if line.last() == Some(&b'\r') { line = &line[..line.len() - 1]; }
            let c: Vec<&[u8]> = line.splitn(4, |&x| x == b'\t').collect();
            if c.len() < 4 { if !line.is_empty() { badln += 1; } continue; }
            rh = kst::fnv(rh, line);
            rh = kst::fnv(rh, b"\n");
            if c[0].is_empty() { crate::bail("reference track: empty contig name"); }
            if cur == usize::MAX || self.rc[cur].name != c[0] {
                if self.rcix.contains_key(c[0]) {
                    crate::bail("reference track: a contig's rows are not together");
                }
                self.rcix.insert(c[0].to_vec(), self.rc.len() as u32);
                self.rc.push(Rcon { name: c[0].to_vec(), st: Vec::new(), end: 0, lab: Vec::new() });
                cur = self.rc.len() - 1;
            }
            let st = crate::bedcoord(c[1], "reference start is not a number");
            let en = crate::bedcoord(c[2], "reference end is not a number");
            let want = if self.rc[cur].st.is_empty() { 0 } else { self.rc[cur].end };
            if st != want { crate::bail("reference track: rows are not sorted and gapless from 0"); }
            if en <= st { crate::bail("reference track: a run of no length"); }
            let l = self.lab(c[3]);
            let r = &mut self.rc[cur];
            r.st.push(st);
            r.lab.push(l);
            r.end = en;
        }
        self.hash ^= rh;
        if badln > 0 { eprintln!("reference: {badln} lines without four columns skipped"); }
        let tot: usize = self.rc.iter().map(|x| x.st.len()).sum();
        eprintln!("reference: {} contigs, {tot} runs", self.rc.len());
    }

    /// The PAF, read for its cg:Z: CIGARs. Everything else is dead weight.
    pub fn load_paf(&mut self, path: &str) {
        let mut lines = crate::Lines::open(path, "cannot open the PAF");
        let (mut nblk, mut badpaf, mut nocg, mut noref) = (0usize, 0u64, 0u64, 0u64);
        while let Some(raw) = lines.next_line() {
            let mut line = raw;
            if line.last() == Some(&b'\r') { line = &line[..line.len() - 1]; }
            if line.is_empty() { continue; }
            let fl: Vec<&[u8]> = line.split(|&c| c == b'\t').collect();
            if fl.len() < 13 { badpaf += 1; continue; }
            let cg = match fl[12..].iter().find(|t| t.starts_with(b"cg:Z:")) {
                Some(t) => &t[5..],
                None => { nocg += 1; continue; }
            };
            if fl[4].len() != 1 || (fl[4][0] != b'+' && fl[4][0] != b'-') {
                crate::bail("PAF strand is not + or -");
            }
            let rev = fl[4][0] == b'-';
            let qs = crate::bedcoord(fl[2], "PAF query start is not a number");
            let qe = crate::bedcoord(fl[3], "PAF query end is not a number");
            let mut tp = crate::bedcoord(fl[7], "PAF target start is not a number");
            if qe < qs { crate::bail("PAF query end is before its start"); }
            let mut qp = if rev { qe } else { qs };
            let tix = match self.rcix.get(fl[5]) { Some(&i) => i, None => { noref += 1; continue; } };
            let mut segs: Vec<Seg> = Vec::new();
            let (mut qcons, mut tcons) = (0u64, 0u64);
            let mut i = 0usize;
            while i < cg.len() {
                if !cg[i].is_ascii_digit() { crate::bail("PAF CIGAR: an operation has no length"); }
                let mut v: u64 = 0;
                while i < cg.len() && cg[i].is_ascii_digit() {
                    if v > (u64::MAX - 9) / 10 { crate::bail("PAF CIGAR: a length overflows"); }
                    v = v * 10 + (cg[i] - b'0') as u64;
                    i += 1;
                }
                if v == 0 { crate::bail("PAF CIGAR: a zero-length operation"); }
                if v > 0xffff_ffff { crate::bail("PAF CIGAR: an operation longer than 4 GB"); }
                if i >= cg.len() { crate::bail("PAF CIGAR: an operation has no length"); }
                let op = cg[i];
                i += 1;
                match op {
                    b'M' | b'=' | b'X' => {
                        if rev && qp < v { crate::bail("PAF CIGAR runs before the query start"); }
                        segs.push(Seg { q: if rev { qp - v } else { qp }, t: tix,
                                        n: v as u32, tp, rev });
                        qp = if rev { qp - v } else { qp + v };
                        tp += v; qcons += v; tcons += v; nblk += 1;
                    }
                    b'I' | b'S' => {
                        if rev && qp < v { crate::bail("PAF CIGAR runs before the query start"); }
                        qp = if rev { qp - v } else { qp + v };
                        qcons += v;
                    }
                    b'D' | b'N' => { tp += v; tcons += v; }
                    b'H' | b'P' => {}
                    _ => crate::bail("PAF CIGAR: unknown operation"),
                }
            }
            if qcons != qe - qs {
                crate::bail("PAF CIGAR consumes a different query length than the record says");
            }
            let te = crate::bedcoord(fl[8], "PAF target end is not a number");
            let ts = tp - tcons;
            if te < ts || tcons != te - ts {
                crate::bail("PAF CIGAR consumes a different target length than the record says");
            }
            if segs.is_empty() { continue; }
            match self.qsix.get(fl[0]) {
                Some(&k) => self.qs[k].1.extend_from_slice(&segs),
                None => {
                    self.qsix.insert(fl[0].to_vec(), self.qs.len());
                    self.qs.push((fl[0].to_vec(), segs));
                }
            }
        }
        for (_, s) in self.qs.iter_mut() {
            s.sort_by_key(|g| (g.q, g.t, g.tp, g.rev, g.n));
        }
        // The binding hash: this PAF and this reference, or nothing.
        self.seghash();
        if badpaf > 0 || nocg > 0 || noref > 0 {
            eprintln!("alignment: skipped {badpaf} malformed PAF lines, {nocg} without a cg:Z: tag, {noref} whose target is not in the reference");
        }
        if nblk == 0 {
            if noref > 0 { crate::bail("alignment: no PAF record maps to a contig of the reference track"); }
            crate::bail("alignment: no usable PAF record, so --aln would silently do nothing. Check the PAF carries cg:Z: CIGAR tags");
        }
        eprintln!("alignment: {} query contigs, {nblk} collinear blocks", self.qs.len());
    }

    /// The alignment's own hash, which binds a track to its PAF.
    fn seghash(&mut self) {
        let mut sum = 0u64;
        for (nm, segs) in &self.qs {
            let mut h = kst::fnv(FNV0, nm);
            h = kst::fnv(h, &[0]);
            for g in segs {
                h = kst::fnv(h, &g.q.to_le_bytes());
                h = kst::fnv(h, &g.t.to_le_bytes());
                h = kst::fnv(h, &g.n.to_le_bytes());
                h = kst::fnv(h, &g.tp.to_le_bytes());
                h = kst::fnv(h, &(g.rev as i32).to_le_bytes());
            }
            if !segs.is_empty() { sum = sum.wrapping_add(h); }
        }
        self.hash ^= sum;
    }

    /// A packed alignment if the file is one, else a PAF, as the C does.
    pub fn load_aln(&mut self, path: &str) {
        if !self.load_packed(path) { self.load_paf(path); }
    }

    /// aln-pack's own format: five coded streams behind a magic and the
    /// reference's contig-order hash. Returns false when it is not one.
    fn load_packed(&mut self, path: &str) -> bool {
        let Ok(d) = std::fs::read(path) else { return false };
        if d.len() < 6 || &d[..6] != b"KSTAL1" { return false; }
        if d.len() < 14 { crate::bail("truncated packed alignment"); }
        let rh = u64::from_le_bytes(d[6..14].try_into().unwrap());
        if rh != refhash(self) {
            crate::bail("this packed alignment was made against another reference, or another contig order");
        }
        let body = &d[14..];
        let mut i = 0usize;
        let mut sec: Vec<Vec<u8>> = Vec::with_capacity(5);
        for _ in 0..5 {
            let (cd, j) = vget(body, i);
            let (rl, j) = vget(body, j);
            let (cl, j) = vget(body, j);
            let cl = cl as usize;
            if cl > body.len() - j { crate::bail("corrupt packed alignment: a stream overruns the file"); }
            let src = &body[j..j + cl];
            i = j + cl;
            sec.push(match cd {
                0 => Vec::new(),
                2 => crate::reader::tok_d(src),
                1 => crate::ffi::zstd_decompress(src, rl as usize)
                        .unwrap_or_else(|| crate::bail("corrupt packed alignment: a stream does not decompress")),
                _ => crate::ffi::paydecomp(src, (cd - 3) as i32, rl as usize)
                        .unwrap_or_else(|| crate::bail("corrupt packed alignment: a stream does not decompress")),
            });
        }
        let (nms, cnt, pos, len, tgt) = (&sec[0], &sec[1], &sec[2], &sec[3], &sec[4]);
        let (mut ci, mut pi, mut li, mut si) = (0usize, 0usize, 0usize, 0usize);
        let nrc = self.rc.len() as u64;
        let mut nb = 0usize;
        let mut at = 0usize;
        while at < nms.len() {
            let Some(k) = memchr::memchr(b'\n', &nms[at..]) else {
                crate::bail("corrupt packed alignment: a name has no newline")
            };
            let name = nms[at..at + k].to_vec();
            at += k + 1;
            let (cn, j) = vget(cnt, ci); ci = j;
            let qi = *self.qsix.entry(name.clone()).or_insert_with(|| {
                self.qs.push((name, Vec::new()));
                self.qs.len() - 1
            });
            let segs = &mut self.qs[qi].1;
            segs.clear();
            segs.reserve(cn as usize);
            let (mut pq, mut pt, mut tx) = (0u64, 0u64, 0u32);
            let (mut rev, mut fresh, mut left) = (false, true, 0u64);
            for _ in 0..cn {
                if left == 0 {
                    let (xr, j) = vget(tgt, si);
                    rev = xr & 1 != 0;
                    if xr >> 1 >= nrc { crate::bail("corrupt packed alignment: a target is not in the reference"); }
                    tx = (xr >> 1) as u32;
                    let (l, j2) = vget(tgt, j); si = j2;
                    if l == 0 { crate::bail("corrupt packed alignment: an empty target run"); }
                    left = l;
                    fresh = true;
                }
                let (q, tp);
                if fresh {
                    let (a, j) = vget(pos, pi);
                    let (b, j2) = vget(pos, j); pi = j2;
                    q = a; tp = b; fresh = false;
                } else {
                    let (a, j) = vget(pos, pi);
                    let (b, j2) = vget(pos, j); pi = j2;
                    q = pq.wrapping_add(unzz(a));
                    tp = pt.wrapping_add(unzz(b));
                }
                let (n, j) = vget(len, li); li = j;
                segs.push(Seg { q, t: tx, n: n as u32, tp, rev });
                pq = q + n; pt = tp + n; left -= 1;
            }
            nb += segs.len();
        }
        eprintln!("packed alignment: {} query contigs, {nb} collinear blocks", self.qs.len());
        self.seghash();
        true
    }

    /// A bare ordinal resolves through the map when one was given.
    pub fn find(&self, nm: &[u8]) -> Option<usize> {
        if !self.ord.is_empty() && !nm.is_empty() && nm.iter().all(|c| c.is_ascii_digit()) {
            let i: usize = std::str::from_utf8(nm).unwrap().parse().ok()?;
            if i >= self.ord.len() { return None; }
            return self.qsix.get(&self.ord[i]).copied();
        }
        self.qsix.get(nm).copied()
    }

    /// One name a line, in the order the annotator produced them.
    pub fn load_ordinals(&mut self, path: &str) {
        let data = crate::slurp(path, "cannot open the file given to --names");
        let mut h = FNV0;
        for line in data.split(|&c| c == b'\n') {
            let mut l = line;
            if l.last() == Some(&b'\r') { l = &l[..l.len() - 1]; }
            if self.ord.is_empty() && (l.first() == Some(&b'>') || l.first() == Some(&b'@')) {
                crate::bail("--names on encode takes a plain list, one name a line");
            }
            if l.is_empty() && self.ord.len() as u64 > 0 { continue; }
            self.ord.push(l.to_vec());
            h = kst::fnv(h, l);
            h = kst::fnv(h, &[0]);
        }
        self.hash ^= h;
        eprintln!("ordinal map: {} names", self.ord.len());
    }

    /// The last segment whose query start is at or before `pos`.
    fn segat(&self, qi: usize, pos: u64) -> Option<&Seg> {
        let s = &self.qs[qi].1;
        if s.is_empty() { return None; }
        let packed = self.cur.load(std::sync::atomic::Ordering::Relaxed);
        let (cq, ci) = ((packed >> 32) as usize, (packed & 0xffff_ffff) as usize);
        let mut lo = usize::MAX;
        if cq == qi {
            let mut l = ci;
            let mut k = 0;
            while l < s.len() && s[l].q <= pos && k < 8 { l += 1; k += 1; }
            while l > 0 && s[l - 1].q > pos && k < 16 { l -= 1; k += 1; }
            if k < 16 && !(l < s.len() && s[l].q <= pos) { lo = l; }
        }
        if lo == usize::MAX {
            lo = s.partition_point(|g| g.q <= pos);
        }
        self.cur.store(((qi as u64) << 32) | lo as u64, std::sync::atomic::Ordering::Relaxed);
        if lo == 0 { return None; }
        let g = &s[lo - 1];
        if pos < g.q + g.n as u64 { Some(g) } else { None }
    }

    fn runat(rc: &Rcon, tpos: u64) -> usize {
        if tpos >= rc.end { return 0; }
        rc.st.partition_point(|&x| x <= tpos)
    }

    pub fn label(&self, qi: usize, pos: u64) -> u32 {
        let Some(g) = self.segat(qi, pos) else { return u32::MAX };
        let tpos = if g.rev { g.tp + (g.q + g.n as u64 - 1 - pos) } else { g.tp + (pos - g.q) };
        let rc = &self.rc[g.t as usize];
        let a = Self::runat(rc, tpos);
        if a > 0 { rc.lab[a - 1] } else { u32::MAX }
    }

    /// The query position where the reference run under `pos` ends.
    pub fn end(&self, qi: usize, pos: u64) -> (u64, bool) {
        let Some(g) = self.segat(qi, pos) else { return (0, false) };
        let rev = g.rev;
        let tpos = if rev { g.tp + (g.q + g.n as u64 - 1 - pos) } else { g.tp + (pos - g.q) };
        let rc = &self.rc[g.t as usize];
        let a = Self::runat(rc, tpos);
        if a == 0 { return (0, rev); }
        let tend = if a < rc.st.len() { rc.st[a] } else { rc.end };
        if rev {
            let ts = rc.st[a - 1];
            if ts < g.tp { return (0, rev); }
            return (g.q + g.n as u64 - (ts - g.tp), rev);
        }
        if tend < g.tp || tend - g.tp > g.n as u64 { return (0, rev); }
        (g.q + (tend - g.tp), rev)
    }

    /// Three probes inside the run decide its label, as the C does.
    pub fn row(&self, qi: usize, st: u64, len: u64) -> u32 {
        if len < 2 { return self.label(qi, st); }
        let a = self.label(qi, st + len / 4);
        let b = self.label(qi, st + len / 2);
        if a == b { return a; }
        if a == self.label(qi, st + 3 * len / 4) { a } else { b }
    }
}

/* --------------------------------------------------------- aln-pack */

/// The packed alignment: five streams, each given its own codec search.
/// A PAF is far larger than what --aln actually reads.
pub fn pack(a: &Aln, out: &str, refhash: u64) {
    use crate::ffi;
    let (mut names, mut cnt, mut pos, mut len, mut seg_) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let mut nb = 0usize;
    for (nm, segs) in &a.qs {
        if segs.is_empty() { continue; }
        names.extend_from_slice(nm);
        names.push(b'\n');
        kst::put(&mut cnt, segs.len() as u64);
        let (mut pq, mut pt, mut px) = (0u64, 0u64, 0u32);
        let mut prev_rev: i32 = -1;
        let mut run = 0u64;
        for (j, g) in segs.iter().enumerate() {
            let newrun = g.t != px || (g.rev as i32) != prev_rev;
            if j == 0 || newrun {
                kst::put(&mut pos, g.q);
                kst::put(&mut pos, g.tp);
            } else {
                kst::put(&mut pos, kst::zz(g.q as i64 - pq as i64));
                kst::put(&mut pos, kst::zz(g.tp as i64 - pt as i64));
            }
            if newrun {
                if j > 0 { kst::put(&mut seg_, run); }
                kst::put(&mut seg_, g.t as u64 * 2 + g.rev as u64);
                px = g.t; prev_rev = g.rev as i32; run = 0;
            }
            run += 1;
            kst::put(&mut len, g.n as u64);
            pq = g.q + g.n as u64;
            pt = g.tp + g.n as u64;
        }
        kst::put(&mut seg_, run);
        nb += segs.len();
    }
    const SNAME: [&str; 5] = ["names", "counts", "positions", "lengths", "targets"];
    let st = [names, cnt, pos, len, seg_];
    let mut all = Vec::new();
    for (k, raw) in st.iter().enumerate() {
        let (mut best, mut bc): (Vec<u8>, i64) = (Vec::new(), -1);
        if !raw.is_empty() {
            if k == 0 {
                if let Some(t) = crate::writer::tok_best(raw, false, None) { best = t; bc = 0x7fff_0000; }
            }
            if bc < 0 {
                best = ffi::zstd(raw, 19);
                bc = 0;
                for c in crate::writer::CAND_FAST {
                    if let Some(t) = ffi::paycomp(raw, c) {
                        if t.len() < best.len() { best = t; bc = c as i64; }
                    }
                }
            }
        }
        eprintln!("  aln {:<9} raw {:>10}  coded {:>10}  codec 0x{:x}",
                  SNAME[k], raw.len(), best.len(), if bc < 0 { 0 } else { bc });
        kst::put(&mut all, if bc < 0 { 0 } else if bc == 0x7fff_0000 { 2 }
                           else if bc == 0 { 1 } else { bc as u64 + 3 });
        kst::put(&mut all, raw.len() as u64);
        kst::put(&mut all, best.len() as u64);
        all.extend_from_slice(&best);
    }
    let mut f = Vec::with_capacity(all.len() + 16);
    f.extend_from_slice(b"KSTAL1");
    f.extend_from_slice(&refhash.to_le_bytes());
    f.extend_from_slice(&all);
    std::fs::write(out, &f).unwrap_or_else(|_| crate::bail("cannot write the packed alignment"));
    let nq = a.qs.iter().filter(|(_, s)| !s.is_empty()).count();
    eprintln!("packed alignment: {nq} queries, {nb} blocks, {} B ({:.2} B/query)",
              f.len(), if nq > 0 { f.len() as f64 / nq as f64 } else { 0.0 });
}

/// The reference contig names, in order: what binds a packed alignment.
pub fn refhash(a: &Aln) -> u64 {
    let mut h = FNV0;
    for c in &a.rc { h = kst::fnv(h, &c.name); h = kst::fnv(h, b"\n"); }
    h
}
