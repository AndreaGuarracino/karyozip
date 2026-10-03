//! Reproducing KaryoScope's smoothing from an unsmoothed track, so the
//! smoothed variant never has to be stored. Every rule is c/karyozip.c's.

use crate::reader::Reader;
use std::io::Write;

pub struct Hier {
    node: Vec<Vec<u8>>,
    par: Vec<i32>,
    anc: Vec<u8>,
    lca: Vec<i16>,
    nlen1: Vec<usize>,   /* the --v1 name length, no compare a row */
    pub root: i32,
}

impl Hier {
    fn id(&mut self, nm: &[u8]) -> i32 {
        if let Some(i) = self.node.iter().position(|n| n == nm) { return i as i32; }
        if self.node.len() >= 4096 { crate::bail("too many hierarchy nodes"); }
        self.node.push(nm.to_vec());
        self.par.push(-1);
        (self.node.len() - 1) as i32
    }

    /// feature_set, child, parent; one feature set at a time.
    pub fn load(path: &str, fs: &[u8]) -> Hier {
        let data = crate::slurp(path, "cannot open the hierarchy file");
        let mut h = Hier { node: Vec::new(), par: Vec::new(), anc: Vec::new(),
                           lca: Vec::new(), nlen1: Vec::new(), root: -1 };
        for raw in data.split(|&c| c == b'\n') {
            let mut line = raw;
            if line.last() == Some(&b'\r') { line = &line[..line.len() - 1]; }
            if line.is_empty() { continue; }
            let cols: Vec<&[u8]> = line.split(|&c| c == b'\t').collect();
            if cols.len() < 3 { continue; }
            if cols.len() > 3 { crate::bail("hierarchy: a row has more than three columns"); }
            if cols[0] == b"feature_set" { continue; }
            if cols[0] != fs { continue; }
            let c = h.id(cols[1]);
            let p = h.id(cols[2]);
            if h.par[c as usize] >= 0 { crate::bail("hierarchy: a child is listed twice"); }
            h.par[c as usize] = p;
        }
        if h.node.is_empty() { crate::bail("feature set not found in the hierarchy"); }
        for a in 0..h.node.len() {
            if h.par[a] < 0 && h.node[a] != b"categorized" {
                crate::bail("hierarchy: a node without a parent that is not categorized");
            }
        }
        for a in 0..h.node.len() {
            let (mut c, mut n) = (a as i32, 0);
            while c >= 0 { n += 1; if n > h.node.len() { crate::bail("the hierarchy has a cycle"); }
                           c = h.par[c as usize]; }
        }
        h.root = h.id(b"categorized");
        h
    }

    /// Ancestor and lowest-common-ancestor tables, rebuilt after the track's
    /// own labels are added as nodes.
    pub fn build(&mut self) {
        let n = self.node.len();
        self.nlen1 = (0..n).map(|a| {
            let (l, nm) = (self.node[a].len(), &self.node[a]);
            if a as i32 != self.root && l > 9 && &nm[l - 9..] == b"_specific" { l - 9 } else { l }
        }).collect();
        self.anc = vec![0u8; n * n];
        self.lca = vec![-1i16; n * n];
        for a in 0..n {
            let mut c = a as i32;
            while c >= 0 { self.anc[a * n + c as usize] = 1; c = self.par[c as usize]; }
        }
        for a in 0..n {
            for b in 0..n {
                let mut got = -1i16;
                let mut c = a as i32;
                while c >= 0 {
                    if self.anc[b * n + c as usize] != 0 { got = c as i16; break; }
                    c = self.par[c as usize];
                }
                self.lca[a * n + b] = got;
            }
        }
    }

}

/// One sequence, in place. Promotes a run to the lowest common ancestor of
/// its neighbours when the noise between them is below it.
fn smooth_one(h: &Hier, st: &[u64], en: &[u64], feat: &mut [i32], maxgap: u64) {
    let n = feat.len();
    if n == 0 { return; }
    let nn = h.node.len();
    // every index below is proven once here: st, en, was are n long, every
    // feature id is a node, dis and par are nn long, anc and lca nn*nn;
    // the loops then read unchecked, as the C's ISANC macro does
    assert!(st.len() == n && en.len() == n && h.par.len() == nn
            && h.anc.len() == nn * nn && h.lca.len() == nn * nn
            && feat.iter().all(|&f| (f as usize) < nn)
            && h.par.iter().all(|&p| p < nn as i32));
    let (anc, par, lca) = (&h.anc[..], &h.par[..], &h.lca[..]);
    macro_rules! at { ($v:expr, $i:expr) => { unsafe { *$v.get_unchecked($i) } } }
    let isanc = |node: i32, pa: i32| at!(anc, node as usize * nn + pa as usize) != 0;
    let mut was = vec![0u8; n];
    let mut dis = vec![0u32; nn];
    let mut rel: Vec<usize> = Vec::with_capacity(n);
    let mut mark = 0u32;
    loop {
        let mut changes = false;
        was.iter_mut().for_each(|x| *x = 0);
        let mut i = 0usize;
        while i < n {
            let ws = i;
            rel.clear();
            rel.push(i);
            let mut firstun: i64 = -1;
            mark += 1;
            let (mut lrf, mut lri, mut lre) = (at!(feat, i), i, at!(en, i));
            let mut j = i;
            while j + 1 < n {
                let nf = at!(feat, j + 1);
                if at!(st, j + 1) - lre > maxgap { break; }
                if isanc(lrf, nf) {
                    if at!(dis, nf as usize) == mark { break; }
                    rel.push(j + 1);
                    lri = j + 1; lrf = nf; lre = at!(en, j + 1); j += 1;
                } else if !isanc(nf, lrf) {
                    let mut c = nf;
                    while c >= 0 { unsafe { *dis.get_unchecked_mut(c as usize) = mark; } c = at!(par, c as usize); }
                    if firstun < 0 && at!(was, j + 1) == 0 { firstun = (j + 1) as i64; }
                    j += 1;
                } else { break; }
            }
            let mut we = lri;
            let peak = at!(feat, we);
            let mut k = we;
            while k + 1 < n {
                let nf = at!(feat, k + 1);
                if at!(st, k + 1) - lre > maxgap { break; }
                if isanc(nf, lrf) {
                    rel.push(k + 1);
                    lri = k + 1; lrf = nf; lre = at!(en, k + 1); k += 1;
                } else if !isanc(lrf, nf) {
                    if isanc(nf, peak) { break; }
                    if firstun < 0 && at!(was, k + 1) == 0 { firstun = (k + 1) as i64; }
                    k += 1;
                } else { break; }
            }
            we = lri;
            let mut nrel = rel.len();
            if nrel >= 2 { nrel -= 1; }   // archive quirk, kept
            for &t in &rel[..nrel] { was[t] = 1; }
            if we > ws {
                let (a, b) = (at!(feat, ws) as usize, at!(feat, we) as usize);
                let l = at!(lca, a * nn + b) as i32;
                if l >= 0 {
                    for x in feat.iter_mut().take(we).skip(ws + 1) {
                        if *x != l && isanc(l, *x) { *x = l; changes = true; }
                    }
                }
            }
            if firstun >= 0 {
                i = if firstun as usize > we && we > ws { we } else { firstun as usize };
            } else if we > i { i = we; }
            else { while i < n && was[i] != 0 { i += 1; } }
        }
        if !changes { break; }
    }
}

/// Merges equal neighbours and prints, dropping the v1 `_specific` suffix
/// when asked.
#[allow(clippy::too_many_arguments)]
fn emit(h: &Hier, out: &mut Vec<u8>, seq: &[u8], st: &[u64], en: &[u64],
        feat: &[i32], novel: &[u8], lo: u64, hi: u64, v1: bool) {
    let n = feat.len();
    if n == 0 { return; }
    let (mut cs, mut ce, mut cf, mut cn) = (st[0], en[0], feat[0], novel[0]);
    for k in 1..=n {
        let same = k < n && feat[k] == cf && st[k] == ce
                   && (cf != h.root || novel[k] == cn);
        if same { ce = en[k]; continue; }
        if (hi == 0 && lo == 0) || (ce > lo && (hi == 0 || cs < hi)) {
            let nm: &[u8] = if cf == h.root && cn != 0 { b"novel" } else { &h.node[cf as usize] };
            let ml = if cf == h.root && cn != 0 { 5 }
                     else if v1 { h.nlen1[cf as usize] } else { nm.len() };
            let mut b = [0u8; 44];      // \tSTART\tEND\t in one piece
            b[0] = b'\t';
            let mut i = 1 + crate::kst::putu64_at(&mut b[1..], cs);
            b[i] = b'\t'; i += 1;
            i += crate::kst::putu64_at(&mut b[i..], ce);
            b[i] = b'\t'; i += 1;
            // one reserve for the whole row, then one pointer through it
            out.reserve(seq.len() + i + ml + 1);
            unsafe {
                let n0 = out.len();
                let p = out.as_mut_ptr().add(n0);
                std::ptr::copy_nonoverlapping(seq.as_ptr(), p, seq.len());
                std::ptr::copy_nonoverlapping(b.as_ptr(), p.add(seq.len()), i);
                let q = p.add(seq.len() + i);
                std::ptr::copy_nonoverlapping(nm.as_ptr(), q, ml);
                *q.add(ml) = b'\n';
                out.set_len(n0 + seq.len() + i + ml + 1);
            }
        }
        if k == n { break; }
        cs = st[k]; ce = en[k]; cf = feat[k]; cn = novel[k];
    }
}

#[allow(clippy::too_many_arguments)]
pub fn run(path: &str, hier: &str, fs: &[u8], seq: Option<&[u8]>,
           lo: u64, hi: u64, v1: bool, maxgap: u64,
           ro: &crate::reader::ReadOpts, feats: Option<&crate::db::Feats>,
           out: &mut (impl Write + Send)) {
    let mut h = Hier::load(hier, fs);
    let mut r = Reader::open_with(path, ro);
    if let Some(f) = feats { f.check_hash(&r); }
    let setk = feats.map(|f| f.sets.iter().position(|s| s == fs)
        .unwrap_or_else(|| crate::bail("smooth: the set is not in features.tsv")));
    // The track's own labels join the hierarchy, so the tables are rebuilt.
    // Under --aln "=" and under --db a non-featureID are resolved per row.
    let (mut labmap, mut novelid) = (Vec::new(), Vec::new());
    for i in 0..r.nlabels() {
        let mut l = r.label_at(i).to_vec();
        if r.has_aln() && l == b"=" { labmap.push(i32::MIN); novelid.push(0); continue; }
        if let Some(f) = feats {
            if l.is_empty() || !l.iter().all(u8::is_ascii_digit) { labmap.push(i32::MIN + 1); novelid.push(0); continue; }
            let id = f.featid(&l);
            l = f.lab(id, setk.unwrap()).to_vec();
        }
        novelid.push(u8::from(l == b"novel"));
        labmap.push(h.id(&l));
    }
    h.build();   // after the labels: they add the unknown nodes

    let (first, last) = match seq {
        Some(s) => { let i = r.find(s); if i < 0 { crate::bail("sequence not found"); }
                     (i as u64, i as u64 + 1) }
        None => (0, r.nreads),
    };
    // One sequence: smoothed whole, printed where it overlaps the range.
    // The phases are the C's serial ones; its threaded path records none.
    let one = |r: &mut Reader, ri: u64, sr: &mut crate::reader::SeqRows,
               buf: &mut Vec<u8>, prof: bool| {
        let s = crate::prof::span(if prof { 6 } else { 10 });
        r.seq_rows(ri, &labmap, &novelid, h.root, sr);
        drop(s);
        if hi > 0 && sr.en.last().is_none_or(|&e| lo >= e) { return; }
        let s = crate::prof::span(if prof { 4 } else { 10 });
        smooth_one(&h, &sr.st, &sr.en, &mut sr.feat, maxgap);
        drop(s);
        let name = r.name_at(ri).to_vec();
        let _s = crate::prof::span(if prof { 5 } else { 10 });
        emit(&h, buf, &name, &sr.st, &sr.en, &sr.feat, &sr.nov, lo, hi, v1);
    };
    let wr = |out: &mut dyn Write, b: &[u8]| out.write_all(b).unwrap_or_else(|_| crate::bail("write failed"));
    let (nthr, chunk) = (rayon::current_num_threads() as u64, 4096u64);
    if nthr > 1 && last - first >= 64 {
        // As the C's smooth_mt: chunks of sequences on worker threads, each
        // with its own reader (sharing a one-block track's block) and one buffer it reuses, taken in order and
        // written in turn, so nothing is held but a buffer per thread. A
        // window of eight chunks per thread held 2 GB against the C's 0.45 GB.
        use std::sync::atomic::{AtomicU64, Ordering};
        let blk = r.share_block();   // one block: decoded once, for every thread, before the names
        let names = r.names_shared();
        let nchunk = (last - first).div_ceil(chunk);
        let next = AtomicU64::new(0);
        let turn = std::sync::Mutex::new((0u64, out));
        let bell = std::sync::Condvar::new();
        let one = &one;
        std::thread::scope(|sc| {
            for _ in 0..nthr {
                let (names, blk, next, turn, bell) = (names.clone(), blk.clone(), &next, &turn, &bell);
                sc.spawn(move || {
                    let mut r = Reader::open_with(path, ro);
                    r.adopt_names(names);
                    if let Some(b) = blk { r.adopt_block(b); }
                    let mut buf = Vec::with_capacity(1 << 20);
                    let mut sr = crate::reader::SeqRows::default();
                    loop {
                        let c = next.fetch_add(1, Ordering::Relaxed);
                        if c >= nchunk { break; }
                        let (a, b) = (first + c * chunk, (first + (c + 1) * chunk).min(last));
                        buf.clear();
                        for ri in a..b { one(&mut r, ri, &mut sr, &mut buf, false); }
                        let mut g = turn.lock().unwrap_or_else(|e| e.into_inner());
                        while g.0 != c { g = bell.wait(g).unwrap_or_else(|e| e.into_inner()); }
                        wr(&mut *g.1, &buf);
                        g.0 += 1;
                        drop(g);
                        bell.notify_all();
                    }
                });
            }
        });
        return;
    }
    let mut buf = Vec::with_capacity(1 << 20);
    let mut sr = crate::reader::SeqRows::default();
    for ri in first..last {
        one(&mut r, ri, &mut sr, &mut buf, true);
        if buf.len() > (1 << 20) - 4096 { wr(out, &buf); buf.clear(); }
    }
    wr(out, &buf);
}
