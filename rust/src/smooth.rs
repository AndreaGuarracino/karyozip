//! Reproducing KaryoScope's smoothing from an unsmoothed track, so the
//! smoothed variant never has to be stored. Every rule is c/karyozip.c's.

use crate::reader::Reader;
use std::io::Write;

pub struct Hier {
    node: Vec<Vec<u8>>,
    par: Vec<i32>,
    anc: Vec<u8>,
    lca: Vec<i16>,
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
                           lca: Vec::new(), root: -1 };
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

    fn is_anc(&self, node: i32, pa: i32) -> bool {
        self.anc[node as usize * self.node.len() + pa as usize] != 0
    }
}

/// One sequence, in place. Promotes a run to the lowest common ancestor of
/// its neighbours when the noise between them is below it.
fn smooth_one(h: &Hier, st: &[u64], en: &[u64], feat: &mut [i32], maxgap: u64) {
    let n = feat.len();
    if n == 0 { return; }
    let nn = h.node.len();
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
            let (mut lrf, mut lri, mut lre) = (feat[i], i, en[i]);
            let mut j = i;
            while j + 1 < n {
                let nf = feat[j + 1];
                if st[j + 1] - lre > maxgap { break; }
                if h.is_anc(lrf, nf) {
                    if dis[nf as usize] == mark { break; }
                    rel.push(j + 1);
                    lri = j + 1; lrf = nf; lre = en[j + 1]; j += 1;
                } else if !h.is_anc(nf, lrf) {
                    let mut c = nf;
                    while c >= 0 { dis[c as usize] = mark; c = h.par[c as usize]; }
                    if firstun < 0 && was[j + 1] == 0 { firstun = (j + 1) as i64; }
                    j += 1;
                } else { break; }
            }
            let mut we = lri;
            let peak = feat[we];
            let mut k = we;
            while k + 1 < n {
                let nf = feat[k + 1];
                if st[k + 1] - lre > maxgap { break; }
                if h.is_anc(nf, lrf) {
                    rel.push(k + 1);
                    lri = k + 1; lrf = nf; lre = en[k + 1]; k += 1;
                } else if !h.is_anc(lrf, nf) {
                    if h.is_anc(nf, peak) { break; }
                    if firstun < 0 && was[k + 1] == 0 { firstun = (k + 1) as i64; }
                    k += 1;
                } else { break; }
            }
            we = lri;
            let mut nrel = rel.len();
            if nrel >= 2 { nrel -= 1; }   // archive quirk, kept
            for &t in &rel[..nrel] { was[t] = 1; }
            if we > ws {
                let l = h.lca[feat[ws] as usize * nn + feat[we] as usize] as i32;
                if l >= 0 {
                    for x in feat.iter_mut().take(we).skip(ws + 1) {
                        if *x != l && h.is_anc(l, *x) { *x = l; changes = true; }
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
            let mut ml = nm.len();
            if v1 && cf != h.root && ml > 9 && &nm[ml - 9..] == b"_specific" { ml -= 9; }
            out.extend_from_slice(seq);
            out.push(b'\t');
            crate::kst::putu64(out, cs);
            out.push(b'\t');
            crate::kst::putu64(out, ce);
            out.push(b'\t');
            out.extend_from_slice(&nm[..ml]);
            out.push(b'\n');
        }
        if k == n { break; }
        cs = st[k]; ce = en[k]; cf = feat[k]; cn = novel[k];
    }
}

#[allow(clippy::too_many_arguments)]
pub fn run(path: &str, hier: &str, fs: &[u8], seq: Option<&[u8]>,
           lo: u64, hi: u64, v1: bool, maxgap: u64,
           ro: &crate::reader::ReadOpts, feats: Option<&crate::db::Feats>, out: &mut impl Write) {
    let mut h = Hier::load(hier, fs);
    h.build();
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
            l = if id == 0 { b"novel".to_vec() } else { f.id2lab[&id][setk.unwrap()].clone() };
        }
        novelid.push(u8::from(l == b"novel"));
        labmap.push(h.id(&l));
    }
    h.build();

    let (first, last) = match seq {
        Some(s) => { let i = r.find_seq(s); if i < 0 { crate::bail("sequence not found"); }
                     (i as u64, i as u64 + 1) }
        None => (0, r.nreads),
    };
    let mut buf = Vec::with_capacity(1 << 20);
    for ri in first..last {
        let _s = crate::prof::span(6);
        let (st, en, mut feat, novel) = r.seq_rows(ri, &labmap, &novelid, h.root);
        drop(_s);
        // The whole sequence is smoothed; a range only cuts the output.
        if hi > 0 && en.last().is_none_or(|&e| lo >= e) { continue; }
        let _s = crate::prof::span(4);
        smooth_one(&h, &st, &en, &mut feat, maxgap);
        drop(_s);
        let name = r.name_at(ri).to_vec();
        let _s = crate::prof::span(5);
        emit(&h, &mut buf, &name, &st, &en, &feat, &novel, lo, hi, v1);
        drop(_s);
        if buf.len() > (1 << 20) - 4096 {
            out.write_all(&buf).unwrap_or_else(|_| crate::bail("write failed"));
            buf.clear();
        }
    }
    out.write_all(&buf).unwrap_or_else(|_| crate::bail("write failed"));
}
