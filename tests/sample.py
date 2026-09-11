#!/usr/bin/env python3
# The wide bench's inputs. Prints counts only, never a row or a name.
#   sample.py pick INVENTORY N STAGED > wide.tsv
#     up to N complete track sets a class, version and backend from the
#     inventory, one a sample, fixed seed; then the staged cells, which alone
#     carry alignments. A line: cell, source directory, alignment or -.
#   sample.py cut SRC OUT PAF ROWS BGZIP
#     every track of SRC to at most about ROWS rows, the same sequences in
#     every track: reads whole, picked by a hash of the name; an assembly, whose
#     few sequences are each bigger than the cut, the same leading window of
#     every sequence, ending where the largest track reaches its share of rows,
#     ends clipped and a file named window left beside the tracks. Ordinal cells renumbered from 0 with their
#     names.txt rewritten; the alignment keeps the picked reads. bgzip level 6,
#     the level of the originals, so the ratio divides by the same kind of file.
import sys, os, zlib, random, subprocess, collections, multiprocessing
STD = {'acrocentric', 'chromosome', 'gene', 'region', 'repeat', 'subtelomeric', 'haplotype'}

def parse(f):                  # db, set, kind of a plain track, or None for a variant
    t = f[:-7]; t = t[:-9] if t.endswith('.features') else t
    db, st, kind = ([''] * 3 + t.split('.'))[-3:]
    ok = db.split('_human')[0] in ('KS', 'HKS') and st in STD and kind in ('presmoothed', 'smoothed')
    return (db, st, kind) if ok and f.endswith('.bed.gz') else None

def pick(inv, n, staged):
    sets = collections.defaultdict(set)
    for line in open(inv):
        path = line.rstrip('\n').split('\t')[1]
        d, f = os.path.split(path); t = parse(f)
        if t and '/KaryoScope/' in d and '_archive' not in d: sets[d].add(t)
    groups = collections.defaultdict(list)
    for d, fs in sets.items():
        db = next(iter(fs))[0]; v = 'v2' if db.endswith('_v2') else 'v1'; b = 'hks' if db.startswith('HKS') else 'kmc'
        lp = d.lower()
        c = ('asm' if '/assembly/' in lp else 'short' if 'illumina' in lp else
             'ont' if '/ont' in lp else 'hifi' if '/hifi' in lp else None)
        full = {(st, k) for st in STD - {'haplotype'} for k in ('presmoothed', 'smoothed')}
        if c and full <= {x[1:] for x in fs} and len({x[0] for x in fs}) == 1:
            smp = d.split('/KaryoScope/')[0].split('/')        # the path after its results directory
            smp = '/'.join(smp[max([i for i, x in enumerate(smp) if 'results' in x] or [len(smp) - 3]) + 1:])
            groups[(c, v, b)].append((smp, d))
    rnd = random.Random(20260910)
    for (c, v, b), ds in sorted(groups.items()):
        rnd.shuffle(ds); seen = set(); k = 0
        for smp, d in ds:
            top = smp.lower().split('/')                 # the sample: the path up to its data type
            top = '/'.join(top[:next((i for i, x in enumerate(top) if x.startswith(('assembly', 'illumina', 'ont', 'hifi'))), len(top))])
            if top in seen or k == n: continue
            seen.add(top); k += 1
            print(f"{c}_{smp.replace('/', '-').replace('.', '-')}.{v}.{b}\t{d}\t-")
    for cell in sorted(os.listdir(os.path.join(staged, 'kstracks'))):
        if not cell.startswith('chm13'):
            print(f"{cell}\t{staged}/kstracks/{cell}\t{staged}/ksaln/{cell.split('.')[0]}.paf")

def lines(path):
    return subprocess.Popen(['zcat', path], stdout=subprocess.PIPE, text=True, bufsize=1 << 20).stdout

def count(path):
    return int(subprocess.run(f"zcat '{path}' | wc -l", shell=True, capture_output=True, text=True).stdout)

def write(job):                                  # one track, filtered, bgzip level 6
    path, dest, bgz, frac, window, end, new = job
    p = subprocess.Popen([bgz, '-l', '6', '-c'], stdin=subprocess.PIPE, stdout=open(dest, 'wb'))
    k, prev, kept = 0, None, False
    for line in lines(path):
        nm, s, e, rest = line.split('\t', 3)
        if nm != prev: prev = nm; kept = window or frac == 1 or zlib.crc32(nm.encode()) < frac * 2 ** 32
        if not kept: continue
        if window:
            if int(s) >= end[nm]: continue
            e = str(min(int(e), end[nm]))
        p.stdin.write(f"{new[nm] if new else nm}\t{s}\t{e}\t{rest}".encode()); k += 1
    p.stdin.close()
    if p.wait(): sys.exit(f'sample: bgzip failed on {path}')
    return k

def cut(src, out, paf, cap, bgz):
    os.makedirs(out, exist_ok=True)
    smp = os.path.basename(out).split('.')[0]
    beds = sorted(f for f in os.listdir(src) if parse(f))
    with multiprocessing.Pool(8) as pool: rows = dict(zip(beds, pool.map(count, [os.path.join(src, f) for f in beds])))
    big = max(beds, key=rows.get)                 # the largest track sets the rate
    frac = min(1.0, cap / rows[big])
    nseq, per, prev = 0, collections.Counter(), None
    for line in lines(os.path.join(src, big)):
        nm = line.split('\t', 1)[0]; per[nm] += 1
        if nm != prev: nseq += 1; prev = nm
    window = frac < 1 and frac * nseq < 100
    end = {}
    if window:                                    # few sequences: each ends where the largest
        want = {nm: int(n * frac) for nm, n in per.items()}   # track reaches its share of rows
        seen = collections.Counter()
        for line in lines(os.path.join(src, big)):
            nm, s, e, _ = line.split('\t', 3); seen[nm] += 1
            if seen[nm] == want[nm] or (want[nm] == 0 and seen[nm] == 1): end[nm] = int(e) if want[nm] else 1
        open(os.path.join(out, 'window'), 'w').close()   # clipped ends: no smoothing oracle
    del per
    first = next(lines(os.path.join(src, big))).split('\t', 1)[0]
    new = names = None
    if first.isdigit() and not window:           # ordinal: read N is line N+1
        olds = sorted({int(l.split('\t', 1)[0]) for l in lines(os.path.join(src, big))
                       if zlib.crc32(l.split('\t', 1)[0].encode()) < frac * 2 ** 32})
        new = {str(o): str(j) for j, o in enumerate(olds)}
        allnames = open(os.path.join(src, 'names.txt')).read().split('\n')
        names = [allnames[o] for o in olds]
        open(os.path.join(out, 'names.txt'), 'w').write('\n'.join(names) + '\n')
    jobs = [(os.path.join(src, f), os.path.join(out, '{}.{}.{}.{}.bed.gz'.format(smp, *parse(f))), bgz, frac, window, end, new) for f in beds]
    with multiprocessing.Pool(8) as pool: kept = pool.map(write, jobs)
    q = 0
    if paf != '-' and not window:
        want = set(names) if names else None
        with open(os.path.join(out, 'aln.paf'), 'w') as o:
            for line in open(paf):
                nm = line.split('\t', 1)[0]
                if (nm in want) if want is not None else zlib.crc32(nm.encode()) < frac * 2 ** 32: o.write(line); q += 1
    print(f"sample: {os.path.basename(out)} {'window' if window else 'reads'} {frac:.6f} of up to "
          f"{rows[big]} rows, {len(beds)} tracks, {min(kept)} to {max(kept)} rows kept, {q} alignment lines", flush=True)

if sys.argv[1] == 'pick': pick(sys.argv[2], int(sys.argv[3]), sys.argv[4])
else: cut(sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), sys.argv[6])
