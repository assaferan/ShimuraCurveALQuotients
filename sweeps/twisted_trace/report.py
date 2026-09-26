#!/usr/bin/env python3
"""Summarise the twisted-trace sweep: results_local.out (the 2026-09-25 local run, combined file)
plus every per-level out/*.out, deduplicated by CurveID.

A violation is COUNTED iff
  tr < -(q+1)                     (tr > q+1 is impossible for any curve: reported as BUG, never counted)
  h involves V3 only if 9 in W    (otherwise V3 need not be defined over Q on C)
  noncomm = 0 on that RES line, unless h is 1 or an AL wQ (twist.m already skips non-commuting ops
                                  per prime; a nonzero noncomm drops that curve's non-AL violations
                                  because the line does not say which op failed to commute)
A curve seen in several files keeps the union of its counted violations (every file is a sound run).
Sanity: every recorded violation must satisfy |tr| <= 2g sqrt(q) and q < 4g^2 (Weil); else WEIL-BREACH.
Writes final_ruled_out.txt: CurveID D N g W q h tr  (min-q counted violation), U and R curves.
"""
import collections, glob, math, os
from fractions import Fraction

T = os.path.dirname(os.path.abspath(__file__)) + '/'
GL = [2555, 2568, 4190, 5635, 5639, 6616, 8495, 7926, 7932]   # group-lemma curves
P = print

def maxgoodprime(n, L):
    return max([p for p in range(2, n + 1) if L % p and all(p % d for d in range(2, int(p ** .5) + 1))], default=0)

def is_al(h): return h == '1' or (h.startswith('w') and h[1:].isdigit())

curves = {}
for l in open(T + 'curves_all.txt'):
    f = l.split()
    if len(f) >= 7 and not l.startswith('#'):
        curves[int(f[1])] = dict(st=f[0], D=int(f[2]), N=int(f[3]), g=int(f[4]), W=f[6])

files = [T + 'results_local.out'] + sorted(glob.glob(T + 'out/*.out'))
recs, bugs, breach, baddim = {}, [], [], []
ndup = 0
levels, done = set(), set()
for fn in files:
    for l in open(fn):
        f = l.split()
        if not f: continue
        if f[0] == 'LEVEL': levels.add((int(f[1]), int(f[2])))
        if f[0] == 'DONE': done.add((int(f[1]), int(f[2])))
        if f[0] == 'BADDIM': baddim.append((os.path.basename(fn), l.strip()))
        if f[0] != 'RES' or len(f) < 21 or f[19] != 'ops': continue       # skip truncated lines
        kv = {f[i]: f[i + 1] for i in range(7, 19, 2)}
        st, cid, D, N, g, W = f[1], int(f[2]), int(f[3]), int(f[4]), int(f[5]), f[6]
        Ws = set(int(x) for x in W.split(','))
        noncomm = int(kv['noncomm'])
        has = len(f) > 24 and f[21] == 'qmax' and f[23] == 'pmax'
        qmax, pmax = (int(f[22]), int(f[24])) if has else (59 ** 2, 59)      # legacy lines: PB = 59
        exh = qmax >= 4 * g * g - 1 and pmax >= maxgoodprime(4 * g * g - 1, D * N)   # Weil range covered
        raw = [] if kv['viol'] in ('-', '[]') else [x.split(':') for x in kv['viol'].rstrip(';').split(';')]
        good = set()
        for q, h, tr in raw:
            q, tr = int(q), Fraction(tr)
            if abs(tr) > 2 * g * math.sqrt(q) + 1e-9 or q >= 4 * g * g: breach.append((cid, q, h, tr))
            if tr > q + 1: bugs.append((cid, q, h, tr)); continue
            if tr >= -(q + 1): continue
            if 'V3' in h and 9 not in Ws: continue
            if noncomm > 0 and not is_al(h): continue
            good.add((q, h, tr))
        if cid in recs:
            ndup += 1
            r = recs[cid]; r['good'] |= good; r['noncomm'] = max(r['noncomm'], noncomm)
            r['weil'] = r['weil'] or exh
            r['src'].add(os.path.basename(fn))
        else:
            recs[cid] = dict(st=st, D=D, N=N, g=g, W=W, good=good, noncomm=noncomm,
                             weil=exh, src={os.path.basename(fn)})

P('=== inputs: %d files (results_local.out + %d out/*.out); %d RES curves %s; %d duplicate RES merged'
  % (len(files), len(files) - 1, len(recs), dict(collections.Counter(r['st'] for r in recs.values())), ndup))
P('levels with DONE:', len(done), '| curves with noncomm>0:', sum(r['noncomm'] > 0 for r in recs.values()))
P('BUGS (tr > q+1):', bugs or 'none')
P('WEIL-BREACH (|tr| > 2g sqrt q or q >= 4g^2):', breach or 'none')
P('BADDIM:', baddim or 'none')

H = {i: r for i, r in recs.items() if r['st'] == 'H'}
Hf = {i: r for i, r in H.items() if r['good']}
P('\n=== 1. Controls (H): %d tested, %d with a violation' % (len(H), len(Hf)))
if Hf:
    P('!!!!!!!! CONTROL FAILURE -- TEST IS UNSOUND OR BUGGY !!!!!!!!')
    for i, r in sorted(Hf.items()): P('  ', i, r['D'], r['N'], r['W'], sorted(r['good']))

def split(d): return 'D>1 %d / D=1 %d' % (sum(r['D'] > 1 for r in d.values()), sum(r['D'] == 1 for r in d.values()))
nU = sum(c['st'] == 'U' for c in curves.values())
U = {i: r for i, r in recs.items() if r['st'] == 'U'}
Uf = {i: r for i, r in U.items() if r['good']}
P('\n=== 2. Undecided (U): tested %d of %d (%s); RULED OUT %d (%s)' % (len(U), nU, split(U), len(Uf), split(Uf)))
h1only = [i for i, r in Uf.items() if all(h == '1' for _, h, _ in r['good'])]
nonal = [i for i, r in Uf.items() if all(not is_al(h) for _, h, _ in r['good'])]
P('  ruled out by h = 1 only (plain trace):', len(h1only))
P('  needing a non-AL h (no AL/1 violation):', len(nonal), sorted(nonal))
P('  ruled out by genus:', dict(sorted(collections.Counter(r['g'] for r in Uf.values()).items())))
nw = [i for i, r in U.items() if not r['good'] and not r['weil'] and r['g'] >= 4]
P('  not ruled out, g >= 4, Weil range q < 4g^2 not exhausted (p <= 59 only):', len(nw),
  'at', len({(U[i]['D'], U[i]['N']) for i in nw}), 'levels')

P('\n=== 3. Group-lemma curves')
for i in GL:
    r = recs.get(i)
    P('  %d %s' % (i, 'not tested' if r is None else '%s (%d,%d) W={%s}: %s' % (
        r['st'], r['D'], r['N'], r['W'], 'RULED OUT %s' % (min(r['good']),) if r['good'] else 'no violation')))
ov = [i for i in GL if i in Uf]
P('  overlap %d; NEW beyond group lemma: %d' % (len(ov), len(Uf) - len(ov)))

P('\n=== 4. Reopened (R)')
for i in sorted(i for i, c in curves.items() if c['st'] == 'R'):
    r = recs.get(i)
    P('  %d %s' % (i, 'not done' if r is None else '(%d,%d) g=%d W={%s}: %s' % (
        r['D'], r['N'], r['g'], r['W'], ('RULED OUT ' + ', '.join('q=%d h=%s tr=%s' % x for x in sorted(r['good'])))
        if r['good'] else 'no violation')))

pend = collections.defaultdict(list)
for i, c in curves.items():
    if c['st'] in 'UR' and i not in recs: pend[(c['D'], c['N'])].append(i)
P('\n=== 5. Pending: %d U/R curves at %d levels' % (sum(map(len, pend.values())), len(pend)))
for lv in sorted(pend, key=lambda x: (x[0] * x[1], x)):
    o = T + 'out/%d_%d.out' % lv
    tag = 'partial/killed/running' if os.path.exists(o) else 'not started'
    P('  (%d,%d) DN=%d  %d curves  %s' % (lv[0], lv[1], lv[0] * lv[1], len(pend[lv]), tag))

with open(T + 'final_ruled_out.txt', 'w') as F:
    F.write('# CurveID D N g W q h tr   (twisted trace, min-q counted violation; U and R curves)\n')
    for i, r in sorted((i, r) for i, r in recs.items() if r['st'] in 'UR' and r['good']):
        q, h, tr = min(r['good'])
        F.write('%d %d %d %d %s %d %s %s\n' % (i, r['D'], r['N'], r['g'], r['W'], q, h, tr))
P('\nwrote', T + 'final_ruled_out.txt')
