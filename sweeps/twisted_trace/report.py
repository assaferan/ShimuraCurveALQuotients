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
Coverage: each RES line says which q it tested (legacy lines without qmax/pmax/pmin: p <= 59,
q <= 59^2).  Per curve the union over all files is compared with the Weil range {q = p^v < 4g^2, p | DN
excluded}; U/R curves with no counted violation and a gap go to curves_supplement.txt, and their
levels, with the smallest missing prime as PMIN, to levels_supplement.txt ("D N PMIN V3ONLY", V3ONLY = 0).
V3 coverage is tracked separately.  When 9 in W, V3 ops are Q-rational (sigma(V3) = V3 W9) and twist.m
now uses them at every good p (RES lines marked "v3all 1").  Older lines (no v3all) used V3 ops only
at p = 1 mod 3, so a curve with 9 in W and a V3 op in its list is complete only once V3 ops have been
tested on the whole Weil range: a line contributes to V3 coverage its full range if v3all, its
p = 1 mod 3 part if not, and a V3ONLY line (v3only 1: only V3 ops, only p = 2 mod 3) only to V3
coverage.  Curves with a V3 gap at p = 2 mod 3 get a levels_supplement.txt row with V3ONLY = 1 and
PMIN = the smallest missing such prime (run.sh -> twist.m V3ONLY:=1 -> out/D_N.v3.out).
Sanity: every recorded violation must satisfy |tr| <= 2g sqrt(q) and q < 4g^2 (Weil); else WEIL-BREACH.
Writes final_ruled_out.txt: CurveID D N g W q h tr  (min-q counted violation), U and R curves.
"""
import collections, glob, math, os
from fractions import Fraction

T = os.path.dirname(os.path.abspath(__file__)) + '/'
GL = [2555, 2568, 4190, 5635, 5639, 6616, 8495, 7926, 7932]   # group-lemma curves
P = print

PR = [p for p in range(2, 4000) if all(p % d for d in range(2, int(p ** .5) + 1))]

def qset(L, qmax, pmin=0, pmax=10 ** 9):
    """All q = p^v <= qmax with p prime, p not dividing L, pmin <= p <= pmax."""
    out = set()
    for p in PR:
        if p > min(qmax, pmax): break
        if p < pmin or L % p == 0: continue
        q = p
        while q <= qmax: out.add(q); q *= p
    return out

def need(g, L): return qset(L, 4 * g * g - 1)      # Weil range: a violation needs q < 4g^2

def pdiv(q): return next(p for p in PR if q % p == 0)

def is_al(h): return h == '1' or (h.startswith('w') and h[1:].isdigit())

curves = {}
for l in open(T + 'curves_all.txt'):
    f = l.split()
    if len(f) >= 7 and not l.startswith('#'):
        curves[int(f[1])] = dict(st=f[0], D=int(f[2]), N=int(f[3]), g=int(f[4]), W=f[6])

files = sorted(glob.glob(T + 'results_*.out')) + sorted(glob.glob(T + 'out/*.out'))
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
        tail = {f[i]: int(f[i + 1]) for i in range(21, len(f) - 1, 2)}
        if 'qmax' in tail:   # tested exactly the good primes pmin <= p <= pmax, powers up to qmax
            cov = qset(D * N, tail['qmax'], tail.get('pmin', 0), tail['pmax'])
        else:                # legacy local-run line: PB = 59
            cov = qset(D * N, 59 ** 2, 0, 59)
        v3only = tail.get('v3only', 0) == 1
        if v3only:                       # only V3 ops, only p = 2 mod 3: V3 coverage only
            v3cov = {q for q in cov if pdiv(q) % 3 == 2}; cov = set()
        elif tail.get('v3all', 0) == 1:  # V3 ops at every p in range
            v3cov = set(cov)
        else:                            # pre-v3all rule: V3 ops at p = 1 mod 3 only
            v3cov = {q for q in cov if pdiv(q) % 3 == 1}
        hasv3 = 9 in Ws and any('V3' in h for h in f[20].split(','))
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
            r['cov'] |= cov; r['v3cov'] |= v3cov; r['hasv3'] |= hasv3
            r['src'].add(os.path.basename(fn))
        else:
            recs[cid] = dict(st=st, D=D, N=N, g=g, W=W, good=good, noncomm=noncomm,
                             cov=cov, v3cov=v3cov, hasv3=hasv3, src={os.path.basename(fn)})

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

gap = {i: sorted(need(r['g'], r['D'] * r['N']) - r['cov']) for i, r in recs.items()
       if r['st'] in 'UR' and not r['good']}
gap = {i: m for i, m in gap.items() if m}
# V3 gap: q < 4g^2 not yet tested with the V3 ops (9 in W, V3 op in the list).  Its p = 1 mod 3 part and
# its p >= 61 part lie inside the general gap (every line covers V3 at p = 1 mod 3; a new general
# supplement is v3all), so only the p = 2 mod 3 part needs a V3ONLY run.
v3gap = {i: sorted(q for q in need(r['g'], r['D'] * r['N']) - r['v3cov'] if pdiv(q) % 3 == 2)
         for i, r in recs.items() if r['st'] in 'UR' and not r['good'] and r['hasv3']}
v3gap = {i: m for i, m in v3gap.items() if m}
lvgap = collections.defaultdict(list)
for i, m in gap.items(): lvgap[(recs[i]['D'], recs[i]['N'])].append(i)
lvv3 = collections.defaultdict(list)
for i, m in v3gap.items(): lvv3[(recs[i]['D'], recs[i]['N'])].append(i)
inc = set(gap) | set(v3gap)
P('\n=== 6. INCOMPLETE COVERAGE: %d U/R curves (no violation) with untested q < 4g^2, at %d levels%s'
  % (len(inc), len(set(lvgap) | set(lvv3)), '' if inc else '  -- all tested curves exhausted'))
P('  all ops (general gap): %d curves at %d levels; V3 ops at p = 2 mod 3 (9 in W): %d curves (%s) at %d levels'
  % (len(gap), len(lvgap), len(v3gap), dict(collections.Counter(recs[i]['st'] for i in v3gap)), len(lvv3)))
if gap:
    P('  missing q by genus:', {g: sorted({q for i, m in gap.items() if recs[i]['g'] == g for q in m})
                                for g in sorted({recs[i]['g'] for i in gap})})
    P('  all missing q prime (no prime powers):', all(len([p for p in PR if q % p == 0]) == 1 and q in PR
                                                       for m in gap.values() for q in m))
with open(T + 'curves_supplement.txt', 'w') as F:
    for l in open(T + 'curves_all.txt'):
        f = l.split()
        if len(f) >= 7 and int(f[1]) in inc: F.write(l)
with open(T + 'levels_supplement.txt', 'w') as F:     # D N PMIN V3ONLY
    rows = [(lv, min(pdiv(gap[i][0]) for i in lvgap[lv]), 0) for lv in lvgap] + \
           [(lv, min(pdiv(v3gap[i][0]) for i in lvv3[lv]), 1) for lv in lvv3]
    for lv, pm, v in sorted(rows, key=lambda x: (x[0][0] * x[0][1], x[0], x[2])):
        F.write('%d %d %d %d\n' % (lv[0], lv[1], pm, v))
P('wrote curves_supplement.txt, levels_supplement.txt')

with open(T + 'final_ruled_out.txt', 'w') as F:
    F.write('# CurveID D N g W q h tr   (twisted trace, min-q counted violation; U and R curves)\n')
    for i, r in sorted((i, r) for i, r in recs.items() if r['st'] in 'UR' and r['good']):
        q, h, tr = min(r['good'])
        F.write('%d %d %d %d %s %d %s %s\n' % (i, r['D'], r['N'], r['g'], r['W'], q, h, tr))
P('\nwrote', T + 'final_ruled_out.txt')
