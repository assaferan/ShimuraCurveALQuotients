#!/usr/bin/env python3
# nmzsolve: replacement for the polymake step of get_integer_prog_solutions,
# powered by a standalone Normaliz binary. Enumerates ALL lattice points of the
# BorcherdsForms eta-quotient polytope and writes them in polymake_solution
# (Magma-evaluable) format.
#   nmzsolve.py <M> <n> <m> <outfile> [k24=12] [sq_disc=0] [cuspidal=0]
# k24 = 24*k (weight k): 12 for the weight-1/2 basis, 36 for weight 3/2, etc.
# Constraints (matching integer_programming_input in BorcherdsForms.m):
#   sum r = 2k;  sum v_2(d) r ≡ 1-sq_disc (2);  sum v_p(d) r ≡ 0 (2) for p>2;
#   sum d r ≡ 0 (24);  sum (M/d) r ≡ 0 (24);
#   ord24_b(r) >= 0 all cusps b|M, except b=1 >= -24m, b=M >= -24n
#   (cuspidal: ord24_b(r) >= 1 at every cusp).
import sys, os, math, subprocess, tempfile

NMZ = os.environ.get('NORMALIZ_BIN') or os.path.join(os.path.dirname(os.path.abspath(__file__)), 'normaliz-3.11.1', 'normaliz')
if not os.path.exists(NMZ):
    NMZ = 'normaliz'   # hope it is on PATH
NMZ_TIMEOUT = int(os.environ.get('NMZ_TIMEOUT', '1800'))

M, n_pole, m_pole, outp = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
k24 = int(sys.argv[5]) if len(sys.argv) > 5 else 12
sq_disc = int(sys.argv[6]) if len(sys.argv) > 6 else 0
cuspidal = int(sys.argv[7]) if len(sys.argv) > 7 else 0
assert k24 % 12 == 0
two_k = k24 // 12          # 2k, the sum of exponents

ds = [d for d in range(1, M+1) if M % d == 0]
nd = len(ds)
def vp(d, p):
    v = 0
    while d % p == 0: v += 1; d //= p
    return v
primes = sorted({p for p in range(2, M+1) if M % p == 0 and all(p % q for q in range(2, p))})
def ord24(delta, b):
    return M * math.gcd(b, delta)**2 // (math.gcd(M, b*b) * delta)

eqs = [[1]*nd + [-two_k]]
congs = []
congs.append([vp(d, 2) for d in ds] + [-(1 - sq_disc), 2])
for p in primes:
    if p == 2: continue
    congs.append([vp(d, p) for d in ds] + [0, 2])
congs.append([d % 24 for d in ds] + [0, 24])
congs.append([(M // d) % 24 for d in ds] + [0, 24])
ineqs = []
for b in ds:
    row = [ord24(d, b) for d in ds]
    if cuspidal:
        c = -1
    else:
        c = 24 * m_pole if b == 1 else (24 * n_pole if b == M else 0)
    ineqs.append(row + [c])

# ---- t-shift fallback for tall weight-1/2 {oo}-pole rungs ------------------
# Enumeration cost grows ~10x per +22 in n at 24-divisor levels, but rung n+k
# is spanned by rung n shifted by weight-0 pole<=k quotients (validated at 308:
# rank 108/108). If a weight-0 shift file and a lower enumerated rung exist in
# the cache, answer tall requests by shifting + thinning instead of normaliz.
import glob as _glob, re as _re
def _load_magma_seq(path):
    return [ [int(x) for x in m.split(',')] for m in _re.findall(r'\\?\[ ([-\d, ]+) \]', open(path).read().replace('\\\n', '').replace('\n', ' ')) ]
def _pole(r):
    s24 = sum(d*x for d, x in zip(ds, r))
    assert s24 % 24 == 0
    return -s24 // 24
if m_pole == 0 and k24 == 12 and cuspidal == 0:
    w0p = f'polymake/tshift_w0_{M}.txt'
    cands = []
    for f in _glob.glob(f'polymake/polymake_solution_{M}_*_0'):
        mm = _re.match(rf'polymake/polymake_solution_{M}_(\d+)_0$', f)
        if mm and int(mm.group(1)) < n_pole:
            cands.append(int(mm.group(1)))
    if os.path.exists(w0p) and cands:
        base_n = max(cands)
        base = _load_magma_seq(f'polymake/polymake_solution_{M}_{base_n}_0')
        shifts = [ t for t in _load_magma_seq(w0p) if any(t) ]
        # a rank-preserving CORE (same span as the full rung, ~dim vectors) beats
        # thinning by pole order: 42 of sol(420,145)'s 105 pivots come from
        # cancellations between same-leading-order vectors, which order-thinning
        # destroys.  If a core file exists, shift it verbatim (no thinning).
        corep = f'polymake/tshift_core_{M}.txt'
        use_core = os.path.exists(corep)
        if use_core:
            core = _load_magma_seq(corep)
            if core: base = core
        if base and shifts:
            tp = max(_pole(t) for t in shifts)                     # shift pole k
            # the cache accumulates rungs THIS fallback wrote (shifted, not
            # enumerated), so the filename's n overstates the reach of `base`
            # once it has been swapped for the core: compute J from the actual
            # max pole of the vectors being shifted, never from base_n.
            base_n = max(_pole(list(r)) for r in base)
            J = max(0, -(-(n_pole - base_n) // tp))                # ceil
            if use_core:
                pts = set()
                for j in range(J + 1):
                    for r in base:
                        v = tuple(a + j*b for a, b in zip(r, shifts[0])) if len(shifts) == 1 else None
                        if v is not None: pts.add(v)
                if len(shifts) > 1:
                    layer = [ tuple(r) for r in base ]
                    pts = set(layer)
                    for j in range(J):
                        layer = [ tuple(a + b for a, b in zip(r, t)) for t in shifts for r in layer ]
                        pts.update(layer)
                pts = sorted(pts)
            else:
                byorder = {}
                layer = [ tuple(r) for r in base ]
                seen = set(layer)
                def _absorb(vs):
                    for v in vs:
                        pl = _pole(list(v))
                        if pl > n_pole: continue
                        b = byorder.setdefault(pl, [])
                        if len(b) < 4:
                            b.append(v)
                _absorb(layer)
                for j in range(J):
                    nxt = []
                    for t in shifts:
                        for r in layer:
                            v = tuple(a + b for a, b in zip(r, t))
                            if v not in seen:
                                seen.add(v); nxt.append(v)
                    _absorb(nxt)
                    layer = nxt
                pts = sorted({ v for b in byorder.values() for v in b })
            residues = { _pole(list(r)) % tp for r in base }
            if len(residues) < tp:
                print(f"# tshift fallback WARNING: only {len(residues)}/{tp} pole residues covered", file=sys.stderr)
            print(f"# tshift fallback: base rung {base_n} ({len(base)} pts) + {len(shifts)} shift(s) x{J} -> {len(pts)} thinned pts for n={n_pole}", file=sys.stderr)
            with open(outp, 'w') as f:
                f.write('[ PowerSequence(IntegerRing()) |\n')
                f.write(',\n'.join('[ ' + ', '.join(str(x) for x in r) + ' ]' for r in pts))
                f.write('\n]\n')
            sys.exit(0)

work = tempfile.mkdtemp(prefix='nmz_')
base = os.path.join(work, f'poly_{M}_{n_pole}_{m_pole}')
with open(base + '.in', 'w') as f:
    f.write(f'amb_space {nd}\n')
    f.write(f'inhom_equations {len(eqs)}\n')
    for r in eqs: f.write(' '.join(map(str, r)) + '\n')
    f.write(f'inhom_congruences {len(congs)}\n')
    for r in congs: f.write(' '.join(map(str, r)) + '\n')
    f.write(f'inhom_inequalities {len(ineqs)}\n')
    for r in ineqs: f.write(' '.join(map(str, r)) + '\n')
    f.write('LatticePoints\n')

try:
    res = subprocess.run([NMZ, '-x=8', base + '.in'], capture_output=True, text=True, timeout=NMZ_TIMEOUT)
except subprocess.TimeoutExpired:
    sys.stderr.write(f"# normaliz timed out after {NMZ_TIMEOUT}s on {base}.in\n")
    sys.exit(2)
if res.returncode != 0:
    sys.stderr.write(res.stdout[-2000:] + res.stderr[-2000:])
    sys.exit(1)

pts = []
with open(base + '.out') as f:
    lines = f.readlines()
i = 0
while i < len(lines):
    if 'lattice points in polytope' in lines[i] and lines[i].rstrip().endswith(':'):
        cnt = int(lines[i].split()[0])
        i += 1
        got = 0
        while got < cnt and i < len(lines):
            parts = lines[i].split()
            if len(parts) == nd + 1 and all(p.lstrip('-').isdigit() for p in parts):
                assert parts[-1] == '1', parts
                pts.append(tuple(int(x) for x in parts[:nd]))
                got += 1
            i += 1
        break
    i += 1

pts = sorted(set(pts))
print(f"# M={M} n={n_pole} m={m_pole} k24={k24} sq_disc={sq_disc} cuspidal={cuspidal}: {len(pts)} lattice points", file=sys.stderr)
with open(outp, 'w') as f:
    f.write('[ PowerSequence(IntegerRing()) |\n')
    f.write(',\n'.join('[ ' + ', '.join(str(x) for x in r) + ' ]' for r in pts))
    f.write('\n]\n')
