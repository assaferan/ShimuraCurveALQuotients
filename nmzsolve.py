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

def _write_solution(outp, pts, derived=None):
    """Write the Magma-evaluable solution file ATOMICALLY.

    ⚠ WHY THIS IS NOT JUST `open(outp,'w')`. BorcherdsForms.m:180 caches every solve at a
    DETERMINISTIC SHARED PATH, polymake/polymake_solution_<M>_<n>_<m>, and reads it back as
    FileExists(...) -> eval Read(...). A plain open() truncates the file and then fills it, so a
    SECOND Magma process running concurrently and needing the same uncached triple can see the
    path exist and read a PARTIAL point list -- valid Magma, fewer points, no error anywhere.
    That is the "a partially-cached base returns a wrong answer rather than an error" failure
    CLAUDE.md flags as critical, and parallel test runs make it reachable.

    os.replace is atomic on POSIX, so a reader sees either the old complete file or the new
    complete one, never a half-written one. Two writers racing is harmless: they compute the same
    answer, so whichever rename lands last wins with identical content. The temp file carries the
    pid so concurrent writers do not clobber each other's scratch.
    """
    tmp = f'{outp}.tmp.{os.getpid()}'
    with open(tmp, 'w') as f:
        if derived:
            # A Magma comment, so `eval Read(file)` still returns the sequence; marks a file that
            # holds a spanning subset built from a lower rung, not the lattice points of the polytope.
            f.write(f'// derived: {derived}\n')
        f.write('[ PowerSequence(IntegerRing()) |\n')
        f.write(',\n'.join('[ ' + ', '.join(str(x) for x in r) + ' ]' for r in pts))
        f.write('\n]\n')
    os.replace(tmp, outp)


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

import glob as _glob, re as _re
def _load_magma_seq(path):
    return [ [int(x) for x in m.split(',')] for m in _re.findall(r'\\?\[ ([-\d, ]+) \]', open(path).read().replace('\\\n', '').replace('\n', ' ')) ]
def _pole(r):
    s24 = sum(d*x for d, x in zip(ds, r))
    assert s24 % 24 == 0
    return -s24 // 24
# ---- a lower two-cusp rung is a subset of a higher cached one -----------------
# (M, n, m) and (M, n', m) differ only in the bound on the pole at oo, so for n' > n the lower rung
# is the points of the higher one whose pole at oo is <= n.  Three conditions: the file name records
# only the level and the pole orders, so this serves only the weight-1/2, default-parity,
# non-cuspidal request; only rungs with a pole at 0 (m > 0), since shifted sets are written under
# the m = 0 names; and a cached file marked derived (a spanning set, not the lattice points) is skipped.
def _is_derived(path):
    with open(path) as f:
        return f.readline().startswith('//')
if m_pole > 0 and k24 == 12 and sq_disc == 0 and cuspidal == 0:
    higher = []
    for f in _glob.glob(f'polymake/polymake_solution_{M}_*_{m_pole}'):
        mm = _re.match(rf'polymake/polymake_solution_{M}_(\d+)_{m_pole}$', f)
        if mm and int(mm.group(1)) > n_pole and not _is_derived(f):
            higher.append(int(mm.group(1)))
    if higher:
        hi = min(higher)
        pts = sorted({ tuple(r) for r in _load_magma_seq(f'polymake/polymake_solution_{M}_{hi}_{m_pole}') if _pole(r) <= n_pole })
        print(f"# lower rung: (M={M}, n={n_pole}, m={m_pole}) read off the cached (M={M}, n={hi}, m={m_pole}): {len(pts)} lattice points", file=sys.stderr)
        _write_solution(outp, pts)
        sys.exit(0)

# ---- shifted spanning sets for tall rungs (polymake/README.md) ------------------
# The products of a lower ENUMERATED rung at the same cusp-0 order with the level's weight-0 shift
# set span the forms of a tall rung; the file written is marked derived (a spanning set, not the
# lattice points), is never used as a base, and is enumerated by the caller when its rank falls
# short.  NMZ_NO_SHIFT=1 bypasses this section; the read of a lower rung off a higher one runs first.
if k24 == 12 and cuspidal == 0 and os.environ.get('NMZ_NO_SHIFT', '') != '1':
    w0p = f'polymake/tshift_w0_{M}.txt'
    cands = []
    for f in _glob.glob(f'polymake/polymake_solution_{M}_*_{m_pole}'):
        mm = _re.match(rf'polymake/polymake_solution_{M}_(\d+)_{m_pole}$', f)
        if mm and int(mm.group(1)) < n_pole and not _is_derived(f):
            cands.append(int(mm.group(1)))
    if os.path.exists(w0p) and cands:
        base_n = max(cands)
        base = _load_magma_seq(f'polymake/polymake_solution_{M}_{base_n}_{m_pole}')
        shifts = [ t for t in _load_magma_seq(w0p) if any(t) ]
        # a rank-preserving CORE (same span as the full rung, ~dim vectors) beats
        # thinning by pole order: 42 of sol(420,145)'s 105 pivots come from
        # cancellations between same-leading-order vectors, which order-thinning
        # destroys.  If a core file exists, shift it verbatim (no thinning).
        corep = f'polymake/tshift_core_{M}.txt'
        use_core = m_pole == 0 and os.path.exists(corep)      # the core is a reduced m = 0 rung
        if use_core:
            core = _load_magma_seq(corep)
            if core: base = core
        if base and shifts:
            tp = max(_pole(t) for t in shifts)                     # shift pole k
            # the core is a reduced rung, so its reach is the largest pole among its vectors
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
                # Every shifted point with pole <= n is kept while the set has at most NMZ_SHIFT_KEEP
                # points (default 20000); thinning to a few points per pole order can lose rank.
                keep_all = int(os.environ.get('NMZ_SHIFT_KEEP', '20000'))
                byorder = {}
                layer = [ tuple(r) for r in base ]
                seen = set(layer)
                def _absorb(vs):
                    for v in vs:
                        pl = _pole(list(v))
                        if pl > n_pole: continue
                        byorder.setdefault(pl, []).append(v)
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
                total = sum(len(b) for b in byorder.values())
                if total > keep_all:
                    # beyond the cap, keep as many points per pole order as the cap allows, at least four
                    per = max(4, keep_all // max(1, len(byorder)))
                    print(f"# tshift fallback: {total} shifted points, keeping {per} per pole order", file=sys.stderr)
                    pts = sorted({ v for b in byorder.values() for v in b[:per] })
                else:
                    pts = sorted({ v for b in byorder.values() for v in b })
            if len(shifts) == 1:
                # with a single shift of pole k, only the residues mod k present in the base are reached
                residues = { _pole(list(r)) % tp for r in base }
                if len(residues) < tp:
                    print(f"# tshift fallback WARNING: only {len(residues)}/{tp} pole residues covered", file=sys.stderr)
            print(f"# tshift fallback: base rung {base_n} ({len(base)} pts) + {len(shifts)} shift(s) x{J} -> {len(pts)} pts for n={n_pole}", file=sys.stderr)
            _write_solution(outp, pts, derived=f'spanning set for (M={M}, n={n_pole}, m={m_pole}) from the rung {base_n} and {len(shifts)} weight-0 shift(s), not the lattice points of the polytope')
            sys.exit(0)

# ---- coordinate order ---------------------------------------------------------
# Normaliz lifts the lattice points coordinate by coordinate in the order of the input.  The "ends"
# order puts the divisors farthest from sqrt(M) first (1, M, 2, M/2, ...); the solution file is
# written back in divisor order, so nothing downstream sees the permutation.  NMZ_ORDER=id keeps
# the divisor order.
order = os.environ.get('NMZ_ORDER', 'ends')
if order == 'ends':
    perm = sorted(range(nd), key=lambda i: -abs(math.log(ds[i]) - math.log(M)/2))
else:
    perm = list(range(nd))
def permute(row):
    return [row[i] for i in perm] + row[nd:]

work = tempfile.mkdtemp(prefix='nmz_')
base = os.path.join(work, f'poly_{M}_{n_pole}_{m_pole}')
with open(base + '.in', 'w') as f:
    f.write(f'amb_space {nd}\n')
    f.write(f'inhom_equations {len(eqs)}\n')
    for r in eqs: f.write(' '.join(map(str, permute(r))) + '\n')
    f.write(f'inhom_congruences {len(congs)}\n')
    for r in congs: f.write(' '.join(map(str, permute(r))) + '\n')
    f.write(f'inhom_inequalities {len(ineqs)}\n')
    for r in ineqs: f.write(' '.join(map(str, permute(r))) + '\n')
    f.write('LatticePoints\n')

try:
    # NMZ_THREADS (default 8): the lifting phase is parallel over partial points
    res = subprocess.run([NMZ, f"-x={os.environ.get('NMZ_THREADS', '8')}", '-k', base + '.in'], capture_output=True, text=True, timeout=NMZ_TIMEOUT)
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
cnt = None
while i < len(lines):
    if 'lattice points in polytope' in lines[i] and lines[i].rstrip().endswith(':'):
        cnt = int(lines[i].split()[0])
        i += 1
        got = 0
        while got < cnt and i < len(lines):
            parts = lines[i].split()
            if len(parts) == nd + 1 and all(p.lstrip('-').isdigit() for p in parts):
                assert parts[-1] == '1', parts
                v = [int(x) for x in parts[:nd]]
                r = [0]*nd
                for pos, j in enumerate(perm): r[j] = v[pos]     # back to divisor order
                pts.append(tuple(r))
                got += 1
            i += 1
        break
    i += 1

# Normaliz can exit 0 and still leave no usable answer: an output file without the
# lattice-points header, or one cut short of the count it announces.  Writing an empty or
# short solution file there would be taken downstream as "the polytope has few points" and
# produce an incomplete eta-quotient basis, so both are errors and no file is written.
if cnt is None:
    sys.stderr.write(f"# normaliz output {base}.out has no 'lattice points in polytope' header\n")
    sys.exit(3)
if len(pts) != cnt:
    sys.stderr.write(f"# normaliz output {base}.out announces {cnt} lattice points but lists {len(pts)}\n")
    sys.exit(4)

pts = sorted(set(pts))
print(f"# M={M} n={n_pole} m={m_pole} k24={k24} sq_disc={sq_disc} cuspidal={cuspidal}: {len(pts)} lattice points", file=sys.stderr)
_write_solution(outp, pts)
