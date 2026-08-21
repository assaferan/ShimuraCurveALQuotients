#!/usr/bin/env python3
# Score closed-form candidates for the THETA-const weights kappa_d(D,N) against
# every base with ground truth. Odd N only (N = 2 is the delicate case, kept out).
#
#   mult_f  =  sum_{d | D} kappa_d * sum_{k >= 1} c_0(-N d k^2)
#
# Candidate shape:  2*kappa_d = prod_{p | D/d} u_p * prod_{p | d} v_p
# with u_p, v_p drawn from Legendre-symbol variants times a small multiplier m_p.
# We enumerate reciprocity splits and m_2 choices and report which variant fits
# every known multiplier on every odd-N base.

import glob, os, itertools
from fractions import Fraction as F

def parse(path):
    info, forms, coefs = None, {}, {}
    for line in open(path):
        t = line.split()
        if not t: continue
        if t[0] == 'BASEINFO': info = (int(t[1]), int(t[2]), int(t[3]), int(t[4]))
        elif t[0] == 'FORM':   forms[int(t[3])] = None if t[4] == '?' else F(t[4])
        elif t[0] == 'COEF':
            k = int(t[3]); idx = (t[4], int(t[5]))
            coefs.setdefault(k, {})[idx] = coefs.setdefault(k, {}).get(idx, F(0)) + F(t[7])
    return info, forms, coefs

def legendre(a, p):
    a %= p
    if a == 0: return 0
    r = pow(a, (p - 1) // 2, p)
    return 1 if r == 1 else -1

def square_free_split(n):
    k = 1; i = 2
    while i * i <= n:
        while n % (i * i) == 0:
            n //= i * i; k *= i
        i += 1
    return n, k

def prime_divisors(n):
    ps = []; i = 2
    while i * i <= n:
        if n % i == 0:
            ps.append(i)
            while n % i == 0: n //= i
        i += 1
    if n > 1: ps.append(n)
    return ps

def theta_indices(coefs_f, D, N, psi=None):
    """[(d, sum_k psi(k)*c)] for one form's oo principal part.
    d = squarefree part of m/N, ANY d (the d|D restriction was an anchor
    artifact -- 14_3's m=15 has d=5 with 5 ndiv 14 and weight 1/2).
    psi(k, N) twists the k-sum (10_7 refuted plain k-independence)."""
    agg = {}
    for (blk, m), c in coefs_f.items():
        if blk != 'oo' or m % N: continue
        d, k = square_free_split(m // N)
        w = c if psi is None else c * psi(k, N)
        agg[d] = agg.get(d, F(0)) + w
    return agg

def score(bases, kappa_fn, verbose=False, psi=None):
    lines, ok = [], True
    for (D, N, forms, coefs) in bases:
        for kf in sorted(forms):
            if forms[kf] is None: continue
            agg = theta_indices(coefs[kf], D, N, psi=psi)
            pred = sum(kappa_fn(D, N, d) * s for d, s in agg.items())
            hit = (pred == forms[kf])
            ok = ok and hit
            if verbose or not hit:
                lines.append(f"    X0^{D}({N}) form {kf}: pred {pred} vs true {forms[kf]}"
                             f" {'OK' if hit else '  <-- MISS'}")
    return ok, lines

def gtsweep_truth(here):
    """Parse gt_*.out gtsweep logs. Primary source: the per-form
    'form K : residual R closed form C TRUE multiplier T' lines (T is the truth,
    valid for BOTH parities). Fallback: the unique residual pair (odd N only,
    where the pipeline applies nothing and residual == multiplier)."""
    import re
    out = {}
    for p in sorted(glob.glob(os.path.join(here, 'gt_*.out'))):
        txt = open(p, errors='replace').read().replace('\n', ' ')
        mb = re.search(r'X0\^(\d+)\((\d+)\)', txt)
        if not mb: continue
        D, N = int(mb.group(1)), int(mb.group(2))
        forms = re.findall(r'form\s+(-\d+)\s*:\s*residual\s+(-?\d+(?:/\d+)?)\s+'
                           r'closed form\s+(-?\d+(?:/\d+)?)\s+TRUE multiplier\s+(-?\d+(?:/\d+)?)',
                           txt)
        if forms:
            out[(D, N)] = {int(k): F(t) for k, _, _, t in forms}
            continue
        mp = re.findall(r'admissible residual pairs[^{]*\{([^}]*)\}', txt)
        if not mp: continue
        pairs = re.findall(r'<\s*(-?\d+)\s*,\s*(-?\d+)\s*>', mp[0])
        if len(pairs) == 1 and N % 2 == 1:
            out[(D, N)] = {-1: F(pairs[0][0]), -2: F(pairs[0][1])}
    return out

def main():
    here = os.path.dirname(__file__)
    gtt = gtsweep_truth(here)
    if gtt:
        print("gtsweep truths found:", {f"X0^{D}({N})": (str(v[-1]), str(v[-2]))
                                        for (D, N), v in sorted(gtt.items())})
    bases = []
    seen = set()
    for p in sorted(glob.glob(os.path.join(here, 'pp*_*.out'))):
        try:
            (D, N, M, det), forms, coefs = parse(p)
        except Exception:
            continue
        if (D, N) in seen or N % 2 == 0: continue   # odd N only
        seen.add((D, N))
        if (D, N) in gtt:
            for row, val in gtt[(D, N)].items():
                if forms.get(row) is None: forms[row] = val
        bases.append((D, N, forms, coefs))
    bases.sort()
    print("odd-N bases:", ", ".join(f"X0^{D}({N})[{sum(1 for v in f.values() if v is not None)}]"
                                    for D, N, f, _ in bases))

    # Generalized shape (support = all N|m, weight through d = sqfree(m/N)):
    #   2*kappa_d = prod_{p | D, p nmid d} u_p  *  prod_{p | d} v_p
    # u_p in {(p/N), (N/p)};  v_p for odd p in {(N/p), (p/N)};
    # v_2 = m2 flat, or m2*(2/N), or m2*(N/2)_Kronecker.  m2 in {3, -3}.
    def kron2(N):  # Kronecker (2/N) for odd N
        return 1 if N % 8 in (1, 7) else -1
    results = []
    axes = itertools.product(('pN', 'Np'), ('Np', 'pN'),
                             ('flat', 'x2N', 'xN2'), (3, -3))
    for u_sym, v_sym, m2mode, m2 in axes:
        def kappa(D, N, d, u_sym=u_sym, v_sym=v_sym, m2mode=m2mode, m2=m2):
            val = F(1, 2)
            for p in prime_divisors(D):
                if d % p:
                    val *= legendre(p, N) if u_sym == 'pN' else legendre(N, p)
            for p in prime_divisors(d):
                if p == 2:
                    val *= m2 * (kron2(N) if m2mode in ('x2N', 'xN2') else 1)
                    # (2/N) and (N/2):=(2/N) agree as Kronecker symbols; keep one
                else:
                    val *= legendre(N, p) if v_sym == 'Np' else legendre(p, N)
            return val
        ok, lines = score(bases, kappa)
        results.append((ok, u_sym, v_sym, m2mode, m2))
        tag = "PASS ALL" if ok else f"{len(lines)} misses"
        print(f"  u={u_sym} v={v_sym} v2={m2mode}*{m2}: {tag}")
        if not ok and len(lines) <= 4:
            print("\n".join(lines))

    # ---- THE MU-LAW (found 2026-08-22 after 6_13 refuted the beta-symbol family):
    #      2*kappa_d = mu(d) * 3^[2|d] * Jacobi( D/gcd(D,d), N )
    def mobius(n):
        m = 1
        i = 2
        while i * i <= n:
            if n % i == 0:
                n //= i; m = -m
                if n % i == 0: return 0
            i += 1
        if n > 1: m = -m
        return m
    def jacobi(a, n):  # n odd > 0
        a %= n; r = 1
        while a:
            while a % 2 == 0:
                a //= 2
                if n % 8 in (3, 5): r = -r
            a, n = n, a
            if a % 4 == 3 and n % 4 == 3: r = -r
            a %= n
        return r if n == 1 else 0
    def kappa_mu(D, N, d):
        g = 1
        for p in prime_divisors(D):
            if d % p == 0: g *= p
        return F(mobius(d) * (3 if d % 2 == 0 else 1) * jacobi(D // g, N), 2)
    ok, lines = score(bases, kappa_mu, verbose=False)
    print(f"\n### MU-LAW  2k_d = mu(d)*3^[2|d]*(D/gcd(D,d) over N), psi = 1: "
          f"{'PASS ALL' if ok else 'MISSES'}")
    print("\n".join(lines))

    # ---- MU2-LAW (found after 10_7 -> (-3,-4) refuted mu-law + k-independence):
    #   w(N d k^2) = kappa_d * psi(k),  psi(k) = lambda(k)*(k/N)  [Liouville * Legendre]
    #   2*kappa_d = prod_{p | D, p ndiv d} g(p,N) * prod_{q | d} h(q,N)
    #   g(2,N) = -1,  g(p,N) = (N/p) odd p;  h(q) = (-1/q) odd q,
    #   h(2,N) = 3*(-1/N)  or  3*(-2/N)   [tied at N = 3, 5]
    def lam(k):  # Liouville
        v = 1; i = 2; n = k
        while i * i <= n:
            while n % i == 0: n //= i; v = -v
            i += 1
        if n > 1: v = -v
        return v
    def psi_mu2(k, N):
        # lambda(k)*(k/N); jacobi = 0 when gcd(k,N) > 1 (untested corner, flagged)
        return lam(k) * jacobi(k, N)
    for h2label, h2 in [("3*(-1/N)", lambda N: 3 * jacobi(N - 1, N)),
                        ("3*(-2/N)", lambda N: 3 * jacobi(N - 2, N))]:
        def kappa_mu2(D, N, d, h2=h2):
            val = F(1, 2)
            for p in prime_divisors(D):
                if d % p:
                    val *= -1 if p == 2 else legendre(N, p)
            for q in prime_divisors(d):
                if q == 2: val *= h2(N)
                else:      val *= (1 if q % 4 == 1 else -1)   # (-1/q)
            return val
        ok, lines = score(bases, kappa_mu2, verbose=True, psi=psi_mu2)
        print(f"\n### MU2-LAW  h(2,N) = {h2label}: {'PASS ALL' if ok else 'MISSES'}")
        print("\n".join(lines if not ok else lines))

    winners = [(u, v, mm, m2) for ok, u, v, mm, m2 in results if ok]
    print("\nold winning variants:", winners if winners else "NONE (refuted by 6_13)")
    for u, v, mm, m2 in winners[:1]:
        def kappa(D, N, d, u=u, v=v, mm=mm, m2=m2):
            val = F(1, 2)
            for p in prime_divisors(D):
                if d % p:
                    val *= legendre(p, N) if u == 'pN' else legendre(N, p)
            for p in prime_divisors(d):
                if p == 2:
                    val *= m2 * (kron2(N) if mm in ('x2N', 'xN2') else 1)
                else:
                    val *= legendre(N, p) if v == 'Np' else legendre(p, N)
            return val
        print(f"\nscorecard u={u} v={v} v2={mm}*{m2}:")
        ok, lines = score(bases, kappa, verbose=True)
        print("\n".join(lines))

if __name__ == '__main__':
    main()
