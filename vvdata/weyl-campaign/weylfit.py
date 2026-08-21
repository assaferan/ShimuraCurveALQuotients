#!/usr/bin/env python3
# The Weyl-vector shortcut (handoff 2026-08-21 evening, next move 1).
#
# Hypothesis under test: the measured m=0 multiplier is a THETA-PAIRING functional
# of the principal part -- the Bruinier-Funke shadow shape
#       mult_f = sum_{d | D} kappa_d * sum_{k>=1} c_0(-N d k^2)        [THETA-const]
# i.e. weight kappa_d at every index m = N*d*k^2, INDEPENDENT of k (a unary theta
# q^{t k^2} has k-independent coefficients), zero elsewhere (incl. the whole 0-block).
#
# Controls, expected to fail:
#   THETA-klin : weight kappa_d * k          (the sigma_1-flavoured variant)
#   SIGMA      : kappa * sum_{e | n} chi(e) e^a,  n = m/N   (per-base constant)
#
# Data: pp2_/pp3_ dumps (COEF/FORM lines in the code's own indexing).
# Everything exact over Q.

import sys, glob, os
from fractions import Fraction as F

def parse(path):
    info, forms, coefs = None, {}, {}
    for line in open(path):
        t = line.split()
        if not t: continue
        if t[0] == 'BASEINFO':
            info = (int(t[1]), int(t[2]), int(t[3]), int(t[4]))
        elif t[0] == 'FORM':
            forms[int(t[3])] = None if t[4] == '?' else F(t[4])
        elif t[0] == 'COEF':
            k = int(t[3]); idx = (t[4], int(t[5]))
            coefs.setdefault(k, {})[idx] = coefs.setdefault(k, {}).get(idx, F(0)) + F(t[7])
    return info, forms, coefs

def rref(M, rhs):
    M = [row[:] for row in M]; rhs = rhs[:]
    rows = len(M); cols = len(M[0]) if rows else 0
    piv = []; r = 0
    for c in range(cols):
        p = next((i for i in range(r, rows) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]; rhs[r], rhs[p] = rhs[p], rhs[r]
        pv = M[r][c]; M[r] = [x / pv for x in M[r]]; rhs[r] = rhs[r] / pv
        for i in range(rows):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
                rhs[i] = rhs[i] - f * rhs[r]
        piv.append(c); r += 1
        if r == rows: break
    return M, rhs, piv

def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]

def square_free_split(n):
    # n = d * k^2 with d squarefree (k maximal)
    k = 1; i = 2
    while i * i <= n:
        while n % (i * i) == 0:
            n //= i * i; k *= i
        i += 1
    return n, k

def theta_decomp(m, D, N):
    # m = N * d * k^2 with d | D  ->  (d, k), else None. (D squarefree here.)
    if m % N: return None
    d, k = square_free_split(m // N)
    return (d, k) if D % d == 0 else None

class Base:
    def __init__(self, path):
        (self.D, self.N, self.M, self.det), self.forms, self.coefs = parse(path)
    def key(self): return (self.D, self.N)
    def known(self): return [k for k in sorted(self.forms) if self.forms[k] is not None]

def solve_family(bases, weightfn, unknown_of_idx, family):
    """Per-base solve: mult_f = sum_idx c_f(idx) * w(idx), w parameterized by unknowns.
    unknown_of_idx(base, blk, m) -> (name, coeff_multiplier) or 'ZERO' (weight forced 0)
    Returns report lines."""
    out = [f"### family {family}"]
    all_sol = {}
    ok_all = True
    for b in bases:
        names = []
        rows, rhs = [], []
        clash = []   # indices forced to weight 0 that carry nonzero coeff on a known form
        for kf in b.known():
            row = {}
            for (blk, m), c in b.coefs[kf].items():
                u = unknown_of_idx(b, blk, m)
                if u is None:  # weight forced 0 (unsupported)
                    continue
                name, wmul = u
                row[name] = row.get(name, F(0)) + c * wmul
            for nm in row:
                if nm not in names: names.append(nm)
            rows.append(row); rhs.append(b.forms[kf])
        names.sort(key=str)
        A = [[row.get(nm, F(0)) for nm in names] for row in rows]
        R, r2, piv = rref(A, rhs)
        incons = [i for i in range(len(R)) if all(x == 0 for x in R[i]) and r2[i] != 0]
        nfree = len(names) - len(piv)
        spare = len(rows) - len(piv)
        tag = "INCONSISTENT" if incons else "consistent"
        if incons: ok_all = False
        out.append(f"  X0^{b.D}({b.N}): eqs {len(rows)}, unknowns {len(names)}, "
                   f"rank {len(piv)}, free {nfree}, spare {spare}  -> {tag}")
        if not incons:
            sol = {}
            free = [c for c in range(len(names)) if c not in piv]
            for r_i, c_i in enumerate(piv): sol[names[c_i]] = r2[r_i]
            for c_i in free: sol[names[c_i]] = None
            all_sol[b.key()] = sol
            pretty = ", ".join(f"{nm}={'FREE' if v is None else v}" for nm, v in sol.items())
            out.append(f"      {pretty}")
    return out, all_sol, ok_all

def main():
    paths = sorted(glob.glob(os.path.join(os.path.dirname(__file__), 'pp*_*.out')))
    bases = []
    seen = set()
    for p in paths:
        try:
            b = Base(p)
        except Exception:
            continue
        if b.key() in seen: continue
        seen.add(b.key()); bases.append(b)
    bases.sort(key=lambda b: (b.D, b.N))
    print("bases loaded:", ", ".join(f"X0^{b.D}({b.N})[{len(b.known())} known]" for b in bases))

    # ---- global support survey: which N|m indices are NOT of theta type N*d*k^2 ----
    print("\n== support survey (N|m indices and their theta decomposition) ==")
    for b in bases:
        idxs = sorted({m for kf in b.known() for (blk, m) in b.coefs[kf] if blk == 'oo' and m % b.N == 0})
        desc = []
        for m in idxs:
            td = theta_decomp(m, b.D, b.N)
            desc.append(f"{m}->" + (f"(d={td[0]},k={td[1]})" if td else "NOT-THETA"))
        print(f"  X0^{b.D}({b.N}): {' '.join(desc) if desc else '(none)'}")

    # ---- THETA-const ----
    def u_theta_const(b, blk, m):
        if blk != 'oo': return None
        td = theta_decomp(m, b.D, b.N)
        if td is None: return None
        return (f"kap[{td[0]}]", F(1))
    rep, sol_tc, ok = solve_family(bases, None, u_theta_const, "THETA-const  w(N d k^2) = kappa_d")
    print("\n" + "\n".join(rep))

    # ---- THETA-klin ----
    def u_theta_klin(b, blk, m):
        if blk != 'oo': return None
        td = theta_decomp(m, b.D, b.N)
        if td is None: return None
        return (f"kap[{td[0]}]", F(td[1]))
    rep, _, _ = solve_family(bases, None, u_theta_klin, "THETA-klin   w(N d k^2) = kappa_d * k")
    print("\n" + "\n".join(rep))

    # ---- SIGMA controls: kappa * sum_{e|n} chi(e) e^a ----
    def kron(a, n):
        # Kronecker symbol (a/n), n>0
        if n == 0: return 0
        res = 1
        if n < 0: n = -n
        while n % 2 == 0:
            n //= 2
            r = a % 8
            if r in (3, 5): res = -res
        a %= n
        while a:
            while a % 2 == 0:
                a //= 2
                if n % 8 in (3, 5): res = -res
            a, n = n, a
            if a % 4 == 3 and n % 4 == 3: res = -res
            a %= n
        return res if n == 1 else 0
    for (a_exp, chi_m, label) in [(1, 1, "sigma1"), (0, 1, "sigma0"),
                                   (1, -4, "sigma1.chi-4"), (0, -4, "sigma0.chi-4"),
                                   (1, 8, "sigma1.chi8"), (1, -8, "sigma1.chi-8")]:
        def u_sigma(b, blk, m, a_exp=a_exp, chi_m=chi_m):
            if blk != 'oo' or m % b.N: return None
            n = m // b.N
            val = sum(kron(chi_m, e) * e**a_exp for e in divisors(n)) if chi_m != 1 \
                  else sum(e**a_exp for e in divisors(n))
            if val == 0: return None
            return ("kap", F(val))
        rep, _, _ = solve_family(bases, None, u_sigma, f"SIGMA {label}")
        print("\n" + "\n".join(rep))

    # ---- LAMBDA control: kappa * sum_{e|n} min(e, n/e) ----
    def u_lam(b, blk, m):
        if blk != 'oo' or m % b.N: return None
        n = m // b.N
        val = sum(min(e, n // e) for e in divisors(n))
        return ("kap", F(val))
    rep, _, _ = solve_family(bases, None, u_lam, "LAMBDA  w = kappa*sum min(e,n/e)")
    print("\n" + "\n".join(rep))

    # ---- kappa_d cross-base table for THETA-const ----
    print("\n== THETA-const kappa_d across bases (for pattern hunting) ==")
    for b in bases:
        sol = sol_tc.get(b.key())
        if sol is None:
            print(f"  X0^{b.D}({b.N}): (inconsistent)")
            continue
        items = ", ".join(f"{nm}={'FREE' if v is None else v}" for nm, v in sorted(sol.items()))
        print(f"  X0^{b.D}({b.N}) [D={b.D} N={b.N}]: {items}")

if __name__ == '__main__':
    main()
