#!/usr/bin/env python3
"""Prop 9.15 (the multiplicity coefficients) in closed form, with a Hurwitz class
number computed from reduced forms rather than from a table.

    -a_E(m) = 12 H(4m) * prod_{p|D} ram_p * prod_{p|N} eich_p,     -4m = D_0 c^2

Self-test (`python3 closedcoef.py`) checks H against known values and then checks
-a_E against the banked genus-combination coefficients of ctwm.log at D=15,39,55.
The point of the second check is that ctwm.log is produced by a different route
(Magma genus sums), so this is not a regression test against the same code.
"""
from fractions import Fraction as F
import math


def legendre(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def kron(a, p):
    """Kronecker symbol (a/p), p prime, p=2 included."""
    if p == 2:
        if a % 2 == 0:
            return 0
        return 1 if a % 8 in (1, 7) else -1
    return legendre(a, p)


def factor(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def isfund(D):
    """Is D a negative fundamental discriminant?"""
    if D >= 0:
        return False
    if D % 4 == 1:
        return all(e == 1 for p, e in factor(-D).items())
    if D % 4 == 0:
        q = D // 4
        if q % 4 in (2, 3):
            return all(e == 1 for p, e in factor(-q).items())
    return False


def decomp(m):
    """-4m = D_0 c^2 with D_0 fundamental."""
    n = -4 * m
    for c in range(int(math.isqrt(4 * m)), 0, -1):
        if n % (c * c) == 0 and isfund(n // (c * c)):
            return n // (c * c), c
    raise ValueError(f"no fundamental decomposition for m={m}")


def H(N):
    """Hurwitz class number: weighted count of reduced positive definite binary
    forms of discriminant -N.  Weight 1/2 on (a,0,a), 1/3 on (a,a,a)."""
    if N == 0:
        return F(-1, 12)
    if N % 4 not in (0, 3):
        return F(0)
    tot = F(0)
    a = 1
    while a * a <= N:
        for b in range(-a + 1, a + 1):
            num = b * b + N
            if num % (4 * a):
                continue
            c = num // (4 * a)
            if c < a:
                continue
            if a == c and b < 0:
                continue
            if b == 0 and a == c:
                tot += F(1, 2)
            elif a == b == c:
                tot += F(1, 3)
            else:
                tot += 1
        a += 1
    return tot


def W(D, N, m):
    """-a_E(m) for the base X_0^D(N), per Prop 9.15."""
    D0, c = decomp(m)
    tot = F(12) * H(4 * m)
    for p in factor(D):
        k = 0
        cc = c
        while cc % p == 0:
            cc //= p
            k += 1
        chi = kron(D0, p)
        sk = sum(p ** i for i in range(k + 1))
        skm = sum(p ** i for i in range(k)) if k else 0
        tot *= F(1 - chi, (p - 1) * (sk - chi * skm))
    for p in (factor(N) if N > 1 else {}):
        k = 0
        cc = c
        while cc % p == 0:
            cc //= p
            k += 1
        chi = kron(D0, p)
        sk = sum(p ** i for i in range(k + 1))
        skm = sum(p ** i for i in range(k)) if k else 0
        tot *= F((p - chi) * p ** k, (p * p - 1) * (sk - chi * skm))
    return tot


if __name__ == "__main__":
    bad = 0
    print("Hurwitz class numbers:")
    for n, exp in [(3, F(1, 3)), (4, F(1, 2)), (7, 1), (8, 1), (11, 1),
                   (12, F(4, 3)), (16, F(3, 2)), (28, 2), (40, 2), (48, F(10, 3))]:
        got = H(n)
        bad += got != exp
        print(f"  H({n:2}) = {str(got):>5}  expected {str(exp):>5}"
              f"  {'ok' if got == exp else '** MISMATCH **'}")

    # ctwm.log: coefficients a(0..12) of the genus-theta combination at N=1.
    ctwm = {15: [None, 0, 0, -4, 0, 0, 0, -12, 0, 0, -6, 0, -10],
            39: [None, 0, 0, 0, 0, 0, -2, -4, 0, 0, 0, 0, 0],
            55: [None, 0, 0, F(-8, 5), 0, F(-6, 5), 0, 0, 0, 0, 0, 0, -4]}
    print("\n-a_E(m) vs the banked genus combination (ctwm.log), N=1:")
    for D, row in sorted(ctwm.items()):
        res = [(m, W(D, 1, m), -row[m]) for m in range(1, 13)]
        hit = sum(1 for _, g, e in res if g == e)
        bad += len(res) - hit
        print(f"  D={D:3}: {hit}/{len(res)}")
        for m, g, e in res:
            if g != e:
                print(f"     m={m}: closed form {g}, ctwm.log {e}")
    print(f"\n{'ALL CHECKS PASS' if not bad else str(bad) + ' MISMATCHES'}")
