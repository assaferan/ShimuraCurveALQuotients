#!/usr/bin/env python3
"""The m=0 local density at a p-scaled plane, by brute-force counting.

Two purposes.

(1) CORRECTION to km0coset2.m's labelling.  It calls  p*diag(1,u), u a
    nonsquare, the "p inert" (anisotropic) case.  That form is
    Q = p(x^2 + u y^2) = p(x^2 - (-u) y^2), which is anisotropic iff -u is a
    NONsquare, i.e. iff p = 1 mod 4.  For p = 3 mod 4 it is isotropic -- a
    hyperbolic plane in disguise -- which is exactly why p = 3, 7 showed the
    same 2(p-1) count as the split case while p = 5 showed 0.  The genuinely
    anisotropic plane is p(x^2 - u y^2), and it has NO nonzero isotropic
    cosets, at every p.  So the paper's reading of the count 2(p-1) as the
    nonzero isotropic cosets of a hyperbolic plane is correct after all.

(2) Whether the anisotropic factor can supply the zero that would make the
    m=0 correction an s-DERIVATIVE rather than a VALUE.  It cannot vanish:
    the density is bounded below by 1.

alpha_k = p^{-k(n-1)} #{x mod p^k : Q(x) = 0 mod p^k}, n = 2, so divide by p^k.
Neither case converges -- which is why the defining series does not terminate
and why the s-variable is unavoidable -- but they diverge differently:

    anisotropic   alpha_k = p, 1, p, 1, ...        (bounded, period 2)
    hyperbolic    alpha_k = (p-1)k + 1             (linear growth)

As formal series in X = p^{-s}, summing over k >= 0:

    anisotropic   sum alpha_k X^k = (1 + pX) / ((1-X)(1+X))     simple pole
    hyperbolic    sum alpha_k X^k = (1 + (p-2)X) / (1-X)^2      double pole

so after the natural (1-X) normalisation the anisotropic place gives
(1+p)/2 at X=1 -- finite and NONZERO -- while the split place at mu=0 still
has a pole (it carries the L-factor).  A bounded, nonvanishing density cannot
produce a zero under any normalisation that is a nonzero constant times a
power of (1-X).  The exact normalising constant is KRY part (i) at m=0 and is
NOT settled here.
"""
from fractions import Fraction as F


def squares(p):
    return {(x * x) % p for x in range(p)}


def nonsquare(p):
    S = squares(p)
    return next(x for x in range(1, p) if x not in S)


def isotropic_cosets(p, form):
    """#{mu != 0 in (Z/p)^2 : Q(mu) = 0 in Q/Z} for the p-scaled plane."""
    return sum(1 for a in range(p) for b in range(p)
               if (a, b) != (0, 0) and form(a, b) % p == 0)


def density(p, k, form):
    pk = p ** k
    N = sum(1 for x in range(pk) for y in range(pk) if (p * form(x, y)) % pk == 0)
    return N, F(N, pk)


def main():
    print("(1) Is  p*diag(1,u)  -- km0coset2.m's \"p inert\" -- actually anisotropic?\n")
    print("   p   u   -u   -u square?   anisotropic?   #nonzero isotropic cosets")
    for p in [3, 5, 7, 11, 13]:
        S, u = squares(p), nonsquare(p)
        mu = (-u) % p
        n = isotropic_cosets(p, lambda a, b, u=u: a * a + u * b * b)
        print(f"  {p:3}  {u:2}  {mu:3}   {str(mu in S):>5}        {str(mu not in S):>5}"
              f"          {n:3}   (2(p-1) = {2*(p-1)})")
    print("\n    => anisotropic only for p = 1 mod 4; the p = 3 mod 4 cases are"
          "\n       hyperbolic planes mislabelled, which explains the 2(p-1) count there.\n")

    print("   the genuinely anisotropic plane  p(x^2 - u y^2):")
    for p in [3, 5, 7, 11, 13]:
        u = nonsquare(p)
        n = isotropic_cosets(p, lambda a, b, u=u: a * a - u * b * b)
        print(f"  {p:3}  nonzero isotropic cosets: {n}")
    print()

    print("(2) m=0 density, ANISOTROPIC  Q = p(x^2 - u y^2), mu = 0")
    for p in [3, 5, 7]:
        u = nonsquare(p)
        row = []
        for k in range(1, 7):
            if p ** k > 3000:
                break
            row.append(str(density(p, k, lambda x, y, u=u: x * x - u * y * y)[1]))
        print(f"  p={p}:  alpha_k = {', '.join(row)}   (bounded, period 2)")
    print()
    print("    m=0 density, HYPERBOLIC  Q = p*x*y, mu = 0")
    for p in [3, 5]:
        row = []
        for k in range(1, 6):
            if p ** k > 3000:
                break
            row.append(str(density(p, k, lambda x, y: x * y)[1]))
        print(f"  p={p}:  alpha_k = {', '.join(row)}   ((p-1)k + 1, linear growth)")
    print("\n    => the anisotropic density is bounded below by 1 and cannot vanish."
          "\n       So no place vanishes at m=0 for a nonzero isotropic mu, and the"
          "\n       m=0 coefficient there is a VALUE, not an s-derivative -- MODULO"
          "\n       the normalising constant, which is KRY part (i) at m=0 and is open.")


if __name__ == "__main__":
    main()
