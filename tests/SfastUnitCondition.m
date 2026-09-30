// Sfast (Popa's |S_N(u,t,n)|, the number of UNIT roots of alpha^2 - t alpha + n mod Nu) used
// to count all roots, via the discriminant and Hensel lifting.  For p | (n, N) the polynomial
// is alpha (alpha - t) mod p and the root alpha = 0 is not a unit, so every trace at
// gcd(n, N) > 1 was off: Popa's full-space formula disagreed with modular symbols on 80 of
// 208 small tuples, and the newform recursion on 72 (all inherited).  Neither paper is wrong:
// Popa's Theorem 4 holds "for all n >= 1", and [Assaf, Cor. 5.5] run on modular-symbol
// full-space traces matched on all 208.  The pipeline never hit it (every caller has p ∤ DN).
//
// (1) the count itself against brute force, where the defect is visible directly;
// (2) Popa's full-space trace and (3) the newform trace against modular symbols at
//     gcd(n, N) > 1, which the old code failed;
// (4) values at gcd(n, N) = 1 pinned from before the change, so the fix is inert there.
//
// Magma's HeckeOperator on a SUBSPACE dies at (N, k, n) = (6, 4, 2) ("incompatible
// coefficient rings"), so the oracle restricts the ambient operator, as tests/trace_formula.m.

_ := ClassNumberLU(-4);   // force ClassNumberData.m before the import (CLAUDE.md)
import "TraceFormula.m" : Sfast, S, TraceFormulaGamma0HeckeAL, TraceFormulaGamma0HeckeALNew;

function ModSymTrace(N, k, n, Q : New := false)
    M := ModularSymbols(N, k, 1);
    C := CuspidalSubspace(M);
    if New then C := NewSubspace(C); end if;
    if Dimension(C) eq 0 then return 0; end if;
    B := Matrix(Basis(VectorSpace(C)));
    return Trace(Solution(B, B*(HeckeOperator(M, n) * AtkinLehner(M, Q))));
end function;

// (1) Sfast vs brute force, with p | (n, N), u > 1, p = 2 and prime powers.  The old code
// fails 50 of the 154 (N, t, n) with u = 1, N in {12, 18, 20, 36, 45} -- e.g. (12, 0, 3): 2 vs 0.
cnt := 0;
for N in [4, 8, 9, 12, 16, 18, 20, 24, 27, 36, 45, 50, 72] do
    for u in Divisors(N) do
        for n in [2, 3, 4, 5, 6, 8, 9, 12, 18, 25] do
            for t in [-Floor(SquareRoot(4*n))..Floor(SquareRoot(4*n))] do
                if (t^2 - 4*n) mod u^2 ne 0 then continue; end if;
                assert Sfast(N, u, t, n) eq #S(N, u, t, n);
                cnt +:= 1;
            end for;
        end for;
    end for;
end for;
assert Sfast(12, 1, 0, 3) eq 0;      // the first case the old code got wrong (it said 2)
printf "  Sfast matches brute force on %o counts\n", cnt;

// (2), (3): full-space and newform traces at gcd(n, N) > 1.
full := 0; new := 0;
for N in [12, 18, 20, 36, 45, 72] do
    for Q in [Q : Q in Divisors(N) | GCD(Q, N div Q) eq 1] do
        for n in [2, 3, 4, 5, 9] do
            if GCD(n, N) eq 1 then continue; end if;
            for k in [2, 4] do
                assert TraceFormulaGamma0HeckeAL(N, k, n, Q) eq ModSymTrace(N, k, n, Q);
                full +:= 1;
                assert TraceFormulaGamma0HeckeALNew(N, k, n, Q) eq ModSymTrace(N, k, n, Q : New);
                new +:= 1;
            end for;
        end for;
    end for;
end for;
printf "  %o full-space and %o newform traces agree with modular symbols at gcd(n, N) > 1\n", full, new;

// (4) inert at gcd(n, N) = 1: TraceDNewALFixed values computed on main at 95cf87b.
W4 := {Integers() | 1, 8, 231, 1848};
assert TraceDNewALFixed(1, 1848, 2, 5, W4) eq -14;
assert TraceDNewALFixed(1, 1848, 2, 25, W4) eq 17;
assert TraceDNewALFixed(6, 307, 2, 25, {Integers() | 1, 6}) eq 36;
assert TraceDNewALFixed(58, 7, 2, 27, {Integers() | 1, 14, 58, 203}) eq -8;
assert TraceDNewALFixed(1, 4522, 2, 9, {Integers() | 1, 4522}) eq 323;
printf "  5 coprime-index pins unchanged\n";
