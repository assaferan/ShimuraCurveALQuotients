// Sfast is Popa's |S_N(u,t,n)|: the number of UNIT roots of alpha^2 - t alpha + n mod Nu.  It used
// to count all roots, through the discriminant and Hensel lifting.  For p | (n, N) the polynomial
// is alpha (alpha - t) mod p and the root alpha = 0 is not a unit, so every trace at gcd(n, N) > 1
// was off.  Neither paper is wrong: Popa's Theorem 4 holds for all n >= 1, and [Assaf, Cor. 5.5]
// run on modular-symbol full-space traces reproduces the newform traces.  The pipeline never hit
// it: every caller has p coprime to DN.
//
// (1) the count against brute force.  On this exact loop the OLD Sfast differs from brute
//     force on 774 of 2677 counts (measured on main at 95cf87b, scratch script
//     vvdata/weyl-campaign/trace-formula-gcd-2026-09-30/old_sfast_count.m on the campaign branch is this loop).
// (2) Popa's full-space trace and (3) the newform trace against Magma's modular symbols at
//     gcd(n, N) > 1 on the grid below.  On this grid main at 95cf87b got 70 of the 176
//     full-space traces wrong and 64 of the 176 newform traces (all inherited from Sfast at
//     sub-levels); measured by vvdata/weyl-campaign/trace-formula-gcd-2026-09-30/grid176.m.
// (4) the D-new, W-averaged trace TraceDNewALFixed at coprime index against modular symbols:
//     the trace of T_n o w on the cuspidal subspace of level DN that is new at every p | D,
//     averaged over w in W, with the Jacquet-Langlands sign (-1)^omega(gcd(w, D)) -- so the
//     change is inert where the pipeline lives.
// (5) TraceDNewALFixed at n sharing a factor with DN, against the same modular-symbols
//     quantity.  The D-new trace used to keep only the n' = 1 term of [Assaf, Cor. 4.27],
//     which is the whole formula exactly when gcd(n, DN) = 1; at gcd > 1 it returned wrong
//     values.  D = 1, N = 30, n = 4: main at 95cf87b gives 3, main with only the Sfast fix of
//     this PR gives -1, modular symbols give -2 (LMFDB 30.2.a.a together with the old forms
//     from 15.2.a.a); (1, 1848, 7) with the V4 below: main gives -20, modular symbols -9.  It
//     now sums Lemma 4.20 block by block over the N' divisible by D, for every n.  Beyond the
//     cases here, 0 of 1548 such traces differed from modular symbols across 23 (D, N) pairs
//     (D up to 26), several W each, n in {2..25}
//     (vvdata/weyl-campaign/trace-formula-gcd-2026-09-30/dnew_general.m on the campaign branch).
//
// Magma's HeckeOperator on a SUBSPACE dies at (N, k, n) = (6, 4, 2) ("incompatible
// coefficient rings"), so the oracles restrict the ambient operator, as tests/trace_formula.m.

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

// D-new, W-averaged trace from modular symbols (weight 2), the quantity TraceDNewALFixed returns.
function DNewTraceModSym(D, N, n, W)
    M := ModularSymbols(D*N, 2, 1);
    C := CuspidalSubspace(M);
    for p in PrimeDivisors(D) do C := NewSubspace(C, p); end for;
    if Dimension(C) eq 0 then return 0; end if;
    B := Matrix(Basis(VectorSpace(C)));
    T := HeckeOperator(M, n);
    tot := 0;
    for w in W do
        tot +:= (-1)^#PrimeDivisors(GCD(w, D)) * Trace(Solution(B, B*(T * AtkinLehner(M, w))));
    end for;
    return tot / #W;
end function;

// (1)
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
assert cnt eq 2677;
assert Sfast(12, 1, 0, 3) eq 0;      // both roots of alpha^2 + 3 mod 3 are 0: no unit root (old code: 2)
printf "  Sfast matches brute force on %o counts\n", cnt;

// (2), (3)
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
assert full eq 176 and new eq 176;
printf "  %o full-space and %o newform traces agree with modular symbols at gcd(n, N) > 1\n", full, new;

// (4): each value below is the modular-symbols quantity described in the header, computed
// here, not a number copied from an earlier run of this code.
W4 := {Integers() | 1, 8, 231, 1848};
for c in [<1, 1848, 5, W4>, <1, 1848, 25, W4>, <6, 307, 25, {Integers() | 1, 6}>,
          <58, 7, 27, {Integers() | 1, 14, 58, 203}>, <1, 4522, 9, {Integers() | 1, 4522}>] do
    D, N, n, W := Explode(c);
    assert GCD(n, D*N) eq 1;
    assert TraceDNewALFixed(D, N, 2, n, W) eq DNewTraceModSym(D, N, n, W);
end for;
printf "  5 coprime-index D-new traces agree with modular symbols\n";

// (5)
for c in [<1, 30, 4, {Integers() | 1}>, <1, 30, 6, {Integers() | 1, 30}>, <1, 1848, 7, W4>,
          <1, 1848, 49, W4>, <6, 7, 3, {Integers() | 1, 6}>, <6, 7, 4, {Integers() | 1, 2, 3, 6}>,
          <10, 9, 6, {Integers() | 1, 10}>, <15, 4, 10, {Integers() | 1, 3, 5, 15}>] do
    D, N, n, W := Explode(c);
    assert GCD(n, D*N) gt 1;
    assert TraceDNewALFixed(D, N, 2, n, W) eq DNewTraceModSym(D, N, n, W);
end for;
assert TraceDNewALFixed(1, 30, 2, 4, {Integers() | 1}) eq -2;      // LMFDB 30.2.a.a + old forms of 15.2.a.a
printf "  8 D-new traces at n sharing a factor with DN agree with modular symbols\n";
