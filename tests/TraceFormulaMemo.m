// Issue #8: TraceDNew used to call the level-N' newform trace once per element of get_ds
// (the level-raising indices d), and TraceFormulaGamma0HeckeAL (Popa's formula) was
// re-evaluated at every sub-level once per level above it and once per w in W.  The fix
// multiplies by #get_ds and memoises Popa on <N, k, n, Q>.  Both are value-preserving by
// construction, so the checks are against an outside source, Magma's modular symbols:
//
// (1) TraceDNewALFixed against the D-new, W-averaged trace: the trace of T_n o w on the
//     cuspidal subspace of level DN that is new at every p | D, averaged over w in W, with the
//     Jacquet-Langlands sign (-1)^omega(gcd(w, D)).  Coprime n only: at gcd(n, DN) > 1 the
//     D-new decomposition in TraceDNew is wrong (D = 1, N = 30, n = 4 gives -1, modular
//     symbols -2), a pre-existing defect that PR #57 makes TraceDNewALFixed refuse.
// (2) the newform trace against modular symbols at non-squarefree levels, n coprime to N,
//     which exercises the alpha = 0 skip (alpha vanishes for p | Q at exponent 1 and off
//     cube-free N/N').  At gcd(n, N) > 1 the trace formula was wrong before #57 (Sfast counted
//     non-unit roots), so that case is tested there, not here.
// (3) TraceDNewQuotient(V2, 1848) = genus(X_0(1848)/V2) = (g + Tr V2)/2 with g = 369 and
//     Tr(V2 | S_2(1848)) = 1, both from modular symbols here.
// (4) the memo layer exists and the recursion has its measured shape: one call at level 1848,
//     #W = 4, on a cleared cache leaves 54 distinct Popa arguments (0 keys = memo layer gone).
//
// None of this guards the SPEED: reverting either change leaves every value equal.

_ := ClassNumberLU(-4);   // force ClassNumberData.m before the import (CLAUDE.md)
import "TraceFormula.m" : TraceFormulaGamma0HeckeALNew;
import "Caching.m" : cached_popa;
import !"Geometry/ModSym/operators.m" : ActionOnModularSymbolsBasis;

// Magma's HeckeOperator on a SUBSPACE dies at (N, k, n) = (6, 4, 2), so restrict the ambient
// operator, as tests/trace_formula.m does.
function ModSymTrace(N, k, n, Q : New := false)
    M := ModularSymbols(N, k, 1);
    C := CuspidalSubspace(M);
    if New then C := NewSubspace(C); end if;
    if Dimension(C) eq 0 then return 0; end if;
    B := Matrix(Basis(VectorSpace(C)));
    return Trace(Solution(B, B*(HeckeOperator(M, n) * AtkinLehner(M, Q))));
end function;

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

// (4) first, on a cleared cache.
W4 := {Integers() | 1, 8, 231, 1848};
CacheClear(cached_popa);
v := TraceDNewALFixed(1, 1848, 2, 5, W4);
b, cache := StoreIsDefined(cached_popa, "cache");
assert b;
printf "  Popa cache after one call at (1, 1848), #W = 4: %o keys\n", #Keys(cache);
assert #Keys(cache) eq 54;

// (1)
W16  := {Integers() | 1, 6, 7, 10, 11, 15, 42, 66, 70, 77, 105, 110, 165, 462, 770, 1155};
cases := [
    <1, 1848, 5, W4>, <1, 1848, 25, W4>, <1, 1848, 125, W4>,
    <6, 307, 5, {Integers() | 1, 6}>, <6, 307, 25, {Integers() | 1, 6}>, <6, 307, 11, {Integers() | 1, 6}>,
    <10, 231, 13, W16>, <10, 231, 169, W16>,
    <14, 15, 11, {Integers() | 1, 5, 42, 210}>, <14, 15, 121, {Integers() | 1, 5, 42, 210}>,
    <14, 15, 1331, {Integers() | 1, 5, 42, 210}>,
    <1, 4522, 3, {Integers() | 1, 4522}>, <1, 4522, 9, {Integers() | 1, 4522}>,
    <58, 7, 3, {Integers() | 1, 14, 58, 203}>, <58, 7, 9, {Integers() | 1, 14, 58, 203}>,
    <58, 7, 27, {Integers() | 1, 14, 58, 203}> ];
for c in cases do
    D, N, n, W := Explode(c);
    assert GCD(n, D*N) eq 1;
    assert TraceDNewALFixed(D, N, 2, n, W) eq DNewTraceModSym(D, N, n, W);
end for;
printf "  %o D-new traces agree with modular symbols\n", #cases;

// (3)
M := ModularSymbols(1848, 2, 1);
C := CuspidalSubspace(M);
g := Dimension(C);
B := Matrix(Basis(VectorSpace(C)));
trV2 := Trace(Solution(B, B*ActionOnModularSymbolsBasis(Eltseq(get_V2(1848)), M)));
assert g eq 369 and trV2 eq 1;
gq := Integers()!(g + trV2) div 2;          // trV2 comes back as a rational from Solution
assert TraceDNewQuotient(get_V2(1848), "V2", 1, {Integers() | 1}, 1, 1848) eq gq;
printf "  genus(X_0(1848)/V2) = (%o + %o)/2 = %o\n", g, trV2, gq;

// (2)
checks := 0;
for N in [12, 18, 20, 36, 45, 72, 108] do
    Qs := [Q : Q in Divisors(N) | GCD(Q, N div Q) eq 1];
    for Q in Qs do
        for n in [n : n in [1, 2, 3, 4, 5, 7, 9, 25] | GCD(n, N) eq 1] do
            for k in [2, 4] do
                assert TraceFormulaGamma0HeckeALNew(N, k, n, Q) eq ModSymTrace(N, k, n, Q : New);
                checks +:= 1;
            end for;
        end for;
    end for;
end for;
printf "  %o newform traces agree with modular symbols at non-squarefree levels\n", checks;
