// Issue #8: TraceDNew used to call the level-N' newform trace once per element of get_ds
// (the level-raising indices d), and TraceFormulaGamma0HeckeAL (Popa's formula) was
// re-evaluated at every sub-level once per level above it and once per w in W.  The fix
// multiplies by #get_ds and memoises Popa on <N, k, n, Q>.  Both are value-preserving by
// construction, so the test is (1) the values TraceDNewALFixed returned BEFORE the change,
// computed on main at 95cf87b, and (2) the newform trace against modular symbols, including
// non-squarefree levels where alpha vanishes and the new `continue` fires.
//
// The pins do not guard the SPEED: reverting either change leaves every value below equal.
// The key count in (3) does pin that the memo layer exists and the recursion has the shape
// it had when measured (54 distinct Popa arguments for one call at level 1848, #W = 4).

_ := ClassNumberLU(-4);   // force ClassNumberData.m before the import (CLAUDE.md)
import "TraceFormula.m" : TraceFormulaGamma0HeckeALNew, get_trace_hecke_AL;
import "Caching.m" : cached_popa;

// (1) values from the unchanged code.  <D, N, n, W, TraceDNewALFixed(D, N, 2, n, W)>
W4   := {Integers() | 1, 8, 231, 1848};
W16  := {Integers() | 1, 6, 7, 10, 11, 15, 42, 66, 70, 77, 105, 110, 165, 462, 770, 1155};
pins := [
    <1, 1848, 5, W4, -14>, <1, 1848, 25, W4, 17>, <1, 1848, 125, W4, 12>,
    <1, 1848, 7, W4, -20>, <1, 1848, 49, W4, 96>,
    <6, 307, 5, {Integers() | 1, 6}, 0>, <6, 307, 25, {Integers() | 1, 6}, 36>,
    <6, 307, 11, {Integers() | 1, 6}, 0>,
    <10, 231, 2, W16, 0>, <10, 231, 4, W16, 41>, <10, 231, 8, W16, 7>,
    <10, 231, 13, W16, 0>, <10, 231, 169, W16, -31>,
    <14, 15, 11, {Integers() | 1, 5, 42, 210}, 0>, <14, 15, 121, {Integers() | 1, 5, 42, 210}, -1>,
    <14, 15, 1331, {Integers() | 1, 5, 42, 210}, 0>,
    <1, 4522, 3, {Integers() | 1, 4522}, -24>, <1, 4522, 9, {Integers() | 1, 4522}, 323>,
    <58, 7, 3, {Integers() | 1, 14, 58, 203}, -2>, <58, 7, 9, {Integers() | 1, 14, 58, 203}, 5>,
    <58, 7, 27, {Integers() | 1, 14, 58, 203}, -8> ];

// (3) first, on a cleared cache: one call at level 1848 with #W = 4 needs 54 distinct
// Popa evaluations (it made ~1260 before the change).  0 keys means the memo layer is gone.
CacheClear(cached_popa);
assert TraceDNewALFixed(1, 1848, 2, 5, W4) eq -14;
b, cache := StoreIsDefined(cached_popa, "cache");
assert b;
printf "  Popa cache after one call at (1, 1848), #W = 4: %o keys\n", #Keys(cache);
assert #Keys(cache) eq 54;

for pin in pins do
    D, N, n, W, v := Explode(pin);
    got := TraceDNewALFixed(D, N, 2, n, W);
    assert got eq v;
end for;
printf "  %o TraceDNewALFixed pins reproduce\n", #pins;

// The non-AL path goes through the same TraceDNewALFixed (with n = 1) plus TraceFormulaGamma0VW.
assert TraceDNewQuotient(get_V2(1848), "V2", 1, {Integers() | 1}, 1, 1848) eq 185;

// (2) newform trace against modular symbols.  Non-squarefree N exercises alpha = 0 (e.g.
// N/N' = 8, 27) where TraceFormulaGamma0HeckeALNew now skips the term.
// ⚠ Only n coprime to N: at gcd(n, N) > 1 the recursion disagrees with modular symbols on
// 72 of 336 tuples of this grid, identically before and after this change (measured on main
// at 95cf87b) -- a pre-existing boundary of the implementation, which assumes (n, DN) = 1
// (TraceDNew's comment says so, and every pipeline caller has p coprime to DN).
checks := 0;
for N in [12, 18, 20, 36, 45, 72, 108] do
    Qs := [Q : Q in Divisors(N) | GCD(Q, N div Q) eq 1];
    for Q in Qs do
        for n in [n : n in [1, 2, 3, 4, 5, 7, 9, 25] | GCD(n, N) eq 1] do
            for k in [2, 4] do
                assert TraceFormulaGamma0HeckeALNew(N, k, n, Q) eq get_trace_hecke_AL(N, k, n, Q : New);
                checks +:= 1;
            end for;
        end for;
    end for;
end for;
printf "  %o newform traces agree with modular symbols at non-squarefree levels\n", checks;
