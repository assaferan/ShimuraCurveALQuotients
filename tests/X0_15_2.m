import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

// tests/X0_15_2.m -- RE-DERIVATION test for X_0^15(2).
//
// Re-runs AllEquationsAboveCovers and compares what the pipeline actually produces, so passing IS
// reproduction -- the stronger claim than GuoYangEquations.m's stored-model comparison.
//
// ⚠ THIS BASE INCLUDES ITS `CRV` (paired) ENTRY, which was impossible before 2026-09-07. The
// helper now CONSTRUCTS the isomorphism for CRV pairs (tests/_crviso.m) instead of calling
// IsIsomorphic, which HANGS on them -- the genus-3 pair at 14_3 ran >50 min, the genus-5 one at
// 26_3 >1 h. The ambient weights are DERIVED (y's weight is half its own equation's degree), which
// tests/CRVStructure.m verifies against every stored entry.
//
// Keys with several stored entries are omitted: they cannot be matched unambiguously.
// The matrix in each cover_data value is a placeholder -- with manual_isomorphism false (the
// default) the helper never reads it.
//
// ⚠ THIS TEST IS WEAKER THAN THE HAND-WRITTEN ONES: `ws_data` IS EMPTY, so it makes ZERO
// involution comparisons. It verifies that each cover is ISOMORPHIC to the stored curve, but not
// that the Atkin-Lehner involutions correspond -- and the involutions are what make these QUOTIENT
// models rather than merely curves. 23 of the 34 X0_*.m tests do check them; the ones generated on
// 2026-09-07 (this file among them) do not.
// ⇒ Closing that needs involution matrices IN OUR MODEL'S COORDINATES. Taking them from the
// pipeline's own `ws` output would be circular; deriving them from Guo-Yang's published
// involutions is independent but is per-base work. Recorded rather than silently accepted.

function load_covers_and_ws_data_15_2()
    _<s> := PolynomialRing(Rationals());


    cover_data := AssociativeArray();
    cover_data[{1,10}] := <HyperellipticCurve(Polynomial(Rationals(), [ -25/48, 0, -17/8, 0, -27/16 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,3,10,30}] := <HyperellipticCurve(Polynomial(Rationals(), [ 5/144, 7/9, -4/3 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1}] := <HyperellipticCurve(Polynomial(Rationals(), [ -7/9, -13/9, -85/36, -103/36, -371/144, -19/12, -47/72, -1/6, -1/48 ])), DiagonalMatrix([1,1,1])>;   // genus 3
    cover_data[{1,2}] := <HyperellipticCurve(Polynomial(Rationals(), [ -5/2, 0, -107/16, 0, 19/8, 0, -3/16 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,2,5,10}] := <HyperellipticCurve(Polynomial(Rationals(), [ -5/2, -107/2, 152, -96 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,5}] := <HyperellipticCurve(Polynomial(Rationals(), [ -5/2, 0, 127/2, 0, -47/2, 0, -75/2 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,6,10,15}] := <HyperellipticCurve(Polynomial(Rationals(), [ -8/9, 8/9 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,2,15,30}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, 8 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,6}] := <HyperellipticCurve(Polynomial(Rationals(), [ -25/6, 0, -347/16, 0, -261/8, 0, -243/16 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,2,3,6}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, -20, -428, 1216, -768 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,30}] := <HyperellipticCurve(Polynomial(Rationals(), [ 5/144, 0, 7/72, 0, -1/48 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,3,5,15}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, -64/9, 64/9 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,5,6,30}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, 5/18, 56/9, -32/3 ])), DiagonalMatrix([1,1,1])>;   // genus 1

    ws_data := AssociativeArray();
    return cover_data, ws_data;
end function;

procedure test_15_2()
    cover_data, ws_data := load_covers_and_ws_data_15_2();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(15, 2, cover_data, ws_data, curves);
    return;
end procedure;

test_15_2();
