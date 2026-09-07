import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

// tests/X0_22_3.m -- RE-DERIVATION test for X_0^22(3).
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

function load_covers_and_ws_data_22_3()
    _<s> := PolynomialRing(Rationals());


    cover_data := AssociativeArray();
    cover_data[{1,3,22,66}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, -729/4, 729/4 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,2,11,22}] := <HyperellipticCurve(Polynomial(Rationals(), [ -11/186624, 35/373248, -131/2985984, 1/110592 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1}] := <HyperellipticCurve(Polynomial(Rationals(), [ -3/4096, 0, -77/9216, 0, -1073/18432, 0, -77/9216, 0, -3/4096 ])), DiagonalMatrix([1,1,1])>;   // genus 3
    cover_data[{1,6,11,66}] := <HyperellipticCurve(Polynomial(Rationals(), [ 9, -9 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,3}] := <HyperellipticCurve(Polynomial(Rationals(), [ -11/2304, 0, 31/4608, 0, -11/4096 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,33}] := <HyperellipticCurve(Polynomial(Rationals(), [ -11/2304, 0, -13/10368, 0, -1/6912 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,3,11,33}] := <HyperellipticCurve(Polynomial(Rationals(), [ -11/2304, 13/4608, -3/4096 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,2}] := <HyperellipticCurve(Polynomial(Rationals(), [ -11/186624, 0, -35/839808, 0, -131/15116544, 0, -1/1259712 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,2,33,66}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, -9/4 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,6}] := <HyperellipticCurve(Polynomial(Rationals(), [ 11/16384, 0, -49/1327104, 0, -23/11943936, 0, -1/3981312 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,6,22,33}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, 11/9216, -13/18432, 3/16384 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,2,3,6}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, 11/746496, -35/1492992, 131/11943936, -1/442368 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,22}] := <HyperellipticCurve(Polynomial(Rationals(), [ -11/9216, 0, 53/18432, 0, -347/147456, 0, 11/16384 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,11}] := <HyperellipticCurve(Polynomial(Rationals(), [ -11/4096, 0, -25/165888, 0, -1/110592 ])), DiagonalMatrix([1,1,1])>;   // genus 1

    ws_data := AssociativeArray();
    return cover_data, ws_data;
end function;

procedure test_22_3()
    cover_data, ws_data := load_covers_and_ws_data_22_3();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(22, 3, cover_data, ws_data, curves);
    return;
end procedure;

test_22_3();
