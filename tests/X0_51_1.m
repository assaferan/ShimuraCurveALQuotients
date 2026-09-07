import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

// tests/X0_51_1.m -- RE-DERIVATION test for X_0^51(1).
//
// Unlike tests/GuoYangEquations.m, which reads the STORED model, this re-runs
// AllEquationsAboveCovers and compares what the pipeline actually produces. Passing therefore
// IS reproduction, which is the stronger claim (see data/models/PROVENANCE.md).
//
// The expected curves below are the committed model's hyperelliptic cover entries, which are
// themselves validated against Guo-Yang's published equations. CRV (paired) entries and empty
// keys are omitted: the helper compares hyperelliptic covers, and a key with several stored
// entries cannot be matched unambiguously.
//
// The second component of each cover_data value is unused here -- with manual_isomorphism false
// (the default) the helper calls IsIsomorphic, so the matrix is a placeholder. ws_data is left
// empty for the same reason: the helper skips involution checks for keys it does not find.
//
// ✅ INVOLUTIONS NOW CHECKED (2026-09-07): `ws_data` carries Guo-Yang's published w_3 and w_51,
// so this test makes 2 involution comparisons as well as 4 curve comparisons. NEGATIVE-CONTROLLED:
// swapping the two matrices makes it FAIL ("no identification matches the involution LABELLING"),
// so the check discriminates rather than merely agreeing.
// The note below applies to the OTHER generated tests, which still have empty ws_data.
// ⚠ THE OTHER GENERATED TESTS ARE WEAKER: their `ws_data` IS EMPTY, so they make ZERO
// involution comparisons. It verifies that each cover is ISOMORPHIC to the stored curve, but not
// that the Atkin-Lehner involutions correspond -- and the involutions are what make these QUOTIENT
// models rather than merely curves. 23 of the 34 X0_*.m tests do check them; the ones generated on
// 2026-09-07 (this file among them) do not.
// ⇒ Closing that needs involution matrices IN OUR MODEL'S COORDINATES. Taking them from the
// pipeline's own `ws` output would be circular; deriving them from Guo-Yang's published
// involutions is independent but is per-base work. Recorded rather than silently accepted.

function load_covers_and_ws_data_51_1()
    _<s> := PolynomialRing(Rationals());

    cover_data := AssociativeArray();
    cover_data[{1,3}] := <HyperellipticCurve(Polynomial(Rationals(), [ -2187/256, -5589/64, -27297/128, 6507/64, -2187/256 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,17}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, 6561/256, 16767/64, 81891/128, -19521/64, 6561/256 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,51}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, -3 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1}] := <HyperellipticCurve(Polynomial(Rationals(), [ -2187/256, 0, 1863/64, 0, -3033/128, 0, -241/64, 0, -27/256 ])), DiagonalMatrix([1,1,1])>;   // genus 3

    // ATKIN-LEHNER INVOLUTIONS, transcribed from Guo-Yang (Compositio 153 (2017), level-one table):
    //     w_3 (x,y) = (-x,  y)        w_51 (x,y) = ( x, -y)
    // ⚠ WHY THESE CARRY OVER UNCHANGED to our coordinates. Our model differs from Guo-Yang's by
    // x -> 3x, y -> (27/16)y (recorded in tests/GuoYangEquations.m). That change is DIAGONAL, and
    // diagonal scalings COMMUTE with sign changes, so an involution of the form (+-x, +-y) has the
    // same matrix in both coordinate systems. This would NOT hold for Guo-Yang's non-linear
    // involutions (e.g. 55_1's w_5(x,y) = (-1/x, y/x^4)), which the helper cannot express as a
    // matrix at all.
    // The ambient of a HyperellipticCurve has coordinates (x, y, z); z is unaffected.
    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][3]  := DiagonalMatrix([-1, 1, 1]);
    ws_data[{1}][51] := DiagonalMatrix([ 1,-1, 1]);
    return cover_data, ws_data;
end function;

procedure test_51_1()
    cover_data, ws_data := load_covers_and_ws_data_51_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(51, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_51_1();
