import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

// tests/X0_26_3.m -- RE-DERIVATION test for X_0^26(3).
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
//
// ⚠ WHAT A PASS DOES AND DOES NOT PIN, measured 2026-09-07 by negative control on X0_14_3.
// The comparison is `IsIsomorphic`, so it pins each cover UP TO ISOMORPHISM -- which is the right
// notion for a model, but it is WEAKER on the genus-0 entries than it looks:
//   * perturbing the GENUS-2 entry at {1,3} (-184 -> -185) makes the test FAIL, as it should;
//   * perturbing a GENUS-0 conic's constant term does NOT -- two conics can be isomorphic with
//     different coefficients, so `IsIsomorphic` correctly still says yes.
// So the genus-0 cover entries are checked only for their isomorphism CLASS. If a conic's actual
// coefficients matter (they do for the CRV constructions, and for `15_4`'s twist), that has to be
// pinned elsewhere -- see tests/CRVFullCurve.m and tests/CRV_15_4.m.
// The second component of each cover_data value is unused here -- with manual_isomorphism false
// (the default) the helper calls IsIsomorphic, so the matrix is a placeholder. ws_data is left
// empty for the same reason: the helper skips involution checks for keys it does not find.
//
// ⚠ THIS TEST IS WEAKER THAN THE HAND-WRITTEN ONES: `ws_data` IS EMPTY, so it makes ZERO
// involution comparisons. It verifies that each cover is ISOMORPHIC to the stored curve, but not
// that the Atkin-Lehner involutions correspond -- and the involutions are what make these QUOTIENT
// models rather than merely curves. 23 of the 34 X0_*.m tests do check them; the ones generated on
// 2026-09-07 (this file among them) do not.
// ⇒ Closing that needs involution matrices IN OUR MODEL'S COORDINATES. Taking them from the
// pipeline's own `ws` output would be circular; deriving them from Guo-Yang's published
// involutions is independent but is per-base work. Recorded rather than silently accepted.

function load_covers_and_ws_data_26_3()
    _<s> := PolynomialRing(Rationals());

    cover_data := AssociativeArray();
    cover_data[{1,6,26,39}] := <HyperellipticCurve(Polynomial(Rationals(), [ -11, 16 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,3}] := <HyperellipticCurve(Polynomial(Rationals(), [ -3/8, 0, -91/64, 0, -33/32, 0, 13/64, 0, -1/8 ])), DiagonalMatrix([1,1,1])>;   // genus 3
    cover_data[{1,2,13,26}] := <HyperellipticCurve(Polynomial(Rationals(), [ -11, 38, -32 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,6}] := <HyperellipticCurve(Polynomial(Rationals(), [ 2197/32768, 0, -699/32768, 0, -25/32768, 0, -1/32768 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,3,13,39}] := <HyperellipticCurve(Polynomial(Rationals(), [ -11/4, 49/4, -291/16, 47/4, -27/4, 4 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,2,3,6}] := <HyperellipticCurve(Polynomial(Rationals(), [ -11/4, 27/4, -75/16, 19/8, -2 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,39}] := <HyperellipticCurve(Polynomial(Rationals(), [ -6591/262144, 0, -25/65536, 0, 387/131072, 0, 7/65536, 0, 1/262144 ])), DiagonalMatrix([1,1,1])>;   // genus 3
    cover_data[{1,6,13,78}] := <HyperellipticCurve(Polynomial(Rationals(), [ 1/4, -1/4, 1/16, -1/8 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,2,39,78}] := <HyperellipticCurve(Polynomial(Rationals(), [ 1/4, -3/4, 9/16, -1/4, 1/4 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,3,26,78}] := <HyperellipticCurve(Polynomial(Rationals(), [ 1, -2 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,78}] := <HyperellipticCurve(Polynomial(Rationals(), [ 1/8, 0, 9/64, 0, -1/32, 0, 1/64 ])), DiagonalMatrix([1,1,1])>;   // genus 2

    ws_data := AssociativeArray();
    return cover_data, ws_data;
end function;

procedure test_26_3()
    cover_data, ws_data := load_covers_and_ws_data_26_3();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(26, 3, cover_data, ws_data, curves);
    return;
end procedure;

test_26_3();
