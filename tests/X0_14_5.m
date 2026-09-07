import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

// tests/X0_14_5.m -- RE-DERIVATION test for X_0^14(5).
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
// ⚠ THIS TEST IS WEAKER THAN THE HAND-WRITTEN ONES: `ws_data` IS EMPTY, so it makes ZERO
// involution comparisons. It verifies that each cover is ISOMORPHIC to the stored curve, but not
// that the Atkin-Lehner involutions correspond -- and the involutions are what make these QUOTIENT
// models rather than merely curves. 23 of the 34 X0_*.m tests do check them; the ones generated on
// 2026-09-07 (this file among them) do not.
// ⇒ Closing that needs involution matrices IN OUR MODEL'S COORDINATES. Taking them from the
// pipeline's own `ws` output would be circular; deriving them from Guo-Yang's published
// involutions is independent but is per-base work. Recorded rather than silently accepted.

function load_covers_and_ws_data_14_5()
    _<s> := PolynomialRing(Rationals());

    cover_data := AssociativeArray();
    cover_data[{1,2,7,14}] := <HyperellipticCurve(Polynomial(Rationals(), [ -11, 216, -1024 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,70}] := <HyperellipticCurve(Polynomial(Rationals(), [ -7/16, 0, 11/8, 0, 1/16 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1}] := <HyperellipticCurve(Polynomial(Rationals(), [ -14375/4096, -22875/1024, -124175/2048, -23565/256, -352037/4096, -3249/64, -38075/2048, -3963/1024, -1439/4096 ])), DiagonalMatrix([1,1,1])>;   // genus 3
    cover_data[{1,2,5,10}] := <HyperellipticCurve(Polynomial(Rationals(), [ -100/289, 0, 57/289, 0, -8/289 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,35}] := <HyperellipticCurve(Polynomial(Rationals(), [ -7/4096, 0, -181/2048, 0, 1/4096 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,2,35,70}] := <HyperellipticCurve(Polynomial(Rationals(), [ 1, -12, 4 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,5,14,70}] := <HyperellipticCurve(Polynomial(Rationals(), [ 1, -8 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,10,14,35}] := <HyperellipticCurve(Polynomial(Rationals(), [ -11, 128 ])), DiagonalMatrix([1,1,1])>;   // genus 0

    ws_data := AssociativeArray();
    return cover_data, ws_data;
end function;

procedure test_14_5()
    cover_data, ws_data := load_covers_and_ws_data_14_5();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(14, 5, cover_data, ws_data, curves);
    return;
end procedure;

test_14_5();
