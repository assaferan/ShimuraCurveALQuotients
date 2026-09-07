import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

// tests/X0_14_3.m -- RE-DERIVATION test for X_0^14(3).
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

function load_covers_and_ws_data_14_3()
    _<s> := PolynomialRing(Rationals());

    cover_data := AssociativeArray();
    cover_data[{1,2,7,14}] := <HyperellipticCurve(Polynomial(Rationals(), [ -11, 104, -128 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,42}] := <HyperellipticCurve(Polynomial(Rationals(), [ -28, 0, 88, 0, 4 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,3,7,21}] := <HyperellipticCurve(Polynomial(Rationals(), [ -176, 2368, -3776, 1024 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,3}] := <HyperellipticCurve(Polynomial(Rationals(), [ 63, 0, -184, 0, -53, 0, -2 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,7}] := <HyperellipticCurve(Polynomial(Rationals(), [ 768, 2304, 4160, 45440/9, 4160, 2304, 768 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,6}] := <HyperellipticCurve(Polynomial(Rationals(), [ 3087/2, 0, 577/2, 0, 17/2, 0, -1/2 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,2,3,6}] := <HyperellipticCurve(Polynomial(Rationals(), [ -11, 236, -1420, 1952, -512 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,2,21,42}] := <HyperellipticCurve(Polynomial(Rationals(), [ 64, -768, 256 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,6,14,21}] := <HyperellipticCurve(Polynomial(Rationals(), [ -11, 16 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,3,14,42}] := <HyperellipticCurve(Polynomial(Rationals(), [ 1, -8 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,21}] := <HyperellipticCurve(Polynomial(Rationals(), [ -343, 0, -26, 0, 1 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,6,7,42}] := <HyperellipticCurve(Polynomial(Rationals(), [ 64, -1280, 6400, -2048 ])), DiagonalMatrix([1,1,1])>;   // genus 1

    ws_data := AssociativeArray();
    return cover_data, ws_data;
end function;

procedure test_14_3()
    cover_data, ws_data := load_covers_and_ws_data_14_3();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(14, 3, cover_data, ws_data, curves);
    return;
end procedure;

test_14_3();
