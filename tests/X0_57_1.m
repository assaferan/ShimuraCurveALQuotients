import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

// tests/X0_57_1.m -- RE-DERIVATION test for X_0^57(1).
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

function load_covers_and_ws_data_57_1()
    _<s> := PolynomialRing(Rationals());

    cover_data := AssociativeArray();
    cover_data[{1,3}] := <HyperellipticCurve(Polynomial(Rationals(), [ -7/9, 16/3, -110/9, 104/9, -95/9, 20, -16 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,19}] := <HyperellipticCurve(Polynomial(Rationals(), [ -63, 180, -144 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,57}] := <HyperellipticCurve(Polynomial(Rationals(), [ 1, -4, 2, 0, 9 ])), DiagonalMatrix([1,1,1])>;   // genus 1

    ws_data := AssociativeArray();
    return cover_data, ws_data;
end function;

procedure test_57_1()
    cover_data, ws_data := load_covers_and_ws_data_57_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(57, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_57_1();
