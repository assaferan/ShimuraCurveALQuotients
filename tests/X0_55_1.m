import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

// tests/X0_55_1.m -- RE-DERIVATION test for X_0^55(1).
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

function load_covers_and_ws_data_55_1()
    _<s> := PolynomialRing(Rationals());

    cover_data := AssociativeArray();
    cover_data[{1,5}] := <HyperellipticCurve(Polynomial(Rationals(), [ -1/45375, -2/45375, -21/15125, -138/15125, -243/3025 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,11}] := <HyperellipticCurve(Polynomial(Rationals(), [ -1/5671875, 4/5671875, -32/1890625, -42/1890625, -1332/1890625, 216/378125, -2187/75625 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1}] := <HyperellipticCurve(Polynomial(Rationals(), [ -1/5808, -1/8712, 1/5808, 1/2178, -1/5808, -1/2178, 1/5808, 1/8712, -1/5808 ])), DiagonalMatrix([1,1,1])>;   // genus 3
    cover_data[{1,55}] := <HyperellipticCurve(Polynomial(Rationals(), [ 1/5, -6/5, 9 ])), DiagonalMatrix([1,1,1])>;   // genus 0

    ws_data := AssociativeArray();
    return cover_data, ws_data;
end function;

procedure test_55_1()
    cover_data, ws_data := load_covers_and_ws_data_55_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(55, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_55_1();
