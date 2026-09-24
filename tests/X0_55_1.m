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
//
// ✅ INVOLUTIONS CHECKED (2026-09-07). Guo-Yang publish, in THEIR coordinates:
//     w_5(x,y)  = (-1/x, y/x^4)          w_55(x,y) = (x, -y)
// Our model is a different presentation, so those matrices do NOT carry over as written. They were
// TRANSPORTED: psi := IsIsomorphic(our stored curve, Guo-Yang's curve) is computed from the two
// EQUATIONS alone, and the involution recorded here is psi^-1 . w_GY . psi, which came out linear
// in the weighted coordinates and so is expressible as a matrix.
// ⚠ WHY THIS IS NOT CIRCULAR: the involutions are Guo-Yang's (external), and psi is derived from
// equations, never from the pipeline's own `ws`. What the harness then checks is that the
// PIPELINE's involution labelled w_m matches Guo-Yang's w_m under some identification -- so an
// error in the pipeline's LABELLING is detectable, which is the whole point.
// ⚠ psi is one element of a torsor under Aut, and another choice would conjugate all the
// transported involutions simultaneously. That is harmless here because the harness searches that
// same torsor (see BorcherdsProducts.m), so the choice cannot cause a false verdict either way.
// Each matrix was verified to be an involution OF OUR CURVE and to equal the transported map.

function load_covers_and_ws_data_55_1()
    _<s> := PolynomialRing(Rationals());

    cover_data := AssociativeArray();
    cover_data[{1,5}] := <HyperellipticCurve(Polynomial(Rationals(), [ -1/45375, -2/45375, -21/15125, -138/15125, -243/3025 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,11}] := <HyperellipticCurve(Polynomial(Rationals(), [ -1/5671875, 4/5671875, -32/1890625, -42/1890625, -1332/1890625, 216/378125, -2187/75625 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1}] := <HyperellipticCurve(Polynomial(Rationals(), [ -1/5808, -1/8712, 1/5808, 1/2178, -1/5808, -1/2178, 1/5808, 1/8712, -1/5808 ])), DiagonalMatrix([1,1,1])>;   // genus 3
    cover_data[{1,55}] := <HyperellipticCurve(Polynomial(Rationals(), [ 1/5, -6/5, 9 ])), DiagonalMatrix([1,1,1])>;   // genus 0

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][5]  := Matrix(3,3,[ 0, 0, 1, 0, 1, 0, -1, 0, 0 ]);
    ws_data[{1}][55] := Matrix(3,3,[ 1, 0, 0, 0, -1, 0, 0, 0, 1 ]);
    return cover_data, ws_data;
end function;

procedure test_55_1()
    cover_data, ws_data := load_covers_and_ws_data_55_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(55, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_55_1();
