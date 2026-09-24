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
// ✅ INVOLUTIONS CHECKED (2026-09-07). Guo-Yang publish, in THEIR coordinates:
//     w_2(x,y) = (-1/x, y/x^4)   w_14(x,y) = (x, -y)
//     w_35(x,y) = ((x+2)/(2x-1), +25y/(2x-1)^4)
//
// ⚠ THE +25 IS A CORRECTION, AND THIS TEST IS WHAT DETERMINES IT. The journal contradicts itself:
// Example 36's text prints +25z/(2x-1)^4, while the TABLE row for X_0^14(5) prints -25y/(2x-1)^4
// (so does arXiv v1's table). Both maps are involutions of the curve, so neither is refutable by
// inspection -- but they differ by w_14, and 35*14/gcd(35,14)^2 = 10, so the two readings are
// w_35 and w_10. The pipeline adjudicates: with +25 labelled w_35 the test passes, with -25
// labelled w_35 it FAILS on the labelling, and with -25 labelled w_10 it passes again. The TABLE
// carries the sign error; Example 36 is right. This is the second Guo-Yang table typo the pipeline
// has determined (models_93_1.m's `-3t` -> `-3s` is the first).
// Both readings are kept below, each under its OWN label, so the run makes 4 involution
// comparisons and a future change cannot quietly relabel one into the other.
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
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2]  := Matrix(3,3,[ -3, 0, 2, 0, 1, 0, -5, 0, 3 ]);
    ws_data[{1}][14] := Matrix(3,3,[ 1, 0, 0, 0, -1, 0, 0, 0, 1 ]);
    ws_data[{1}][35] := Matrix(3,3,[ -5, 0, 2, 0,  25, 0, -10, 0, 5 ]);   // Example 36's +25
    ws_data[{1}][10] := Matrix(3,3,[ -5, 0, 2, 0, -25, 0, -10, 0, 5 ]);   // the table's -25 IS w_10
    return cover_data, ws_data;
end function;

procedure test_14_5()
    cover_data, ws_data := load_covers_and_ws_data_14_5();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(14, 5, cover_data, ws_data, curves);
    return;
end procedure;

test_14_5();
