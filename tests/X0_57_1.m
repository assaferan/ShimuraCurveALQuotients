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
//
// ✅ INVOLUTIONS CHECKED (2026-09-07), and the W={1} CRV pair is now compared too. Guo-Yang
// publish, in THEIR coordinates:
//     w_19(s,x,y) = (s,x,-y)   w_57(s,x,y) = (s,-x,y)      [Guo-Yang's variables are (s,x,y)]
// Our model is a different presentation, so those do NOT carry over as written. They were
// TRANSPORTED: psi := construct_crv_isomorphism(our stored pair, Guo-Yang's pair) is computed from
// the two EQUATIONS alone, and the matrix recorded here is psi^-1 . w_GY . psi, which came out
// linear in the weighted coordinates. Script: tests/_gyinvol.m / the CRV variant beside it.
// ⚠ WHY THIS IS NOT CIRCULAR: the involutions are Guo-Yang's (external) and psi comes from
// equations, never from the pipeline's own `ws`. The harness then checks that the PIPELINE's
// involution labelled w_m matches Guo-Yang's w_m, so a labelling error is detectable.
// ⚠ psi is one element of a torsor under Aut; another choice conjugates all the transported
// involutions simultaneously, and the harness searches that same torsor, so it cannot cause a
// false verdict. Each matrix was verified to PRESERVE our curve (ideal membership) and to be an
// involution projectively -- `M^2 = identity` is the wrong test on a weighted ambient, where
// M^2 must act as (x_i) -> (lambda^{w_i} x_i).

function load_covers_and_ws_data_57_1()
    _<s> := PolynomialRing(Rationals());

    P3_1<x,y,s,z> := WeightedProjectiveSpace(Rationals(), [1,2,1,1]);

    cover_data := AssociativeArray();
    cover_data[{1}] := <Curve(P3_1, [ y^2 - 9*s^4 - 2*s^2*z^2 + 4*s*z^3 - z^4, x^2 + 144*s^2 - 180*s*z + 63*z^2 ]), DiagonalMatrix([1,1,1,1])>;   // genus 3, CRV pair
    cover_data[{1,3}] := <HyperellipticCurve(Polynomial(Rationals(), [ -7/9, 16/3, -110/9, 104/9, -95/9, 20, -16 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,19}] := <HyperellipticCurve(Polynomial(Rationals(), [ -63, 180, -144 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,57}] := <HyperellipticCurve(Polynomial(Rationals(), [ 1, -4, 2, 0, 9 ])), DiagonalMatrix([1,1,1])>;   // genus 1

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][19] := Matrix(4,4,[ 1,0,0,0,  0,-1,0,0,  0,0,1,0,  0,0,0,1 ]);
    ws_data[{1}][57] := Matrix(4,4,[ -1,0,0,0, 0, 1,0,0,  0,0,1,0,  0,0,0,1 ]);
    return cover_data, ws_data;
end function;

procedure test_57_1()
    cover_data, ws_data := load_covers_and_ws_data_57_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(57, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_57_1();
