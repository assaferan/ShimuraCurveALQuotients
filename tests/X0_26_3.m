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
// ✅ INVOLUTIONS CHECKED (2026-09-07), and the W={1} CRV pair is now compared too. Guo-Yang
// publish, in THEIR coordinates:
//     w_2(x,y,z) = (-x,-y,-z)   w_3(x,y,z) = (x,-y,-z)   w_26(x,y,z) = (x,-y,z)
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

function load_covers_and_ws_data_26_3()
    _<s> := PolynomialRing(Rationals());

    P3_1<x,y,s,z> := WeightedProjectiveSpace(Rationals(), [1,3,1,1]);

    cover_data := AssociativeArray();
    cover_data[{1}] := <Curve(P3_1, [ y^2 - 1/64*s^6 + 1/32*s^4*z^2 - 9/64*s^2*z^4 - 1/8*z^6, x^2 + 8*s^2 + 3*z^2 ]), DiagonalMatrix([1,1,1,1])>;   // genus 5, CRV pair
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
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2]  := Matrix(4,4,[ -1,0,0,0,  0,-1,0,0,  0,0,-1,0,  0,0,0,1 ]);
    ws_data[{1}][3]  := Matrix(4,4,[ -1,0,0,0,  0,-1,0,0,  0,0, 1,0,  0,0,0,1 ]);
    ws_data[{1}][26] := Matrix(4,4,[  1,0,0,0,  0,-1,0,0,  0,0, 1,0,  0,0,0,1 ]);
    return cover_data, ws_data;
end function;

procedure test_26_3()
    cover_data, ws_data := load_covers_and_ws_data_26_3();
    curves := GetHyperellipticCandidates();
    // ⚠ base_label := 8103 IS REQUIRED, not cosmetic: models_26_3.m deliberately stores the V_4
    // that Guo-Yang use, which a DEFAULT run does not produce (it gives a different, equally valid
    // one). Without it the W={1} CRV pair the pipeline emits is a genuinely different presentation
    // and the isomorphism assertion fails. See models_26_3.m's header.
    // ⚠ model_drift_ok: this test pins a NON-ZERO base_label, and AllEquationsAboveCovers gates
    // EquationsByRebase on `base_label eq 0` (EquationsCovers.m:1061). So it cannot reproduce the
    // cover keys that the rebase FILLED on a default run -- [1,2] and [1,13], which
    // data/models/models_26_3.m records as "previously EMPTY ... now filled, unlocked by
    // EquationsByRebase". MISSING keys only: a key this test DOES produce must still be the
    // committed curve, and model_drift_ok does not silence that.
    // ⚠ THE CONTROL GROUP is what makes this a diagnosis rather than an excuse: 14_3, 21_2 and
    // 6_17 also pin a base_label and all three PASS -- 14_3's empties were fixed by the coprime
    // filter flip, not the rebase, and the other two never had any. The gate costs exactly the
    // rebase-filled keys and nothing else.
    test_AllEquationsAboveCoversSingleCurve(26, 3, cover_data, ws_data, curves : base_label := 8103);
    return;
end procedure;

test_26_3();
