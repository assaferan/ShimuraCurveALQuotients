// tests/_offline/X0_111_1.m -- RE-DERIVATION test for X_0^111(1).
//
// ⚠ WHY. Like 93_1, this base was validated only against a COMMITTED FILE: GuoYangEquations.m
// compares data/models/models_111_1.m to Guo-Yang's published degree-16 curve, and ModelChecks
// checks it structurally. Neither RUNS THE PIPELINE, so "we can no longer produce this model" was
// invisible. These two were the only Guo-Yang bases with a stored full curve and no X0_* test.
//
// ⚠ STRONGER THAN THE 93_1 TEST, and worth saying why: 111_1's W={1} curve is HYPERELLIPTIC and
// Guo-Yang PRINT it, so the external anchor is the FULL CURVE itself rather than a quotient. At
// 93_1 the full curve is a genus-5 CRV pair, which forced that test to anchor on the [1,93] and
// [1,3] quotients and leave the top curve unpinned.
//
// OFFLINE because it is slow (the 93_1 sibling took 14.1 h; this base is genus 7 at W={1}).
//     NORMALIZ_BIN=... magma -b filename:=tests/_offline/X0_111_1.m run_tests.m < /dev/null
// ⚠ NORMALIZ_BIN MUST BE SET, or a fresh polytope solve fails SILENTLY -- "no solutions" rather
// than an error -- which is how X0_10_19 once spent 84 min in CI verifying nothing.
//
// WHAT IS CHECKED:
//   [1] EXTERNAL: the re-derived W={1} curve against Guo-Yang's published
//       y^2 = -(19x^8-44x^7-16x^6+55x^5+37x^4-55x^3-16x^2+44x+19)(x^8-3x^5-x^4+3x^3+1),
//       degree 16, genus 7. Transcribed from tests/GuoYangEquations.m, where the same polynomial
//       already validates the STORED model -- so this file adds the re-derivation, not the reading.
//   [2] DRIFT: the helper cross-checks every committed cover key against the same run, so
//       [1,3], [1,37] and [1,111] are compared too without being named here.
//
// ⚠ NOT CHECKED: involutions. ws_data is deliberately EMPTY -- Guo-Yang's w_m for 111_1 are
// transcribed nowhere in this repo, and taking them from our own `ws` would be circular (see the
// tests/_gyinvol.m header). A pass here says the CURVES are re-derived, not their labelling.
import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_111_1()
    P<x> := PolynomialRing(Rationals());

    // [Guo-Yang, Table A.1] D = 111, N = 1
    gy_f := -(19*x^8 - 44*x^7 - 16*x^6 + 55*x^5 + 37*x^4 - 55*x^3 - 16*x^2 + 44*x + 19)
            * (x^8 - 3*x^5 - x^4 + 3*x^3 + 1);

    cover_data := AssociativeArray();
    // second component: the `scales` matrix, used ONLY under manual_isomorphism/algebra_map, which
    // this test does not pass. The identity is a placeholder, not a claim about coordinates.
    cover_data[{1}] := <HyperellipticCurve(gy_f), IdentityMatrix(Rationals(), 3)>;

    ws_data := AssociativeArray();      // deliberately empty -- see the header
    return cover_data, ws_data;
end function;

procedure test_111_1()
    cover_data, ws_data := load_covers_and_ws_data_111_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(111, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_111_1();
