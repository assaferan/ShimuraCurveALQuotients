// tests/_offline/X0_93_1.m -- RE-DERIVATION test for X_0^93(1).
//
// ⚠ WHY THIS FILE EXISTS. Before it, 93_1 was validated only against a COMMITTED FILE: the
// gy93_* block at the end of tests/GuoYangEquations.m compares data/models/models_93_1.m to
// Guo-Yang's published row, and ModelChecks checks it structurally. Neither RUNS THE PIPELINE, so
// "we can no longer produce this model" was invisible -- and 93_1 is exactly the base where that
// matters, because it regenerates ONLY since the vx fix (`n_oo`, BorcherdsForms.m:771). A silent
// regression of that fix would have left every committed artifact looking fine.
//
// ⚠ OFFLINE BECAUSE IT IS SLOW, NOT BECAUSE IT IS BROKEN: the model took 50927 s (14.1 h) on
// lovelace with no non-default flags (see the model file header). Far past any CI budget, and past
// GitHub's job limit, so it lives here with X0_10_19.m / X0_39_2.m / X0_87_1.m.
//     NORMALIZ_BIN=... magma -b filename:=tests/_offline/X0_93_1.m run_tests.m < /dev/null
// ⚠ NORMALIZ_BIN MUST BE SET. Without it a fresh polytope solve fails SILENTLY -- "no solutions"
// rather than an error (CLAUDE.md) -- which is precisely how X0_10_19 spent 84 min in CI verifying
// nothing.
//
// WHAT IS CHECKED, and at what strength:
//   [1] EXTERNAL, one comparison: the re-derived [1,93] quotient against GUO-YANG'S PUBLISHED
//       curve, y^2 = (3s^3-7s^2-3s-1)(3s^3+s^2-3s-9). ⚠ That is the TYPO-CORRECTED reading: the
//       journal prints `-3t`, and `-3t -> -3s` is the repair this repo DETERMINED (three other
//       plausible repairs give genus-2 curves that are NOT isomorphic; GuoYangEquations.m keeps all
//       three refutations live so the conclusion cannot decay into "some reading works"). The
//       journal version independently confirms it. So this entry is external evidence, not ours.
//   [2] DRIFT, automatically: test_AllEquationsAboveCoversSingleCurve now cross-checks EVERY
//       committed cover key against the same run, so [1,3] and [1,31] are compared too without
//       being named here.
//
// ⚠ WHAT IS **NOT** CHECKED, so nobody reads more into a pass than is there:
//   * NO INVOLUTIONS. ws_data is deliberately EMPTY. Guo-Yang's w_m for 93_1 are transcribed
//     nowhere in this repo, and they cannot join tests/GuoYangQuotientOracle.m's table because that
//     table encodes a single hyperelliptic f per base while 93_1's full curve is a genus-5 CRV PAIR.
//     Inventing them from our own `ws` would be circular (see the _gyinvol.m header).
//   * NOT THE FULL CURVE. The W={1} entry is that genus-5 CRV pair; the model cross-check SKIPS CRV
//     entries (the file records defining polynomials as strings but not the ambient weights) and
//     reports the skip in its count. A direct IsIsomorphic there is the 10 h+ regime measured at
//     26_3. So three of the four cover keys are pinned and the fourth is their fibre product --
//     strong, but not a full-curve proof. models_93_1.m's header says the same.
import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_93_1()
    P<s> := PolynomialRing(Rationals());

    // [Guo-Yang, Table A.1] for D = 93, N = 1, reading the printed `-3t` as `-3s`.
    gy_A := 3*s^3 - 7*s^2 - 3*s - 1;
    gy_B := 3*s^3 + s^2 - 3*s - 9;

    cover_data := AssociativeArray();
    // The second component is the `scales` matrix, used ONLY under manual_isomorphism/algebra_map,
    // which this test does not pass -- the identity is a placeholder, not a claim about coordinates.
    cover_data[{1,93}] := <HyperellipticCurve(gy_A*gy_B), IdentityMatrix(Rationals(), 3)>;

    ws_data := AssociativeArray();      // deliberately empty -- see the header
    return cover_data, ws_data;
end function;

procedure test_93_1()
    cover_data, ws_data := load_covers_and_ws_data_93_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(93, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_93_1();
