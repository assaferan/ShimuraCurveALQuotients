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
//   [3] ✅ THE FULL CURVE AND ITS INVOLUTIONS, both added 2026-09-23 -- see the two blocks in
//       load_covers_and_ws_data_93_1 below. Guo-Yang publish 93_1 as a PAIR, so the W={1} entry is
//       built from THEIR equations in P(1,3,1,1) and compared through the helper's
//       construct-the-isomorphism branch rather than a 10 h+ IsIsomorphic; their w_3 and w_31 come
//       from the same journal row.
//
// ⚠ WHAT WAS **NOT** CHECKED BEFORE THAT, kept because it explains the file's shape and because
// both statements were TRUE WHEN WRITTEN and are now retired:
//   * NO INVOLUTIONS -- "Guo-Yang's w_m for 93_1 are transcribed nowhere in this repo". They are
//     now: journal, Table A.1, printed page 34. The separate point stands, that they cannot join
//     tests/GuoYangQuotientOracle.m's table, which encodes a single hyperelliptic f per base while
//     93_1's full curve is a genus-5 CRV pair. Taking them from our own `ws` would still be
//     circular (see the _gyinvol.m header); they are taken from the paper.
//   * NOT THE FULL CURVE -- three of four cover keys pinned, the fourth their fibre product. The
//     model cross-check does still SKIP our stored CRV entry (models_93_1.m records defining
//     polynomials as strings without the ambient weights, and reports the skip in its count); what
//     changed is that the expected curve no longer has to come from that stored entry at all,
//     since Guo-Yang printed the pair. ⚠ models_93_1.m's header still says "not a full-curve
//     proof" and should be updated once this file has actually been RUN green.
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

    // ---- THE FULL CURVE, W = {1} -----------------------------------------------------------
    // ✅ ADDED 2026-09-23. This was the last X0_* test with NO top-curve key at all.
    //
    // Guo-Yang publish 93_1 as a PAIR, so the expected curve is built from THEIR equations rather
    // than from our stored strings -- external evidence, and strictly stronger than comparing a
    // committed file against itself:
    //     y^2 = (3s^3-7s^2-3s-1)(3s^3+s^2-3s-9)        x^2 = -4s^2-6s-9          [genus 5]
    // Homogenised by z, in P(1,3,1,1) with coordinates (x,y,z,s): y has weight 3 because its
    // equation is degree 6, exactly as 21_2's does (X0_82_1.m's y has weight 2 for a degree-4 one).
    //
    // ⚠ NO manual_isomorphism, DELIBERATELY, and this is the whole point of the entry. The helper
    // has a branch for exactly this shape -- both sides non-CrvHyp with 2 defining polynomials --
    // which CONSTRUCTS the isomorphism (Mobius map from the hyperelliptic quotient, constant
    // squares both sides, certified by IsIsomorphism) in hundredths of a second, instead of the
    // 10 h+ IsIsomorphic measured at 26_3. Pinning a matrix by hand would also work but is the
    // brittle option: tests/BorcherdsProducts.m records that under CMNONCOPRIME=1 the pipeline
    // re-presents 10_13's curve and the hardcoded map stops being a map at all. Only two tests
    // pin a matrix (X0_82_1.m, X0_10_19.m); this one should not become the third.
    //
    // ⚠ RESIDUAL RISK, stated because it is not checkable without a 14 h run: the construct branch
    // fires only if the RE-DERIVED W={1} also arrives as a non-CrvHyp with exactly 2 defining
    // polynomials. If the pipeline ever presents it otherwise, the helper falls through to a plain
    // IsIsomorphic and this test becomes the 10 h+ case rather than failing. If this file's runtime
    // suddenly jumps, that is the cause -- not a regression in the mathematics.
    P3<x,y,z,s3> := WeightedProjectiveSpace(Rationals(), [1,3,1,1]);
    gy_A3 := 3*s3^3 - 7*s3^2*z - 3*s3*z^2 -   z^3;
    gy_B3 := 3*s3^3 +   s3^2*z - 3*s3*z^2 - 9*z^3;
    cover_data[{1}] := <Curve(P3, [ y^2 - gy_A3*gy_B3,  x^2 + 4*s3^2 + 6*s3*z + 9*z^2 ]),
                        IdentityMatrix(Rationals(), 4)>;

    // ---- ATKIN-LEHNER INVOLUTIONS ----------------------------------------------------------
    // From the JOURNAL, Table A.1, printed page 34 -- the same row as the equation above:
    //     w_3 (s,x,y) = (s, -x, -y)        w_31(s,x,y) = (s,  x, -y)
    // ✅ This RETIRES this file's own stated blocker. The header used to say Guo-Yang's w_m for
    // 93_1 "are transcribed nowhere in this repo"; that was true when written and is not now.
    // They are in GUO-YANG's coordinates, which is right here because cover_data[{1}] above is
    // THEIR curve, so no transport is needed (unlike 39_2, where our model sits elsewhere).
    //
    // ⚠ VERIFIED BY SUBSTITUTION INTO THE DEFINING IDEAL, NOT by IsIsomorphism. On this weighted
    // ambient IsIsomorphism reports FALSE for both of these although both manifestly preserve the
    // equations -- the toric-ambient breakage of Magma #123, which hits exactly our P(1,w,1,1)
    // pairs. Substituting (x,y,z,s) -> (-x,-y,z,s) and (x,y,z,s) -> (x,-y,z,s) leaves BOTH
    // defining polynomials literally unchanged, while a control substitution x -> x+z breaks them.
    // ⇒ A FAILING CHECK NEEDS ITS OBJECT CHECKED AS MUCH AS A PASSING ONE: taken at face value that
    // IsIsomorphism result would have rejected two correct matrices.
    // Precedent that the helper copes with diagonal ws_data on a weighted ambient: X0_82_1.m does
    // exactly this in P(1,2,1,1) and is green.
    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][3]  := DiagonalMatrix([-1,-1,1,1]);
    ws_data[{1}][31] := DiagonalMatrix([ 1,-1,1,1]);
    return cover_data, ws_data;
end function;

procedure test_93_1()
    cover_data, ws_data := load_covers_and_ws_data_93_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(93, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_93_1();
