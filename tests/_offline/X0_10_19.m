// ⚠ MOVED OUT OF CI 2026-09-06 -- IT WAS VERIFYING NOTHING THERE.
//
// In GitHub CI it made ZERO curve comparisons: `AllEquationsAboveCovers` produced the expected
// `W={1}` key with NO BASES, so the comparison loop never ran. It had been passing green on that
// basis, and cost 5016 s (84 min) per run to do it. The vacuity guard added to
// `test_AllEquationsAboveCoversSingleCurve` on 2026-09-06 is what exposed it.
//
// ⚠⚠ THE ORIGINAL DIAGNOSIS WAS WRONG, AND SO WAS "IT PASSES LOCALLY". This header used to open
// "This test passes LOCALLY with real comparisons" and blame CI's missing `NORMALIZ_BIN`:
//     "CI never sets NORMALIZ_BIN ... a fresh polytope solve fails SILENTLY ... Locally the
//      variable is set and the same cover is found."
// REFUTED BY DIRECT MEASUREMENT, 2026-09-27. The 2026-09-07 tree (17f407c) was checked out into a
// worktree and this test run against it LOCALLY, WITH `NORMALIZ_BIN` SET (3330 s). It failed:
//     NO EVIDENCE -- the test made ZERO curve comparisons, so it verified nothing.
//     Keys MATCHED but produced no bases, so there was nothing to compare against.
// -- exactly the CI symptom. So `NORMALIZ_BIN` was never the cause here, and the cover was NOT
// being found locally either.
//
// ⇒ THE ACTUAL ROOT CAUSE: the W={1} cover was NOT PRODUCIBLE BY ANY RUN until `EquationsByRebase`
// landed on 2026-09-08 (3684b75). models_10_19.m's `[1]` key was EMPTY from its 2026-07-08
// generation until e95cca4 (2026-09-09), whose message says the four empty keys
// ({1},{1,2},{1,10},{1,190}) were "populated, unlocked by EquationsByRebase" -- and that stage
// "adopts ONLY keys that were empty", so before it existed this key could not be filled at all.
//
// ⇒ SO THIS TEST HAS NEVER VERIFIED ITS CURVE, AND ITS RED STATUS IS NOT A REGRESSION. Timeline:
// born 2025-12-15 in a commit saying outright "we can't yet find a curve"; switched to
// `manual_isomorphism` three days later; vacuously GREEN until the vacuity guard (9653f12);
// RED-because-vacuous until 2026-09-08; RED-because-the-matrix-is-not-a-map once the cover finally
// became producible. **The pinned matrix has never once been exercised successfully** -- it was
// written speculatively and nothing has ever checked it.
// ⇒ DO NOT go looking for the commit that "broke" this. A comparison started running for the first
// time and failed immediately. Same shape as the CI vacuity: the instrument, not the curve.
//
// ⚠ MEASURED, so the scope is known: every OTHER X0_* job in CI reports full coverage
// (`X0_10_11` 1/1, `X0_26_1` 3/3, `X0_6_17` 1/1), so `10_19` is the only test affected. This is
// not a general CI collapse.
//
// ⇒ Installing Normaliz in CI is still worth doing, but it will NOT make this file green and it is
// no longer "the real fix" for this test -- the transport below is.
//
import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_10_19()
    P3<x,y,z,s> := WeightedProjectiveSpace(Rationals(), [1,3,1,1]);
     // D = 10, N = 19   (this said "D = 82" until 2026-09-25 -- a copy-paste from another file)
    cover_data := AssociativeArray();
    cover_data[{1}] := <Curve(P3, [y^2 + 8*x^6 - 57*x^4*s^2 + 40*x^2*s^4 - 16*s^6, z^2 - 5*x^2 + 32*s^2]), Matrix([[0,0,1,0],[0,1/8,0,0],[-1/8,0,0,-1/8],[-1/4,0,0,0]])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2] := DiagonalMatrix([-1,1,1,1]);
    ws_data[{1}][5] := DiagonalMatrix([1,-1,-1,1]);
    ws_data[{1}][38] := DiagonalMatrix([1,-1,1,1]);

    return cover_data, ws_data;
end function;

procedure test_10_19()
    cover_data, ws_data := load_covers_and_ws_data_10_19();
    curves := GetHyperellipticCandidates();
    // ⚠⚠ THIS TEST IS RED, measured 2026-09-25 (3500 s), and the cause is NOT what the header
    // above says.  The header's claim that it "passes LOCALLY with real comparisons" dates from
    // 2026-09-06 and is stale; run with NORMALIZ_BIN set it fails at BorcherdsProducts.m:93 with
    // "Polynomials do not define a map into the codomain" -- the pinned matrix in cover_data[{1}]
    // is no longer a map.
    //
    // ⚠ DO NOT "FIX" IT BY DROPPING manual_isomorphism.  That was tried on 2026-09-25 and cost an
    // hour: the failure merely moves to "the produced cover matches NONE of the 1 acceptable
    // curve(s)".  The CRV branch cannot bridge these two, and the reason is visible without
    // running anything -- the two pairs present the curve over DIFFERENT intermediate quotients:
    //
    //     committed model (data/models/models_10_19.m, key [1]):
    //         y^2 + 1/320000 s^6 + ... + 475/2097152 z^6     sextic in (s,z)
    //         x^2 + 1/128 s^2 - 125/2048 z^2                 conic variable is x
    //     expected here:
    //         y^2 + 8x^6 - 57x^4 s^2 + 40x^2 s^4 - 16 s^6    sextic in (x,s)
    //         z^2 - 5x^2 + 32 s^2                            conic variable is z
    //
    // The roles of x and z are swapped.  `construct_crv_isomorphism` declines exactly this case,
    // so the pinned matrix was the only bridge between the two presentations -- and it is the
    // bridge that rotted, not the matrix's arithmetic.
    //
    // ✅ THE EXPECTED CURVE IS PUBLISHED -- IT IS AN ORACLE, NOT A SNAPSHOT.  Guo-Yang,
    // Compositio Math. 153 (2017), EXAMPLE 37 is exactly X = X_0^10(19).  With s the Hauptmodul of
    // X/W_{10,19} normalised by s(tau_-8) = 0, s(tau_-40) = infinity, s(tau_-3) = 1, they print
    //
    //     X/<w_2, w_95> :  y^2 = -8s^3 + 57s^2 - 40s + 16     (Cremona E190A1)
    //     X/w_190       :  y^2 = -8x^6 + 57x^4  - 40x^2 + 16
    //
    // and the sextic in cover_data[{1}] above, homogenised, IS that second equation verbatim
    // (y^2 = -8x^6 + 57x^4 s^2 - 40x^2 s^4 + 16 s^6, i.e. s = 1).  So this key is anchored in the
    // published text and MUST NOT be replaced by the committed model's presentation -- doing so
    // would throw the oracle away and leave the pipeline compared with itself.
    //
    // ⚠ TRANSCRIPTION GAP: tests/GuoYangEquations.m carries NO row for 10_19, because its rows come
    // from Tables A.1/A.2 and this base's equations are in the BODY, in Example 37.  Remark 38
    // (X_0^10(19) is not hyperelliptic over Q) is why it is not in the hyperelliptic tables at all.
    // Same shape as the GuoYangTable1 find: ask of every source what else it publishes that we do
    // not read.
    //
    // ✅✅ THE TRANSCRIPTION IS VERIFIED -- THERE IS NO TYPO HERE.  Checked 2026-09-27 against the
    // JOURNAL PDF itself (Compositio 153 (2017), Example 37, printed page 28), verbatim: both
    // equations, all four coefficients of each, the three involutions
    //     w_2 (x,y,z) = (-x, y, z)   w_5 (x,y,z) = (x,-y,-z)   w_19(x,y,z) = (-x,-y, z)
    // and the normalisation s(tau_-8) = 0, s(tau_-40) = infinity, s(tau_-3) = 1.  ⇒ DO NOT spend
    // time hunting a transcription error to explain this file's RED status; the expected curve is
    // right and the failure is purely one of PRESENTATION (see the structural analysis below).
    //
    // Corroborated independently of the PDF, so the check does not rest on one reading:
    //   * y^2 = -8s^3+57s^2-40s+16 IS Cremona 190a1, conductor 190 = D*N -- and Guo-Yang assert
    //     E190A1 themselves, so their equation is confirmed too, not merely our copy of it.
    //   * the conic z^2 = 5x^2-32 is POINTLESS over Q and its quaternion algebra ramifies at
    //     exactly {2,5} = the primes dividing D = 10.  (Remark 38 is precisely this point.)
    //   * the pair has genus 5 = GenusShimuraCurve(10,19).
    //   * tests/GuoYangQuotients_10_19.m derives 15 quotients from these same two polynomials and
    //     matches them against the pipeline: green, 0.28 s, 0 keys empty.  A coefficient slip would
    //     break that oracle.
    // ✅✅ AND THE w_m LABELS ARE VERIFIED TOO (2026-09-27) -- UNLIKE 39_2, THIS BASE CAN BE
    // ARBITRATED.  The reason is that Example 37 states its CM DISCRIMINANTS in prose, which is
    // exactly the information a bare Table A.1 row withholds; that is why 39_2 needed CM values
    // computed from scratch and this one does not.  NumFixedPointsByCMOrder(10,19,m) -- class
    // numbers and Ogg's local embedding numbers, no Guo-Yang input at all -- gives
    //     w_2  : 4 fixed pts, ALL of disc  -8        w_5, w_19, w_95 : NO fixed points
    //     w_10 : 4 fixed pts, ALL of disc -40
    //     w_38 : 12 fixed pts, ALL of disc -152
    //     w_190: 4 fixed pts, ALL of disc -760
    // The discriminant is -4m and is INJECTIVE here, so it identifies the involution.  Checking
    // their two ramification statements against it:
    //   * "X/<w_5,w_38> -> X/W_{10,19} is ramified at the CM point of disc -8 and the CM point of
    //     disc -40".  That cover is quotient by the coset w_2*{1,5,38,190} = {w_2,w_10,w_19,w_95},
    //     whose fixed points are exactly w_2 (-8) and w_10 (-40); w_19 and w_95 have none. ✓
    //     ⚠ Note w_190 is NOT in that coset (w_2*w_190 = w_95, not the other way round) -- and that
    //     is load-bearing: w_190 HAS 4 fixed points, so if it were in the coset their claim would
    //     be false.  Getting this composition backwards is the easy mistake here.
    //   * "X/w_38 -> X/<w_5,w_38> is ramified at THE TWO CM points of disc -760".  That cover is
    //     quotient by the coset {w_5, w_190}: w_5 has none, w_190 has 4 on X, which become 2 on
    //     X/w_38. ✓ Both the discriminant AND the count.
    //   * their normalisation s(tau_-8) = 0, s(tau_-40) = infinity puts the Hauptmodul's zero and
    //     pole precisely at those two branch points. ✓
    // ⇒ So the labels are pinned by the ONE instrument that can pin them.  Corroborating, and
    // showing the existing oracle is not vacuous: all seven Ogg-predicted quotient genera equal the
    // genus of the curve Guo-Yang's own matrices produce, and within each genus class the quotients
    // are pairwise NON-ISOMORPHIC ({w_190,w_2,w_10} at genus 2, {w_5,w_19,w_95} at genus 3).  Genus
    // alone cannot separate those triples -- the curves can -- so GuoYangQuotients_10_19.m's 15
    // comparisons WOULD have caught a swap.  Internally too: ws_data[38] below is w_2*w_19, and the
    // composite w_190 = (x,y,-z) leaves y^2 = f(x) of genus 2, the equation they print for X/w_190.
    // ⚠ One slip in THEIR text, harmless and already read correctly here: they print "s(tau_760)"
    // where the argument requires discriminant -760 ("the two CM points of discriminant -760").
    //
    // ⇒ SO THE FIX IS A TRANSPORT, not a re-pin: carry the published presentation onto the
    // pipeline's current one and record the resulting matrix, exactly as tests/_gyinvol.m does for
    // involutions.  The three ws_data matrices below must be transported through the SAME map, or
    // they will be checked against the wrong object.
    //
    // ✅ THE HAUPTMODUL BRIDGE IS FOUND AND VERIFIED (2026-09-25).  Example 37 states Guo-Yang's
    // normalisation outright -- s is the Hauptmodul of X/W_{10,19} with s(tau_-8) = 0,
    // s(tau_-40) = infinity, s(tau_-3) = 1 -- so the two presentations differ by the Mobius map
    // between their Hauptmodul normalisations, and that map is DETERMINED, not searched for.
    // ValuesAtCMPoints gives ours:
    //
    //     disc      -8     -40     -3        -760
    //     ours       0       1      infinity  32/27
    //     Guo-Yang   0       infinity  1      32/5
    //
    // The first three force  s_GY = s_ours / (s_ours - 1).  The FOURTH IS A CHECK, NOT AN INPUT:
    // that map sends 32/27 to 32/5, which is exactly what Example 37 prints for s(tau_-760).
    // So the bridge is confirmed by a value it was not fitted to.
    //
    // ✅ THE TWO MODELS ARE THE SAME CURVE -- measured 2026-09-27, so a pipeline defect here is
    // RULED OUT.  Counting points of each (sextic, conic) pair over F_p and INCLUDING the fibre
    // over t = infinity (both polynomials have even degree, so that fibre is governed by their
    // leading coefficients) gives identical totals at all 15 good primes p <= 61:
    //     8 16 8 20 20 16 40 24 60 32 44 52 60 72 56
    // ⚠ The affine count alone is NOT an invariant and says the opposite -- it differs at 6 of 13
    // primes, always by exactly +-4, purely from the boundary.  Do not repeat that: count the
    // fibre at infinity or the comparison is presentation-dependent.
    //
    // ⚠⚠ WHY THE TRANSPORT IS STILL NOT DONE, and why three attempts failed the same way.  The
    // obstruction is structural, not arithmetic: the two pairs present the curve over DIFFERENT
    // intermediate quotients (ours has conic variable x over the (z,s)-line; Guo-Yang's has conic
    // variable z over the x-line).  So there is no Mobius map between the two BASE lines to find,
    // and any parametrisation assuming one is empty by construction.  Measured: solving
    //     nu^2 * P(z) = -8u^6 + 57u^4 v^2 - 40u^2 v^4 + 16v^6,   mu^2 * Q(z) = 5u^2 - 32v^2
    // for u = az+b, v = cz+d has NO rational solution under any pin of a, b, c or d -- including
    // with mu^2, nu^2 as the unknowns rather than mu, nu, so a non-square scaling is not the issue.
    // ⚠⚠ base_label WAS SWEPT, AND IT DOES NOT RESOLVE THIS BASE (2026-09-27).  On a CRV base the
    // FIRST move is `magma -b Dd:=10 Nn:=19 tests/_basesweep.m` -- the base decides which V_4 the
    // pair presents, and it is what found 5394 for 14_3 and 8103 for 26_3 in seconds after exactly
    // this symptom.  Here it comes back negative:
    //
    //     W={1} is label 4089; default base 4091;  candidates 4097, 4102, 4103
    //       4097  no pair produced
    //       4102  genus 5, construct_crv_isomorphism to GY: FALSE
    //       4103  no pair produced
    //
    // So this base is in the 21_2 category, not the 14_3 / 26_3 one: the candidate pool is three,
    // two yield no pair, and the one that does is not Guo-Yang's V_4.  ⇒ Do not re-run the sweep
    // expecting a different answer; the base run alone is 3371 s.
    //
    // ⇒ WHAT IS LEFT is the slow route, as at 21_2: a general weight-respecting map on P(1,3,1,1)
    // with the sextic identity holding MODULO the conic.  ⚠ Three parametrisations that assume a
    // Mobius between the two BASE lines are already known to come back empty -- that assumption is
    // what the differing V_4 rules out, so do not spend a fourth on it.
    //
    // ⚠⚠ TWO FAILED ATTEMPTS BEFORE THIS, BOTH THE SAME ERROR -- assuming a coordinate
    // correspondence instead of deriving one.  (1) Solving for a linear map required the sextic
    // identity to hold in the polynomial ring, when the two equations define the curve TOGETHER and
    // it need only hold MODULO the conic.  (2) That same slip made me argue the published z must be
    // proportional to the pipeline's conic variable X; the conic isomorphism shows the published
    // **s** is, while z corresponds to the pipeline's S.  The resulting "no solution, dimension -1"
    // was an artifact of the ansatz and says nothing about the curves.  Independently: the two
    // conics ARE in the same class -- both pointless, both ramified at {2,5} -- so the pipeline is
    // not building a wrong object here.
    test_AllEquationsAboveCoversSingleCurve(10, 19, cover_data, ws_data, curves : manual_isomorphism);
    return;
end procedure;

test_10_19();