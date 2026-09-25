// ⚠ MOVED OUT OF CI 2026-09-06 -- IT WAS VERIFYING NOTHING THERE.
//
// This test passes LOCALLY with real comparisons, but in GitHub CI it made ZERO curve
// comparisons: `AllEquationsAboveCovers` produced the expected `W={1}` key with NO BASES, so the
// comparison loop never ran. It had been passing green on that basis, and cost 5016 s (84 min) per
// run to do it. The vacuity guard added to `test_AllEquationsAboveCoversSingleCurve` on
// 2026-09-06 is what exposed it.
//
// ROOT CAUSE: **CI never sets `NORMALIZ_BIN`.** `CLAUDE.md` is explicit that without it a fresh
// polytope solve fails SILENTLY -- "you get 'no solutions' rather than an error" -- so any cover
// needing a solve beyond the committed cache simply comes back empty. Locally the variable is set
// and the same cover is found.
//
// ⚠ MEASURED, so the scope is known: every OTHER X0_* job in CI reports full coverage
// (`X0_10_11` 1/1, `X0_26_1` 3/3, `X0_6_17` 1/1), so `10_19` is the only test affected. This is
// not a general CI collapse.
//
// ⇒ THE REAL FIX is to install Normaliz in CI and set `NORMALIZ_BIN`, after which this file should
// move back to `tests/`. Until then it lives here, where it runs meaningfully.
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
    // ⚠ WHAT REMAINS is mechanical but not done: turning that Mobius map into the 4x4 coordinate
    // matrix this test needs.  Guo-Yang's model has s = x^2 (compare their X/<w_2,w_95> equation
    // y^2 = -8s^3+57s^2-40s+16 with X/w_190's y^2 = -8x^6+57x^4-40x^2+16), so the remaining step is
    // to express our model's base coordinate in terms of our s and compose.
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