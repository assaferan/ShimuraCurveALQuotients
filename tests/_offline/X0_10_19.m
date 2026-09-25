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
    // ⇒ A REAL FIX HAS TO RE-DERIVE BOTH SIDES TOGETHER, and there is nothing external to anchor
    // them: `tests/GuoYangEquations.m` carries NO equation for this base (the journal's only word
    // on it is Remark 38, that X_0^10(19) is not hyperelliptic over Q), and neither the expected
    // curve nor the three ws_data matrices here cite a source.  So both are unattributed snapshots
    // of pipeline output, and re-pinning them against today's pipeline would make this test
    // compare the pipeline with itself.  What is worth preserving is the INVOLUTION check, which
    // needs the curve and the matrices to be re-derived in ONE consistent presentation.
    test_AllEquationsAboveCoversSingleCurve(10, 19, cover_data, ws_data, curves : manual_isomorphism);
    return;
end procedure;

test_10_19();