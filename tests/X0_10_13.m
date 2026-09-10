import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_10_13()
    P3<x,y,z,s> := WeightedProjectiveSpace(Rationals(), [1,2,1,1]);
     // D = 82
    cover_data := AssociativeArray();
    cover_data[{1}] := <Curve(P3, [y^2 - 5*x^4 + 74*x^2*s^2 - 325*s^4, z^2 + 2*x^2 + 25*s^2]), Matrix([[0,0,1,0],[0,1/8,0,0],[-1/8,0,0,-1/8],[-1/4,0,0,0]])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2] := DiagonalMatrix([1,-1,-1,1]);
    // Here there is a mistake in [GY] - this can be checked by looking at
    // the fixed point under w_{10}. Using their w_5, we get that it is defined over Q(sqrt{13}),
    // while it should be defined over Q(sqrt{5}, sqrt{-2})
    ws_data[{1}][5] := DiagonalMatrix([-1,-1,1,1]);
    ws_data[{1}][65] := DiagonalMatrix([1,-1,1,1]);

    return cover_data, ws_data;
end function;

procedure test_10_13()
    cover_data, ws_data := load_covers_and_ws_data_10_13();
    curves := GetHyperellipticCandidates();
    // manual_isomorphism DROPPED 2026-09-07: the helper now CONSTRUCTS the isomorphism for
    // CRV pairs (tests/_crviso.m) instead of calling IsIsomorphic, which hangs on them. The
    // pinned matrix was brittle -- it stopped being a map at all when CMNONCOPRIME=1 changed
    // the presentation -- while the construction survives re-presentation and still PROVES
    // the isomorphism (it exhibits a map and certifies it with IsIsomorphism).
    // ⚠ model_drift_ok: this test pins a NON-ZERO base_label, and AllEquationsAboveCovers gates
    // EquationsByRebase on `base_label eq 0` (EquationsCovers.m:1061). So it cannot reproduce the
    // cover keys that the rebase FILLED on a default run -- [1,2], [1,5] and [1,26], which
    // data/models/models_10_13.m records as "previously EMPTY ... now filled, unlocked by
    // EquationsByRebase". MISSING keys only: a key this test DOES produce must still be the
    // committed curve, and model_drift_ok does not silence that.
    // ⚠ THE CONTROL GROUP is what makes this a diagnosis rather than an excuse: 14_3, 21_2 and
    // 6_17 also pin a base_label and all three PASS -- 14_3's empties were fixed by the coprime
    // filter flip, not the rebase, and the other two never had any. The gate costs exactly the
    // rebase-filled keys and nothing else.
    test_AllEquationsAboveCoversSingleCurve(10, 13, cover_data, ws_data, curves : model_drift_ok := true, base_label := 4069);
    return;
end procedure;

test_10_13();