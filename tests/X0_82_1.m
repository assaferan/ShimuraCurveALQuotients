import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_82_1()
    P3<x,y,z,s> := WeightedProjectiveSpace(Rationals(), [1,2,1,1]);
     // D = 82
    cover_data := AssociativeArray();
    // cover_data[{1}] is Guo-Yang's PUBLISHED curve (Table A.1: y^2 = 4s^4+4s^3+s^2-2s+1,
    // x^2 = -19s^2+18s-11), so this key is a genuine oracle.  The pipeline presents the same curve
    // over a different base coordinate; the isomorphism is now CONSTRUCTED, not pinned.
    //
    // ⚠⚠ THE PREVIOUS NOTE HERE WAS WRONG, AND IT WAS MINE (2026-09-27, same day).  It said the
    // hand-pinned matrix was LOAD-BEARING, that `construct_crv_isomorphism` "cannot match" the two,
    // and that the matrix would have to be re-derived by hand if the pipeline ever re-presented
    // this curve.  The real cause was a BUG IN THE HELPER, not a fact about these curves:
    // `construct_crv_isomorphism` tested whether Evaluate(f_expected, mu)/f_ours is a constant
    // square, but under t -> num/den the model y^2 = f(t) becomes f(mu)*den^deg.  Those agree only
    // when den is CONSTANT -- an AFFINE base change -- so every genuine Mobius was silently
    // declined.  Fixed in tests/_crviso.m; the affine path is unchanged.
    // ⇒ Here the correct base map is t -> (2-2t)/(4-3t), a true Mobius, which is exactly why this
    // base tripped the bug while 10_19 and 21_2 (both affine) did not.  With the factor restored
    // the ratios are 1024 and 1, whose square roots 32 and 1 are EXACTLY the y- and x-scales of the
    // matrix that used to be pinned here -- so the helper re-derives the hand-written bridge.
    // The constructed map is (x,y,s,z) -> (x/2, 8y, 3s/2-2z, s-z), projectively the pinned one.
    //
    // ⚠ AND THE SWEEP VERDICT RECORDED HERE WAS A FALSE NEGATIVE FROM THE SAME BUG.  It said
    // `magma -b Dd:=82 Nn:=1 tests/_basesweep.m` is NEGATIVE, the only candidate base being 10898
    // (whose pair is the committed model).  The sweep calls construct_crv_isomorphism, so its
    // "FALSE" carried no information about the curves.  ⇒ Treat every recorded
    // "construct_crv_isomorphism: FALSE" from before this fix as UNVERIFIED, not as a result.
    //
    // ⚠ NOT the 10_19 situation, and worth stating so nobody tries that remedy here: D*N = 82 has
    // only FOUR Atkin-Lehner involutions, so there is exactly ONE V_4 and both presentations must
    // use it.  Branch genera are distinct (w_2 -> 2, w_41 -> 0, w_82 -> 1), so the branch-to-
    // involution assignment is forced and identical on both sides; both conics are pointless and
    // ramified at {2}; the y-branch Jacobians are isomorphic over Q with conductor 82.  There is
    // nothing to re-present.
    cover_data[{1}] := <Curve(P3, [y^2 - 4*s^4- 4*s^3*z-s^2*z^2+2*s*z^3-z^4, x^2 +19*s^2-18*s*z+11*z^2]), Matrix([[1,0,0,0],[0,32,0,0],[0,0,-3,-2],[0,0,4,2]])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2] := DiagonalMatrix([-1,-1,1,1]);
    ws_data[{1}][41] := DiagonalMatrix([1,-1,1,1]);

    return cover_data, ws_data;
end function;

procedure test_82_1()
    cover_data, ws_data := load_covers_and_ws_data_82_1();
    curves := GetHyperellipticCandidates();
    // ⚠ NO manual_isomorphism since 2026-09-27: the construct branch now derives the map itself
    // (0.03 s), so the matrix in cover_data above is an unused placeholder, kept only because it
    // records the value the helper independently reproduces.  This removes the brittleness the old
    // note warned about -- a re-presentation no longer breaks a hardcoded map, because there is none.
    test_AllEquationsAboveCoversSingleCurve(82, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_82_1();