import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_82_1()
    P3<x,y,z,s> := WeightedProjectiveSpace(Rationals(), [1,2,1,1]);
     // D = 82
    cover_data := AssociativeArray();
    // ⚠⚠ THE PINNED MATRIX BELOW IS LOAD-BEARING, AND base_label CANNOT RESCUE IT (2026-09-27).
    // cover_data[{1}] is Guo-Yang's PUBLISHED curve (Table A.1: y^2 = 4s^4+4s^3+s^2-2s+1,
    // x^2 = -19s^2+18s-11), so this key is a genuine oracle -- but the pipeline presents the same
    // curve differently, and the second component is a hand-pinned coordinate change bridging the
    // two under manual_isomorphism.
    // `magma -b Dd:=82 Nn:=1 tests/_basesweep.m` was run and is NEGATIVE: the only candidate base
    // is 10898, whose pair is the committed model
    //     y^2 - 169/1024 s^4 + 183/256 s^3 z - 301/256 s^2 z^2 + 7/8 s z^3 - 1/4 z^4 ,
    //     x^2 + 67 s^2 - 164 s z + 108 z^2
    // and construct_crv_isomorphism cannot match it to Guo-Yang's.  So if the pipeline ever
    // re-presents this curve, the pinned matrix stops being a map and this test fails exactly as
    // tests/_offline/X0_10_19.m did -- and sweeping the base will NOT fix it; the matrix has to be
    // re-derived.  Recorded so that hour is not spent twice.
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
    test_AllEquationsAboveCoversSingleCurve(82, 1, cover_data, ws_data, curves : manual_isomorphism);
    return;
end procedure;

test_82_1();