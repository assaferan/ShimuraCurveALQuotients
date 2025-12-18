import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_82_1()
    P3<x,y,z,s> := WeightedProjectiveSpace(Rationals(), [1,2,1,1]);
     // D = 82
    cover_data := AssociativeArray();
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
    test_AllEquationsAboveCoversSingleCurve(82, 1, cover_data, ws_data, curves : manual_isomoprhism);
    return;
end procedure;

test_82_1();