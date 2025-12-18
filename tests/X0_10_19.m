import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_10_19()
    P3<x,y,z,s> := WeightedProjectiveSpace(Rationals(), [1,3,1,1]);
     // D = 82
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
    test_AllEquationsAboveCoversSingleCurve(10, 19, cover_data, ws_data, curves);
    return;
end procedure;

test_10_19();