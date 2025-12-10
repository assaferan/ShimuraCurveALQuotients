import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_6_17()
    P3<x,y,z,s> := WeightedProjectiveSpace(Rationals(), [1,2,1,1]);
     // D = 82
    cover_data := AssociativeArray();
    cover_data[{1}] := <Curve(P3, [y^2 - 17*s^4 + 10*s^2*z^2 - 9*z^4, x^2 + 3*s^2 + 16*z^2]), DiagonalMatrix([3,4,1,1])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2] := DiagonalMatrix([1,1,1,-1]);
    ws_data[{1}][3] := DiagonalMatrix([-1,-1,1,1]);
    ws_data[{1}][34] := DiagonalMatrix([1,-1,1,1]);

    return cover_data, ws_data;
end function;

procedure test_6_17()
    cover_data, ws_data := load_covers_and_ws_data_6_17();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(6, 17, cover_data, ws_data, curves : base_label := 1557);
    return;
end procedure;

test_6_17();