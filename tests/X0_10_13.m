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
    test_AllEquationsAboveCoversSingleCurve(10, 13, cover_data, ws_data, curves : base_label := 4069, manual_isomorphism);
    return;
end procedure;

test_10_13();