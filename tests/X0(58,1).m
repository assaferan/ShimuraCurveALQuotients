import "BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_58_1()
    _<s> := PolynomialRing(Rationals());

    // D = 58
    cover_data := AssociativeArray();
    cover_data[{1}] := <-2*s^6-78*s^4-862*s^2-1682, Matrix([[0,0,-1],[0,-32,0],[-1,0,0]])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2] := DiagonalMatrix([-1,-1,1]);
    ws_data[{1}][29] := DiagonalMatrix([1,-1,1]);

    return cover_data, ws_data;
end function;

procedure test_58_1()
    cover_data, ws_data := load_covers_and_ws_data_58_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(58, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_58_1();