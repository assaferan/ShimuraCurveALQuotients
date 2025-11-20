import "BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_62_1()
    _<s> := PolynomialRing(Rationals());

    // Verifying [GY, Table A.1, p. 34]
    // D = 62
    cover_data := AssociativeArray();
    cover_data[{1}] := <-64*s^8-99*s^6-90*s^4-43*s^2-8, DiagonalMatrix([1,4,1])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2] := DiagonalMatrix([-1,1,1]);
    ws_data[{1}][62] := DiagonalMatrix([1,-1,1]);

    return cover_data, ws_data;
end function;

procedure test_62_1()
    cover_data, ws_data := load_covers_and_ws_data_62_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(62, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_62_1();