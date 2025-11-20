import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_38_1()
    _<s> := PolynomialRing(Rationals());

    // Verifying [GY, Table A.1, p. 33]
    // D = 38
    cover_data := AssociativeArray();
    cover_data[{1}] := <-16*s^6-59*s^4-82*s^2-19, DiagonalMatrix([1,16,1])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2] := DiagonalMatrix([-1,-1,1]);
    ws_data[{1}][38] := DiagonalMatrix([1,-1,1]);

    return cover_data, ws_data;
end function;

procedure test_38_1()
    cover_data, ws_data := load_covers_and_ws_data_38_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(38, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_38_1();