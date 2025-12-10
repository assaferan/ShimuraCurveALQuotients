import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_134_1()
    _<s> := PolynomialRing(Rationals());

    // Verifying [GY, Table A.1, p. 35]
    // D = 134
    cover_data := AssociativeArray();
    cover_data[{1}] := <HyperellipticCurve(-16*s^14-347*s^12-2518*s^10-13341*s^8-91876*s^6+32859*s^4-2518*s^2-67), DiagonalMatrix([-1,-8192,1])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2] := DiagonalMatrix([-1,-1,1]);
    ws_data[{1}][134] := DiagonalMatrix([1,-1,1]);

    return cover_data, ws_data;
end function;

procedure test_134_1()
    cover_data, ws_data := load_covers_and_ws_data_134_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(134, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_134_1();