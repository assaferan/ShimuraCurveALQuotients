import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_74_1()
    _<s> := PolynomialRing(Rationals());

     // D = 74
    cover_data := AssociativeArray();
    cover_data[{1}] := <HyperellipticCurve(-2*s^10+47*s^8-328*s^6+946*s^4-4158*s^2-1369), DiagonalMatrix([1,256,1])>;

    ws_data:= AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2] := DiagonalMatrix([-1,-1,1]);
    ws_data[{1}][74] := DiagonalMatrix([1,-1,1]);

    return cover_data, ws_data;
end function;

procedure test_74_1()
    cover_data, ws_data := load_covers_and_ws_data_74_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(74, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_74_1();