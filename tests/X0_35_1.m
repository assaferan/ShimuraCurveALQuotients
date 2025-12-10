import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_35_1()
    _<x> := PolynomialRing(Rationals());

    // verifying [Guo-Yang, Table A.1, p. 33]
    // D = 35
    cover_data := AssociativeArray();
    cover_data[{1}] := <-(x^2+7)*(7*x^6+51*x^4+197*x^2+1), DiagonalMatrix([-1/5, 16, 1])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][5] := DiagonalMatrix([-1,-1,1]);
    ws_data[{1}][35] := DiagonalMatrix([1,-1,1]);

    return cover_data, ws_data;
end function;

procedure test_35_1()
    cover_data, ws_data := load_covers_and_ws_data_35_1();
    curves := GetHyperellipticCandidates();
   
    test_AllEquationsAboveCoversSingleCurve(35, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_35_1();