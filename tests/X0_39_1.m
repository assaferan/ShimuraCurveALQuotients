import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_39_1()
    _<x> := PolynomialRing(Rationals());

    // verifying [Guo-Yang, Table A.1, p. 33]
    // D = 39
    cover_data := AssociativeArray();
    cover_data[{1}] := <-(x^4-x^3-x^2+x+1)*(7*x^4-23*x^3+5*x^2+23*x+7), Matrix([[0,0,1],[0,3,0],[-1,0,0]])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][13] := Matrix([[0,0,1],[0,1,0],[-1,0,0]]);
    ws_data[{1}][39] := DiagonalMatrix([1,-1,1]);

    return cover_data, ws_data;
end function;

procedure test_39_1()
    cover_data, ws_data := load_covers_and_ws_data_39_1();
    curves := GetHyperellipticCandidates();
    
    test_AllEquationsAboveCoversSingleCurve(39, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_39_1();