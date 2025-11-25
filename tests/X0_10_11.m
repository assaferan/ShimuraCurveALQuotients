import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_10_11()
    _<x> := PolynomialRing(Rationals());

    // verifying [Guo-Yang, Table A.2, p. 37]
    // D = 10, N = 11
    cover_data := AssociativeArray();
    cover_data[{1}] := <-8*x^12-35*x^10+30*x^8+277*x^6+120*x^4-560*x^2-512, Matrix([[-1,0,0],[0,1,0],[-1,0,1]])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][10] := Matrix([[0,0,1],[0,-8,0],[-2,0,0]]);
    ws_data[{1}][22] := Matrix([[0,0,1],[0,8,0],[2,0,0]]);
    ws_data[{1}][110] := DiagonalMatrix([1, -1, 1]);

    return cover_data, ws_data;
end function;

procedure test_10_11()
    cover_data, ws_data := load_covers_and_ws_data_10_11();
    curves := GetHyperellipticCandidates();
    
    test_AllEquationsAboveCoversSingleCurve(10, 11, cover_data, ws_data, curves);
    return;
end procedure;

test_10_11();