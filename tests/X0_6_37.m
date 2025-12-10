import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_6_37()
    _<x> := PolynomialRing(Rationals());

    // verifying [Guo-Yang, Table A.2, p. 36]
    // D = 6, N = 37
    cover_data := AssociativeArray();
    cover_data[{1}] := <-4096*x^12-18480*x^10-40200*x^8-51595*x^6-40200*x^4-18480*x^2-4096, Matrix([[1,0,0], [0,4,0], [2,0,1]])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2] := DiagonalMatrix([-1,1,1]);
    ws_data[{1}][3] := Matrix([[0,0,1],[0,1,0],[1,0,0]]);
    ws_data[{1}][222] := DiagonalMatrix([1,-1,1]);

    return cover_data, ws_data;
end function;

procedure test_6_37()
    cover_data, ws_data := load_covers_and_ws_data_6_37();
    curves := GetHyperellipticCandidates();
    
    test_AllEquationsAboveCoversSingleCurve(6, 37, cover_data, ws_data, curves);
    return;
end procedure;

test_6_37();
