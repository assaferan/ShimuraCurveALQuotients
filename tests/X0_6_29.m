import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_6_29()
    _<x> := PolynomialRing(Rationals());

    // verifying [Guo-Yang, Table A.2, p. 36]
    // D = 6, N = 29
    cover_data := AssociativeArray();
    cover_data[{1}] := <HyperellipticCurve(-64*x^12+813*x^10-3066*x^8+4597*x^6-12264*x^4+13008*x^2-4096), Matrix([[0,0,1], [0,2592,0], [2,0,1]])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2] := DiagonalMatrix([-1,1,1]);
    ws_data[{1}][3] := Matrix([[0,0,1],[0,8,0],[-2,0,0]]);
    ws_data[{1}][174] := DiagonalMatrix([1, -1, 1]);

    return cover_data, ws_data;
end function;

procedure test_6_29()
    cover_data, ws_data := load_covers_and_ws_data_6_29();
    curves := GetHyperellipticCandidates();
    
    test_AllEquationsAboveCoversSingleCurve(6, 29, cover_data, ws_data, curves);
    return;
end procedure;

test_6_29();