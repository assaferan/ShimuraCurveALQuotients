import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_6_31()
    _<x> := PolynomialRing(Rationals());

    // verifying [Guo-Yang, Table A.2, p. 36]
    // D = 6, N = 31
    cover_data := AssociativeArray();
    cover_data[{1}] := <-243*x^12+11882*x^10-177701*x^8+803948*x^6-1599309*x^4+962442*x^2-177147, Matrix([[-2,0,1], [0,65536,0], [-1,0,1]])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2] := Matrix([[0,0,1],[0,-27,0],[3,0,0]]);
    ws_data[{1}][3] := DiagonalMatrix([-1, 1, 1]);
    ws_data[{1}][186] := DiagonalMatrix([1, -1, 1]);

    return cover_data;
end function;

procedure test_6_31()
    cover_data := load_covers_and_ws_data_6_31();
    curves := GetHyperellipticCandidates();
    
    test_AllEquationsAboveCoversSingleCurve(6, 31, cover_data, ws_data, curves);
    return;
end procedure;

test_6_31();