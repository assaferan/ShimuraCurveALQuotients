import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_6_11()
    _<s> := PolynomialRing(Rationals());

    // verifying [Guo-Yang, Table A.2, p. 36]
    // D = 6, N = 11
    cover_data := AssociativeArray();
    cover_data[{1}] := <-19*s^8-166*s^7-439*s^6-166*s^5+612*s^4+166*s^3-439*s^2+166*s-19, DiagonalMatrix([1, 3, 1])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][3] := Matrix([[1,0,1],[0,-4,0],[1,0,-1]]);
    ws_data[{1}][2] := Matrix([[0,0,-1],[0,-1,0],[1,0,0]]);
    ws_data[{1}][66] := DiagonalMatrix([1, -1, 1]);

    return cover_data, ws_data;
end function;

procedure test_6_11()
    cover_data, ws_data := load_covers_and_ws_data_6_11();
    curves := GetHyperellipticCandidates();
    
    test_AllEquationsAboveCoversSingleCurve(6, 11, cover_data, ws_data, curves);
    return;
end procedure;

test_6_11();