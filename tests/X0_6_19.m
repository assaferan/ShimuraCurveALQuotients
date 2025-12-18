import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_6_19()
    _<s> := PolynomialRing(Rationals());

    // verifying [Guo-Yang, Table A.2, p. 36]
    // D = 6, N = 19
    cover_data := AssociativeArray();
    cover_data[{1}] := <HyperellipticCurve(-19*s^8+210*s^6-625*s^4+210*s^2-19), Matrix([[1,0,-1], [0,256,0], [1,0,1]])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2] := Matrix([[0,0,-1], [0,-1,0], [1,0,0]]);
    ws_data[{1}][3] := Matrix([[0,0,1], [0,1,0], [1,0,0]]);
    ws_data[{1}][114] := DiagonalMatrix([1, -1, 1]);

    return cover_data, ws_data;
end function;

procedure test_6_19()
    cover_data, ws_data := load_covers_and_ws_data_6_19();
    curves := GetHyperellipticCandidates();
    
    test_AllEquationsAboveCoversSingleCurve(6, 19, cover_data, ws_data, curves);
    return;
end procedure;

test_6_19();