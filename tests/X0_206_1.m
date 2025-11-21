import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_206_1()
    _<s> := PolynomialRing(Rationals());

    // D = 206
    cover_data := AssociativeArray();
    cover_data[{1}] := <-8*s^20+13*s^18+42*s^16+331*s^14+220*s^12-733*s^10-6646*s^8-19883*s^6-28840*s^4-18224*s^2-4096, Matrix([[0,0,-1],[0,1,0],[-1,0,0]])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2] := DiagonalMatrix([-1,1,1]);
    ws_data[{1}][206] := DiagonalMatrix([1,-1,1]);

    return cover_data, ws_data;
end function;

procedure test_206_1()
    cover_data, ws_data := load_covers_and_ws_data_206_1();
    curves := GetHyperellipticCandidates();
    // Not yet working
    // test_AllEquationsAboveCoversSingleCurve(206, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_206_1();