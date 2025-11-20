import "BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_94_1()
    _<s> := PolynomialRing(Rationals());

    // D = 94
    cover_data := AssociativeArray();
    cover_data[{1}] := <-8*s^8+69*s^6-234*s^4+381*s^2-256, DiagonalMatrix([-4,800,1])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2] := DiagonalMatrix([-1,1,1]);
    ws_data[{1}][94] := DiagonalMatrix([1,-1,1]);

    return cover_data, ws_data;
end function;

procedure test_94_1()
    cover_data, ws_data := load_covers_and_ws_data_94_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(94, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_94_1();