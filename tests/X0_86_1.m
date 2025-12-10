import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_86_1()
    _<s> := PolynomialRing(Rationals());

     // D = 86
    cover_data := AssociativeArray();
    cover_data[{1}] := <HyperellipticCurve(-16*s^10+245*s^8-756*s^6-1506*s^4-740*s^2-43), DiagonalMatrix([1,256,1])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2] := DiagonalMatrix([-1,-1,1]);
    ws_data[{1}][86] := DiagonalMatrix([1,-1,1]);

    return cover_data, ws_data;
end function;

procedure test_86_1()
    cover_data, ws_data := load_covers_and_ws_data_86_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(86, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_86_1();