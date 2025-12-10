import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_15_1()
    _<s> := PolynomialRing(Rationals());

    // verifying [Guo-Yang, Example 32, p. 22-24]
    // D = 15
    cover_data := AssociativeArray();
    cover_data[{1}] := <HyperellipticCurve(-1/3*(s^2+243)*(s^2+3)), DiagonalMatrix([9, 4*27, 1])>;
    cover_data[{1,3}] := <HyperellipticCurve(-1/3*(s+243)*(s+3)), DiagonalMatrix([-81/2, 4*27, 1])>;
    cover_data[{1,15}] := <HyperellipticCurve(s), DiagonalMatrix([-1/2, 1, 1]) >;

    ws_data := AssociativeArray();

    return cover_data, ws_data;
end function;

procedure test_15_1()
    cover_data, ws_data := load_covers_and_ws_data_15_1();
    curves := GetHyperellipticCandidates();
    // Not yet working
    test_AllEquationsAboveCoversSingleCurve(15, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_15_1();