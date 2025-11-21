import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_26_1()
    _<s> := PolynomialRing(Rationals());

    // verifying [Guo-Yang, Example 33, p. 24-25]
    cover_data := AssociativeArray();
    cover_data[{1,13}] := <-2*s^3+19*s^2-24*s-169, DiagonalMatrix([1,1,1])>;
    cover_data[{1,26}] := <s, DiagonalMatrix([1,1,1])>;
    cover_data[{1}] := <-2*s^6+19*s^4-24*s^2-169, DiagonalMatrix([1,1,1])>;

    ws_data:= AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2] := DiagonalMatrix([-1,-1,1]);
    ws_data[{1}][26] := DiagonalMatrix([1,-1,1]);

    return cover_data, ws_data;
end function;

procedure test_26_1()
    cover_data, ws_data := load_covers_and_ws_data_26_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(26, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_26_1();