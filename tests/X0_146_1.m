import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_146_1()
    _<s> := PolynomialRing(Rationals());

    // Verifying [GY, Example 35, p. 27]
    // D = 146
    cover_data := AssociativeArray();
    cover_data[{1}] := <-11*s^16+82*s^15-221*s^14+214*s^13+133*s^12-360*s^11-170*s^10+676*s^9
                                 -150*s^8-676*s^7-170*s^6+360*s^5+133*s^4-214*s^3-221*s^2-82*s-11, Matrix([[-1,0,-1],[0,128,0],[-1,0,0]])>;
    cover_data[{1,73}] := <-11*s^8+82*s^7-309*s^6+788*s^5-1413*s^4+1858*s^3-1803*s^2+1240*s-688, Matrix([[-1,0,0],[0,128,0],[1,0,1]])>;
    cover_data[{1,146}] := <s^2 + 4, Matrix([[1,0,0],[0,1,0],[-1,0,1]])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][146] := DiagonalMatrix([1,-1,1]);
    ws_data[{1}][73] := Matrix([[0,0,1],[0,1,0],[-1,0,0]]);

    return cover_data, ws_data;
end function;

procedure test_146_1()
    cover_data, ws_data := load_covers_and_ws_data_146_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(146, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_146_1();