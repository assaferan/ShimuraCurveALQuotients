import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_194_1()
    _<s> := PolynomialRing(Rationals());

    // D = 194
    cover_data := AssociativeArray();
    cover_data[{1}] := <-19*s^20-92*s^19-286*s^18-592*s^17-921*s^16-1016*s^15
                                 -872*s^14+460*s^13+1545*s^12+1752*s^11+34*s^10-1752*s^9
                                 +1545*s^8-460*s^7-872*s^6+1016*s^5-921*s^4+592*s^3-286*s^2+92*s-19, DiagonalMatrix([-1,1,1])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][97] := Matrix([[0,0,1],[0,-1,0],[-1,0,0]]);
    ws_data[{1}][194] := DiagonalMatrix([1,-1,1]);

    return cover_data, ws_data;
end function;

procedure test_194_1()
    cover_data, ws_data := load_covers_and_ws_data_194_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(194, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_194_1();