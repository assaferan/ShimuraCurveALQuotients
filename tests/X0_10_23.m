import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_10_23()
    _<x> := PolynomialRing(Rationals());

    // verifying [Guo-Yang, Table A.2, p. 37]
    // D = 10, N = 23
    cover_data := AssociativeArray();
    cover_data[{1}] := <-43*x^20+318*x^19-1071*x^18+3014*x^17-10540*x^16+28266*x^15-72217*x^14+81478*x^13-62765*x^12-68732*x^11+18840*x^10+68732*x^9\
-62765*x^8-81478*x^7-72217*x^6-28266*x^5-10540*x^4-3014*x^3-1071*x^2-318*x-43, Matrix([[0,0,1/2],[0,1/1024,0],[1/2,0,1]])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2] := Matrix([[2,0,1],[0,-5^5,0],[1,0,-2]]);
    ws_data[{1}][5] := Matrix([[0,0,1],[0,-1,0],[-1,0,0]]);
    ws_data[{1}][230] := DiagonalMatrix([1, -1, 1]);

    return cover_data, ws_data;
end function;

procedure test_10_23()
    cover_data, ws_data := load_covers_and_ws_data_10_23();
    curves := GetHyperellipticCandidates();
    
    test_AllEquationsAboveCoversSingleCurve(10, 23, cover_data, ws_data, curves);
    return;
end procedure;

test_10_23();