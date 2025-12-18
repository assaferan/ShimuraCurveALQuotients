
procedure test_AllEquationsAboveCoversSingleCurve(D, N, cover_data, ws_data, curves : algebra_map := false, base_label := 0, manual_isomorphism := false)
    // no longer needed as we now have a test for each curve
    // printf "testing equations of covers of X0*(%o;%o)...", D, N;
    assert exists(Xstar){X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)};
    covers, ws := AllEquationsAboveCovers(Xstar, curves : base_label := base_label);
    for label in Keys(covers) do
        X := curves[label];
        is_def, datum := IsDefined(cover_data, X`W);
        if not is_def then continue; end if;
        C_ex, scales := Explode(datum);
        P<[x]> := AmbientSpace(C_ex);
        for base in Keys(covers[label]) do
            C := covers[label][base];
            if manual_isomorphism then
                if algebra_map then
                    phi := scales;
                else
                    phi := map<C -> C_ex | Eltseq(Vector(x)*ChangeRing(scales, Universe(x)))>;
                end if;
                is_isom := IsIsomorphism(phi);
            else
                is_isom, phi := IsIsomorphic(C, C_ex);
            end if;
            assert is_isom;
            ws_def, ws_ex := IsDefined(ws_data, X`W);
            if not ws_def then continue; end if;
            for Q in Keys(ws_ex) do
                w_alg := AlgebraMap(phi)*AlgebraMap(ws[label][base][Q])*AlgebraMap(phi^(-1));
                phi1 := map< C_ex -> C_ex | [w_alg(x[j]) : j in [1..#x]]>;
                phi2 := map< C_ex -> C_ex | Eltseq(Vector(x)*ChangeRing(ws_ex[Q], Universe(x)))>;
                assert phi1 eq phi2;
            end for;
        end for;
    end for;
    // no longer needed as we now have a test for each curve
    // printf "Done\n";
    return;
end procedure;

function load_covers_and_ws_data()
    _<s> := PolynomialRing(Rationals());
    
    cover_data := AssociativeArray();
    ws_data := AssociativeArray();

    // This is back to not working
    // verifying [Guo-Yang, Example 32, p. 22-24]
    /*
    cover_data[<15,1>] := AssociativeArray();
    cover_data[<15,1>][{1,3}] := <-1/3*(s+243)*(s+3), DiagonalMatrix([-243, 4*27, 1])>;
    cover_data[<15,1>][{1,15}] := <s, DiagonalMatrix([-3, 1, 1]) >;
    cover_data[<15,1>][{1}] := <-1/3*(s^2+243)*(s^2+3), DiagonalMatrix([9, 4*27, 1])>;
    */

    // verifying [Guo-Yang, Example 33, p. 24-25]
    cover_data[<26,1>] := AssociativeArray();
    cover_data[<26,1>][{1,13}] := <-2*s^3+19*s^2-24*s-169, DiagonalMatrix([1,1,1])>;
    cover_data[<26,1>][{1,26}] := <s, DiagonalMatrix([1,1,1])>;
    cover_data[<26,1>][{1}] := <-2*s^6+19*s^4-24*s^2-169, DiagonalMatrix([1,1,1])>;
    ws_data[<26,1>] := AssociativeArray();
    ws_data[<26,1>][{1}] := AssociativeArray();
    ws_data[<26,1>][{1}][2] := DiagonalMatrix([-1,-1,1]);
    ws_data[<26,1>][{1}][26] := DiagonalMatrix([1,-1,1]);

    // Verifying [GY, Table A.1, p. 33]
    // D = 38
    cover_data[<38,1>] := AssociativeArray();
    cover_data[<38,1>][{1}] := <-16*s^6-59*s^4-82*s^2-19, DiagonalMatrix([1,16,1])>;
    ws_data[<38,1>] := AssociativeArray();
    ws_data[<38,1>][{1}] := AssociativeArray();
    ws_data[<38,1>][{1}][2] := DiagonalMatrix([-1,-1,1]);
    ws_data[<38,1>][{1}][38] := DiagonalMatrix([1,-1,1]);

    // D = 58
    cover_data[<58,1>] := AssociativeArray();
    cover_data[<58,1>][{1}] := <-2*s^6-78*s^4-862*s^2-1682, Matrix([[0,0,-1],[0,-32,0],[-1,0,0]])>;
    ws_data[<58,1>] := AssociativeArray();
    ws_data[<58,1>][{1}] := AssociativeArray();
    ws_data[<58,1>][{1}][2] := DiagonalMatrix([-1,-1,1]);
    ws_data[<58,1>][{1}][29] := DiagonalMatrix([1,-1,1]);

    // Verifying [GY, Table A.1, p. 34]
    // D = 62
    cover_data[<62,1>] := AssociativeArray();
    cover_data[<62,1>][{1}] := <-64*s^8-99*s^6-90*s^4-43*s^2-8, DiagonalMatrix([1,4,1])>;
    ws_data[<62,1>] := AssociativeArray();
    ws_data[<62,1>][{1}] := AssociativeArray();
    ws_data[<62,1>][{1}][2] := DiagonalMatrix([-1,1,1]);
    ws_data[<62,1>][{1}][62] := DiagonalMatrix([1,-1,1]);

    // D = 74
    cover_data[<74,1>] := AssociativeArray();
    cover_data[<74,1>][{1}] := <-2*s^10+47*s^8-328*s^6+946*s^4-4158*s^2-1369, DiagonalMatrix([1,256,1])>;
    ws_data[<74,1>] := AssociativeArray();
    ws_data[<74,1>][{1}] := AssociativeArray();
    ws_data[<74,1>][{1}][2] := DiagonalMatrix([-1,-1,1]);
    ws_data[<74,1>][{1}][74] := DiagonalMatrix([1,-1,1]);

    // D = 86
    cover_data[<86,1>] := AssociativeArray();
    cover_data[<86,1>][{1}] := <-16*s^10+245*s^8-756*s^6-1506*s^4-740*s^2-43, DiagonalMatrix([1,256,1])>;
    ws_data[<86,1>] := AssociativeArray();
    ws_data[<86,1>][{1}] := AssociativeArray();
    ws_data[<86,1>][{1}][2] := DiagonalMatrix([-1,-1,1]);
    ws_data[<86,1>][{1}][86] := DiagonalMatrix([1,-1,1]);

    // D = 94
    cover_data[<94,1>] := AssociativeArray();
    cover_data[<94,1>][{1}] := <-8*s^8+69*s^6-234*s^4+381*s^2-256, DiagonalMatrix([-4,800,1])>;
    ws_data[<94,1>] := AssociativeArray();
    ws_data[<94,1>][{1}] := AssociativeArray();
    ws_data[<94,1>][{1}][2] := DiagonalMatrix([-1,1,1]);
    ws_data[<94,1>][{1}][94] := DiagonalMatrix([1,-1,1]);

    // Verifying [GY, Table A.1, p. 35]
    // D = 134
    cover_data[<134,1>] := AssociativeArray();
    cover_data[<134,1>][{1}] := <-16*s^14-347*s^12-2518*s^10-13341*s^8-91876*s^6+32859*s^4-2518*s^2-67, DiagonalMatrix([-1,-8192,1])>;
    ws_data[<134,1>] := AssociativeArray();
    ws_data[<134,1>][{1}] := AssociativeArray();
    ws_data[<134,1>][{1}][2] := DiagonalMatrix([-1,-1,1]);
    ws_data[<134,1>][{1}][134] := DiagonalMatrix([1,-1,1]);

    // verifying [GY, Example 35, p. 27]
    // D = 146
    cover_data[<146,1>] := AssociativeArray();
    cover_data[<146,1>][{1}] := <-11*s^16+82*s^15-221*s^14+214*s^13+133*s^12-360*s^11-170*s^10+676*s^9
                                 -150*s^8-676*s^7-170*s^6+360*s^5+133*s^4-214*s^3-221*s^2-82*s-11, Matrix([[-1,0,-1],[0,128,0],[-1,0,0]])>;
    cover_data[<146,1>][{1,73}] := <-11*s^8+82*s^7-309*s^6+788*s^5-1413*s^4+1858*s^3-1803*s^2+1240*s-688, Matrix([[-1,0,0],[0,128,0],[1,0,1]])>;
    cover_data[<146,1>][{1,146}] := <s^2 + 4, Matrix([[1,0,0],[0,1,0],[-1,0,1]])>;
    ws_data[<146,1>] := AssociativeArray();
    ws_data[<146,1>][{1}] := AssociativeArray();
    ws_data[<146,1>][{1}][146] := DiagonalMatrix([1,-1,1]);
    ws_data[<146,1>][{1}][73] := Matrix([[0,0,1],[0,1,0],[-1,0,0]]);

    // D = 194
    cover_data[<194,1>] := AssociativeArray();
    cover_data[<194,1>][{1}] := <-19*s^20-92*s^19-286*s^18-592*s^17-921*s^16-1016*s^15
                                 -872*s^14+460*s^13+1545*s^12+1752*s^11+34*s^10-1752*s^9
                                 +1545*s^8-460*s^7-872*s^6+1016*s^5-921*s^4+592*s^3-286*s^2+92*s-19, DiagonalMatrix([-1,1,1])>;
    ws_data[<194,1>] := AssociativeArray();
    ws_data[<194,1>][{1}] := AssociativeArray();
    ws_data[<194,1>][{1}][97] := Matrix([[0,0,1],[0,-1,0],[-1,0,0]]);
    ws_data[<194,1>][{1}][194] := DiagonalMatrix([1,-1,1]);

    // D = 206 - Not working yet!!!
    /*
    cover_data[<206,1>] := AssociativeArray();
    cover_data[<206,1>][{1}] := <-8*s^20+13*s^18+42*s^16+331*s^14+220*s^12-733*s^10-6646*s^8-19883*s^6-28840*s^4-18224*s^2-4096, 
                                  Matrix([[0,0,-1],[0,1,0],[-1,0,0]])>;
    ws_data[<206,1>] := AssociativeArray();
    ws_data[<206,1>][{1}] := AssociativeArray();
    ws_data[<206,1>][{1}][2] := DiagonalMatrix([-1,1,1]);
    ws_data[<206,1>][{1}][206] := DiagonalMatrix([1,-1,1]);
    */
    return cover_data, ws_data;
end function;

procedure test_AllEquationsAboveCovers()
    // we split it into smaller tests to allow parallelization
    /*
    cover_data, ws_data := load_covers_and_ws_data();
    curves := GetHyperellipticCandidates();
    for DN in Sort([x : x in Keys(cover_data)]) do
        D,N := Explode(DN);
        test_AllEquationsAboveCoversSingleCurve(D, N, cover_data[DN], ws_data[DN], curves);
    end for;
    */
    return;
end procedure;

// we split it into smaller tests to allow parallelization
test_AllEquationsAboveCovers();