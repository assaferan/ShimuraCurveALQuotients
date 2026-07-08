
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
