import !"Geometry/ModSym/operators.m" : ActionOnModularSymbolsBasis;
// Following [FH] Section 2

function get_Vmu(mu, N, SDN_new_basis, MDN, get_S)
    Smu_MDN := ActionOnModularSymbolsBasis([mu,1,0,mu], MDN);
    if (get_S) then
        Smu_SN := Solution(SDN_new_basis, SDN_new_basis * Smu_MDN);
        return Smu_SN;
    end if;
    W_MDN := AtkinLehnerOperator(MDN, mu^Valuation(N, mu));
    Vmu_MDN := Smu_MDN*W_MDN*Smu_MDN^(-1);
    Vmu_SN := Solution(SDN_new_basis, SDN_new_basis * Vmu_MDN);
    return Vmu_SN;
end function;

// Start with an implmentation based on modular symbols
function CheckModularNonALInvolutionModSym(X)
    assert X`g ne 0;
    MDN := ModularSymbols(X`D*X`N, 2, 0);
    SDN := CuspidalSubspace(MDN);
    // We want the D-new subspace
    ps := PrimeDivisors(X`D);
    SDN_new := SDN;
    for p in ps do
        SDN_new := NewSubspace(SDN_new, p);
    end for;
    SDN_new_basis := Matrix([Representation(v) : v in Basis(SDN_new)]);
    Vs := [];
    V_names := [];
    if (X`N mod 4 eq 0) and &and[IsOdd(w) : w in X`W] then
        Append(~Vs, get_Vmu(2, X`N, SDN_new_basis, MDN, true));
        Append(~V_names, "S2");
    end if;
    if (X`N mod 8 eq 0) then
        Append(~Vs, get_Vmu(2, X`N, SDN_new_basis, MDN, false));
        Append(~V_names, "V2");
    end if;
    if (Valuation(X`N, 3) eq 2) then
        not_commute := false;
        if (9 notin X`W) then
            not_commute := exists(w){w : w in X`W | exists(p){p : p in PrimeDivisors(w) | (p^Valuation(w,p) mod 3) eq 2}};
        end if;
        if not not_commute then
            Append(~Vs, get_Vmu(3, X`N, SDN_new_basis, MDN, false));
            Append(~V_names, "V3");
        end if;
    end if;
    all_vs := Vs cat [v1*v2 : v1, v2 in Vs | v1 ne v2];
    all_names := V_names cat [v1 cat " " cat v2 : v1, v2 in V_names | v1 ne v2];
    als := [Q : Q in Divisors(X`D*X`N) | GCD(Q, X`D*X`N div Q) eq 1]; 
    ws := X`W diff {1};
    other_ws := {w : w in als | w notin ws};
    for idx->V_SN in all_vs do
        bad_ws := {};
        if "S2" in all_names[idx] then
            bad_ws join:= {w : w in other_ws | IsEven(w)};
        end if;
        if ("V3" in all_names[idx]) and (9 notin X`W) then
            bad_ws join:= {w : w in other_ws | exists(p){p : p in PrimeDivisors(w) | (p^Valuation(w,p) mod 3) eq 2}};
        end if;
        for other_w in (other_ws diff bad_ws) do
            other_W_MDN := AtkinLehnerOperator(MDN, other_w);
            other_W_SN := Solution(SDN_new_basis, SDN_new_basis * other_W_MDN);
            fixed_subspace := Kernel(Matrix(V_SN * other_W_SN) - 1);
            for w in ws do
                W_MDN := AtkinLehnerOperator(MDN, w);
                W_SN := Solution(SDN_new_basis, SDN_new_basis * W_MDN);
                al_sign := (X`D mod w eq 0) select -1 else 1;
                kerw := Kernel(Matrix(W_SN) - al_sign);
                fixed_subspace meet:= kerw;
            end for;
            d := Dimension(fixed_subspace);
            // print "d = ", d;
            assert IsEven(d);
            g := d div 2;
            // print "g = ", g;
            if (g eq 0) then
                return true, all_names[idx] cat " " cat Sprintf("W%o", other_w);
            end if;
        end for;
    end for;    
    return false;
end function;