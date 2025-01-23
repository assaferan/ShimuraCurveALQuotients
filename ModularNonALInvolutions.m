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
intrinsic CheckModularNonALInvolutionModSym(X::ShimuraQuot) -> RngIntElt, MonStgElt, RngIntElt
{Returns 1 if any of the non-AL modular involutions is hyperelliptic, in which case also returns the hyperelliptic involution,
returns 0 if the curve is non-hyperelliptic, and the involution with too many fixed points.
Otherwise, returns -1.}
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
            // not_commute := exists(w){w : w in X`W | exists(p){p : p in PrimeDivisors(w) | (p^Valuation(w,p) mod 3) eq 2}};
            not_commute := exists(w){w : w in X`W | (w div 3^Valuation(w, 3)) mod 3 eq 2};
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
            name := all_names[idx] cat " " cat Sprintf("W%o", other_w);
            if (g eq 0) then
                return 1, name, _;
            end if;
            fix := 2*X`g - 4*g + 2;
            if IsEven(X`g) and fix gt 2 then
                return 0, name, fix;
            elif IsOdd(X`g) and fix gt 4 then
                return 0, name, fix;
            end if;
        end for;
    end for;    
    return -1, _, _;
end intrinsic;

intrinsic FilterByNonALInvolutions(~curves::SeqEnum[ShimuraQuot])
{Update curves in which the non-AL modular involutions are hyperelliptic.}
    for i->X in curves do
        if (i mod 500 eq 0) then
            vprintf ShimuraQuotients, 1: "i = %o/%o\n", i, #curves;
        end if; 
        if assigned X`IsSubhyp then continue; end if;
        is_hyp, inv, num := CheckModularNonALInvolutionModSym(X);
        if (is_hyp eq 1) then
            curves[i]`IsSubhyp := true;
            // we assume that if the genus is 0 or 1 this is already assigned
            curves[i]`IsHyp := true;
            curves[i]`TestInWhichProved := "ModularNonALInvolution " cat inv;
        end if;
        if (is_hyp eq 0) then
            curves[i]`IsSubhyp := false;
            // we assume that if the genus is 0 or 1 this is already assigned
            curves[i]`IsHyp := false;
            curves[i]`TestInWhichProved := "ModularNonALInvolution " cat inv cat " " cat Sprintf("%o", num);
        end if;
    end for;
    return;
end intrinsic;