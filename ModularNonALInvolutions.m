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


intrinsic CanApplyTraceFormula(vname, Q) -> BoolElt, RngIntElt
    {From the data of the name of the V and the Al operator, do we know a trace formula that is fast}
    can_apply_trace_formula := false;
    if vname eq "V2" then
        p := 2;
        if Valuation(Q,2) eq 0 then
            return true, p;
        end if;
    elif vname eq "V3" then
        p := 3;
        if (Valuation(Q, 3) eq 0) and (Q mod 3 eq 1) then
            return true, p;
        end if;
    end if;
        return false, _;
end intrinsic;

intrinsic TraceDNewQuotient(V::AlgMatElt, vname::MonStgElt, Q::RngIntElt, Ws::SetEnum, D::RngIntElt, N::RngIntElt) -> RngIntElt
    {Return the trace on the trace on the subspace of X fixed by Ws and V*W_Q}
    sum := 0;
    for w in Ws do
        neww := al_mul(Q, w, N*D);
        trformula,p := CanApplyTraceFormula(vname, neww);
        sgn := (-1)^#PrimeDivisors(GCD(neww,D));
        if trformula then
            subtrace := TraceFormulaGamma0VWDNew(p, neww, D, N, 2);
            sum +:= sgn*subtrace; 
            vprintf ShimuraQuotients, 2: "trace for V = %o, Q = %o, w = %o, is %o\n", vname, Q, w, subtrace;
        else
            wmat := al_matrix(neww, D*N);
            g := V*wmat;
            subtrace := TraceFormulaGamma0gDNew(Eltseq(g), D, N, 2);
            sum +:= sgn*subtrace; 
            vprintf ShimuraQuotients, 2: "trace for V = %o, Q = %o, w = %o, is %o, \n", vname, Q, w, subtrace;
            //print sgn;
        end if;
    end for;
    //print sum;
    sum *:= 1/#Ws;

    tr := TraceDNewALFixed(D,N,2, 1, Ws);

    return Integers()!(sum + tr) div 2 ;

end intrinsic;

intrinsic CheckModularNonALInvolutionTrace(X::ShimuraQuot) -> RngIntElt, MonStgElt, RngIntElt
{Returns 1 if any of the non-AL modular involutions is hyperelliptic, in which case also returns the hyperelliptic involution,
returns 0 if the curve is non-hyperelliptic, and the involution with too many fixed points.
Otherwise, returns -1.}
    assert X`g ne 0;
    // We want the D-new subspace

    Vs := [];
    V_names := [];
    if (X`N mod 4 eq 0) and &and[IsOdd(w) : w in X`W] then
        Append(~Vs, Matrix(Integers(),2,2,[2,1,0,2]));
        Append(~V_names, "S2");
    end if;
    if (X`N mod 8 eq 0) then
        V2 := get_V2(X`D*X`N);
        Append(~Vs, V2);
        Append(~V_names, "V2");
    end if;
    if (Valuation(X`N, 3) eq 2) then
        not_commute := false;
        if (9 notin X`W) then
            // not_commute := exists(w){w : w in X`W | exists(p){p : p in PrimeDivisors(w) | (p^Valuation(w,p) mod 3) eq 2}};
            not_commute := exists(w){w : w in X`W | (w div 3^Valuation(w, 3)) mod 3 eq 2};
        end if;
        if not not_commute then
            V3 := get_V3(X`D*X`N);
            Append(~Vs, V3);
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
            //two cases
            //we have a trace formula:
            can_apply_trace_formula, p := CanApplyTraceFormula(all_names[idx], other_w);
            if can_apply_trace_formula then
                trV := TraceFormulaGamma0VWDNew(p,other_w,X`D,X`N,2);
            else
                W := al_matrix(other_w, X`D*X`N);
                g := V_SN*W;
                trV:= TraceFormulaGamma0gDNew(Eltseq(g),X`D, X`N, 2); 
            end if;

            tr := TraceDNewQuotient(V_SN, all_names[idx], other_w,X`W,X`D, X`N);
            // print "g = ", g;
            name := all_names[idx] cat " " cat Sprintf("W%o", other_w);
            vprintf ShimuraQuotients, 2: "trace is %o\n", tr;
            if (tr eq 0) then
                return 1, name, _;
            end if;
            fix := 2*X`g - 4*tr + 2;
            if IsEven(X`g) and fix gt 2 then
                return 0, name, fix;
            elif IsOdd(X`g) and fix gt 4 then
                return 0, name, fix;
            end if;
        end for;
    end for;    
    return -1, _, _;
end intrinsic;


// Start with an implmentation based on modular symbols
intrinsic CheckModularNonALInvolutionModSym(X::ShimuraQuot) -> RngIntElt, MonStgElt, RngIntElt
{Returns 1 if any of the non-AL modular involutions is hyperelliptic, in which case also returns the hyperelliptic involution,
returns 0 if the curve is non-hyperelliptic, and the involution with too many fixed points.
Otherwise, returns -1.}
    assert X`g ne 0;
    has_modularnonALinvolutions := false;
    if (X`N mod 4 eq 0) or (Valuation(X`N, 3) eq 2) then has_modularnonALinvolutions eq true; end if;
    if not has_modularnonALinvolutions then
        return -1, _, _;
    end if;
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
            //print "g = ", g;
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