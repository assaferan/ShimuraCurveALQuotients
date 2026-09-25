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
    {Return the trace of V W_Q on the subspace of X fixed by Ws and V*W_Q}
    sum := 0;
    for w in Ws do
        neww := AtkinLehnerMul(Q, w, N*D);
        trformula,p := CanApplyTraceFormula(vname, neww);
        sgn := (-1)^#PrimeDivisors(GCD(neww,D));
        vprintf ShimuraQuotients, 2: "sgn = %o\n", sgn;
        if trformula then
            subtrace := TraceFormulaGamma0VWDNew(p, neww, D, N, 2);
            /*
            import "tests/trace_formula.m" : checkTraceVWDNew;  
            checkTraceVWDNew(p, neww, D, N, 2);
            */
            sum +:= sgn*subtrace; 
            vprintf ShimuraQuotients, 2: "trace for V = %o, Q = %o, w = %o, is %o\n", vname, Q, w, subtrace;
        else
            wmat := al_matrix(neww, D*N);
            g := V*wmat;
            VV, QQ := ModularInvolution(vname, D*N);
            QQ := AtkinLehnerMul(QQ, neww, D*N);
            assert VV eq V;
            g_subspaces := [ModularInvolution(vname, d*N)*al_matrix(GCD(neww, d*N), d*N) : d in Divisors(D)];
            subtrace := TraceFormulaGamma0gDNew(g_subspaces, QQ, D, N, 2);
            /*
            import "tests/trace_formula.m" : checkTracegDNew;  
            checkTracegDNew(Eltseq(g), g_subspaces, QQ, D, N, 2);
            */
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

// ---------------------------------------------------------------------------------------
// Which non-AL modular involutions v*W_o are tested on Y = X_0(D,N)/W.
//
// Both checks below compute g' = genus(Y/<v W_o>) and conclude from it (g' = 0, or g = 3 and
// g' = 2, or the Riemann-Hurwitz count fix = 2g - 4g' + 2).  That is valid only if v W_o is an
// INVOLUTION of Y, i.e. (v W_o)^2 lies in Q^x Gamma_0(DN) W and v W_o normalizes Q^x Gamma_0(DN) W.
// Relations used ([FH] Lemma 1, p.111, together with direct matrix computation in the
// normalizer of Gamma_0(DN) modulo Q^x Gamma_0(DN); here 2^a || N and
// eps(d) = 1 iff the 3-free part of d is 2 mod 3, which is multiplicative in d):
//   * S2 commutes with W_{p^v}, p odd, but S2 W_{2^a} has order 3 (a = 2) or 4 (a >= 3);
//     so S2 descends to Y iff every w in W is odd, and S2 W_o is an involution iff o is odd.
//   * V2 (8 | N) commutes with every W_Q.
//   * V3 (9 || N): V3 W_d V3^-1 = W_9^eps(d) W_d, so V3 descends iff 9 in W or eps(w) = 0 for
//     all w in W, and (V3 W_o)^2 = W_9^eps(o).
//   * S2 V3 = V3 S2.  V3 V2 V3^-1 = W_9^eps(2^a) V2 (conjugate V3 W_{2^a} V3^-1 by S2), so
//     (V2 V3 W_o)^2 = W_9^eps(2^a o): V2 behaves like W_{2^a} for this purpose.
//   * S2 V2 is never needed: S2 V2 W_o is an involution only when 2^a | o, and then
//     S2 V2 W_o = W_{2^a} (S2 W_{o/2^a}) W_{2^a}^-1, which is conjugate on Y (by the automorphism
//     w_{2^a}) to S2 W_{o/2^a} and so gives an isomorphic quotient, already tested.
// Unordered pairs are taken once: S2 V3 = V3 S2, and V3 V2 W_o = V2 V3 W_{9^e o} is covered as
// o runs over the Atkin-Lehner divisors.

function eps3(d)
    // 1 iff the 3-free part of d is 2 mod 3 ([FH] Lemma 1); multiplicative in d.
    return ((d div 3^Valuation(d, 3)) mod 3 eq 2) select 1 else 0;
end function;

intrinsic ModularNonALInvolutionCandidates(D::RngIntElt, N::RngIntElt, W::SetEnum)
    -> SeqEnum, SeqEnum, SeqEnum, SeqEnum
{The non-AL modular involutions v*W_o tested on X_0(D,N)/W by CheckModularNonALInvolutionTrace
and CheckModularNonALInvolutionModSym.  Returns the base names V_names (subset of S2, V2, V3
that descend to the quotient), the index sequences into V_names of the v that are tested
(singletons, then the unordered pairs S2 V3 and V2 V3), their names, and for each v the set
of Atkin-Lehner divisors o (o notin W, or o = 1) for which v*W_o is an involution of the
quotient.  See [FH] Lemma 1.}
    DN := D*N;
    V_names := [];
    if (N mod 4 eq 0) and &and[IsOdd(w) : w in W] then
        Append(~V_names, "S2");
    end if;
    if (N mod 8 eq 0) then
        Append(~V_names, "V2");
    end if;
    if (Valuation(N, 3) eq 2) then
        not_commute := false;
        if (9 notin W) then
            // not_commute := exists(w){w : w in W | exists(p){p : p in PrimeDivisors(w) | (p^Valuation(w,p) mod 3) eq 2}};
            not_commute := exists(w){w : w in W | eps3(w) eq 1};
        end if;
        if not not_commute then
            Append(~V_names, "V3");
        end if;
    end if;
    idx_sets := [[i] : i in [1..#V_names]];
    // S2 V2 is omitted: it never gives an involution not already tested (see above).
    idx_sets cat:= [[i, j] : j in [i+1..#V_names], i in [1..#V_names]
                            | {V_names[i], V_names[j]} ne {"S2", "V2"}];
    all_names := [&cat[(k eq 1 select "" else " ") cat V_names[I[k]] : k in [1..#I]] : I in idx_sets];
    als := [Q : Q in Divisors(DN) | GCD(Q, DN div Q) eq 1];
    ws := W diff {1};
    other_ws := {w : w in als | w notin ws};
    good_ws := [];
    for nm in all_names do
        parts := Split(nm, " ");
        bad_ws := {};
        if "S2" in parts then
            bad_ws join:= {w : w in other_ws | IsEven(w)};
        end if;
        if ("V3" in parts) and (9 notin W) then
            // (v W_o)^2 = W_9^eps(o * 2^a) if V2 is a factor of v, W_9^eps(o) otherwise.
            twist := ("V2" in parts) select 2^Valuation(N, 2) else 1;
            bad_ws join:= {w : w in other_ws | eps3(twist*w) eq 1};
        end if;
        Append(~good_ws, other_ws diff bad_ws);
    end for;
    return V_names, idx_sets, all_names, good_ws;
end intrinsic;

intrinsic IsModularInvolutionOnQuotient(g::AlgMatElt, W::SetEnum, L::RngIntElt) -> BoolElt
{For g in the normalizer of Gamma_0(L) (an integer or rational 2x2 matrix) and W a group of
Atkin-Lehner divisors of L, decide whether g induces an involution of X_0(L)/W: g is not in
Q^x Gamma_0(L) W, g^2 is, and g normalizes Q^x Gamma_0(L) W.}
    M2Q := MatrixAlgebra(Rationals(), 2);
    function in_QGamma(M)
        d := Determinant(M);
        if d le 0 then return false; end if;
        ok, c := IsSquare(d);
        if not ok then return false; end if;
        M1 := M / c;
        if not &and[IsIntegral(x) : x in Eltseq(M1)] then return false; end if;
        return (Integers()!M1[2,1]) mod L eq 0;
    end function;
    Wm := [M2Q!al_matrix(w, L) : w in W];
    in_grp := func<M | exists{w : w in Wm | in_QGamma(M * w^-1)}>;
    g := M2Q!g;
    if in_grp(g) then return false; end if;
    if not in_grp(g*g) then return false; end if;
    return &and[in_grp(g * w * g^-1) : w in Wm];
end intrinsic;

intrinsic CheckModularNonALInvolutionTrace(X::ShimuraQuot) -> RngIntElt, MonStgElt, RngIntElt
{Returns 1 if any of the non-AL modular involutions is hyperelliptic, in which case also returns the hyperelliptic involution,
returns 0 if the curve is non-hyperelliptic, and the involution with too many fixed points.
Otherwise, returns -1.}
    assert X`g ne 0;
    // We want the D-new subspace

    // Only v*W_o that are involutions of X/W are tested; see ModularNonALInvolutionCandidates.
    V_names, idx_sets, all_names, good_ws := ModularNonALInvolutionCandidates(X`D, X`N, X`W);
    Vs := [];
    for vname in V_names do
        if vname eq "S2" then
            Append(~Vs, Matrix(Integers(),2,2,[2,1,0,2]));
        elif vname eq "V2" then
            V2 := get_V2(X`D*X`N);
            Append(~Vs, V2);
        elif vname eq "V3" then
            V3 := get_V3(X`D*X`N);
            Append(~Vs, V3);
        end if;
    end for;
    all_vs := [&*[Vs[i] : i in I] : I in idx_sets];
    for idx->V_SN in all_vs do
        for other_w in good_ws[idx] do
            // Guard: never draw a conclusion from an element that is not an involution of X/W.
            if not IsModularInvolutionOnQuotient(V_SN*al_matrix(other_w, X`D*X`N), X`W, X`D*X`N) then
                vprintf ShimuraQuotients, 1: "WARNING: %o W%o is not an involution on %o; skipping\n", all_names[idx], other_w, X;
                continue;
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
intrinsic NonALModSymMaxLevel() -> RngIntElt
{Maximum level D*N for which the modular non-AL involution check runs ModularSymbols(D*N).
Above this, CheckModularNonALInvolutionModSym returns "undetermined" instead, because the
modular-symbols computation at large highly-composite levels costs many hours and GB; those
few large-level curves are determined/pruned by the other stages. Tunable.}
    return 3000;
end intrinsic;

intrinsic CheckModularNonALInvolutionModSym(X::ShimuraQuot) -> RngIntElt, MonStgElt, RngIntElt
{Returns 1 if any of the non-AL modular involutions is hyperelliptic, in which case also returns the hyperelliptic involution,
returns 0 if the curve is non-hyperelliptic, and the involution with too many fixed points.
Otherwise, returns -1.}
    assert X`g ne 0;
    has_modularnonALinvolutions := false;
    if (X`N mod 4 eq 0) or (Valuation(X`N, 3) eq 2) then has_modularnonALinvolutions := true; end if;
    if not has_modularnonALinvolutions then
        vprintf ShimuraQuotients, 2: "The curve %o has no non-AL modular involutions\n", X;
        return -1, _, _;
    end if;
    // Skip curves whose level D*N is too large for an affordable ModularSymbols computation.
    // At D*N ~ 7000-8000 the modular-symbols space plus its new-subspace/Atkin-Lehner linear
    // algebra runs for many hours and many GB; such curves are left undetermined here and are
    // handled (and pruned) by the cheaper filtering stages instead.
    if X`D*X`N gt NonALModSymMaxLevel() then
        vprintf ShimuraQuotients, 2: "The curve %o exceeds the non-AL ModSym level cap\n", X;
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
    // Only v*W_o that are involutions of X/W are tested; see ModularNonALInvolutionCandidates.
    V_names, idx_sets, all_names, good_ws := ModularNonALInvolutionCandidates(X`D, X`N, X`W);
    Vs := [];
    Vmats := [];
    for vname in V_names do
        if vname eq "S2" then
            Append(~Vs, get_Vmu(2, X`N, SDN_new_basis, MDN, true));
            Append(~Vmats, Matrix(Integers(),2,2,[2,1,0,2]));
        elif vname eq "V2" then
            Append(~Vs, get_Vmu(2, X`N, SDN_new_basis, MDN, false));
            Append(~Vmats, get_V2(X`D*X`N));
        elif vname eq "V3" then
            Append(~Vs, get_Vmu(3, X`N, SDN_new_basis, MDN, false));
            Append(~Vmats, get_V3(X`D*X`N));
        end if;
    end for;
    all_vs := [&*[Vs[i] : i in I] : I in idx_sets];
    all_vmats := [&*[Vmats[i] : i in I] : I in idx_sets];
    ws := X`W diff {1};
    // The (chi-twisted) W-invariant subspace; the same for every candidate.
    W_fixed := VectorSpace(Rationals(), Nrows(SDN_new_basis));
    for w in ws do
        W_MDN := AtkinLehnerOperator(MDN, w);
        W_SN := Solution(SDN_new_basis, SDN_new_basis * W_MDN);
        al_sign := (-1)^#PrimeDivisors(GCD(w, X`D));
        W_fixed meet:= Kernel(Matrix(W_SN) - al_sign);
    end for;
    W_fixed_basis := BasisMatrix(W_fixed);
    for idx->V_SN in all_vs do
        for other_w in good_ws[idx] do
            // Guard: never draw a conclusion from an element that is not an involution of X/W.
            if not IsModularInvolutionOnQuotient(all_vmats[idx]*al_matrix(other_w, X`D*X`N), X`W, X`D*X`N) then
                vprintf ShimuraQuotients, 1: "WARNING: %o W%o is not an involution on %o; skipping\n", all_names[idx], other_w, X;
                continue;
            end if;
            other_W_MDN := AtkinLehnerOperator(MDN, other_w);
            other_W_SN := Solution(SDN_new_basis, SDN_new_basis * other_W_MDN);
            // Atkin-Lehner operators at primes dividing D act with opposite signs on the
            // Shimura and modular Jacobians, so the Shimura-side invariants are the
            // chi(m) = (-1)^#PD(gcd(m,D)) eigenspaces here.  chi is a *character* of the
            // AL group, so the sign of each W_w is chi(w) (not chi(w*other_w)), and the
            // V W_{other} involution is taken in its chi(other_w) eigenspace.
            chiQ := (-1)^#PrimeDivisors(GCD(other_w, X`D));
            op := Matrix(V_SN * other_W_SN);
            // Second guard, on the operator itself: it must preserve the W-invariant subspace
            // and square to the identity there.
            if (Nrows(W_fixed_basis) gt 0) and
               ((RowSpace(W_fixed_basis * op) notsubset W_fixed) or (W_fixed_basis * op * op ne W_fixed_basis)) then
                vprintf ShimuraQuotients, 1: "WARNING: %o W%o does not act as an involution on the W-invariants of %o; skipping\n", all_names[idx], other_w, X;
                continue;
            end if;
            fixed_subspace := Kernel(op - chiQ) meet W_fixed;
            d := Dimension(fixed_subspace);
            // print "d = ", d;
            assert IsEven(d);
            g := d div 2;
            //print "g = ", g;
            name := all_names[idx] cat " " cat Sprintf("W%o", other_w);
            if (g eq 0) or ((X`g eq 3) and (g eq 2)) then
                return 1, name, _;
            end if;
            fix := 2*X`g - 4*g + 2;
            // If X`g is even, fix = 2 mod 4, and so fix ne 2 is equivalent to fix gt 2
            if IsEven(X`g) and fix gt 2 then
                return 0, name, fix;
            // If X`g is odd fix = 0 mod 4, so fix gt 4 is equivalent to fix notin {0,4}
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