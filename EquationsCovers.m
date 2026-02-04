
intrinsic RationalConstraintsOnEquations(schofer_table::SchoferTable, curves::SeqEnum[ShimuraQuot]) -> SeqEnum
{Impose linear constraints on equations given by rational CM points}

    keys_fs := schofer_table`Keys_fs;
    table := schofer_table`Values;
    ds := schofer_table`Discs;
    k_idxs := schofer_table`K_idxs;
    s_idx := Index(keys_fs,-1);
    genus_list := [curves[keys_fs[i]]`g : i in k_idxs];
    kernels :=[* *];
    for j->idx in k_idxs do
        g := genus_list[j];
        require #ds ge 2*g+3 : "We don't have enough values computed to determine the curve equation of curve", keys_fs[idx];
        //add one because one value will be the infinite value and useless
        //now remove the infinite column
        inf_idx_y2 := Index(table[idx], Infinity());
        inf_idx_s := Index(table[s_idx], Infinity());
        require inf_idx_y2 eq inf_idx_s : "y^2 and s have poles in different places";
        y2vals := Remove(table[idx], inf_idx_y2);
        svals := Remove(table[s_idx],inf_idx_s);
        rat_svals := [s :  i->s in svals   | Type(s) ne RngUPolElt];
        rat_sidxs := [i :  i->s in svals   | Type(s) ne RngUPolElt];
        rat_y2vals := [ y2vals[i]  : i in rat_sidxs];

        M := [];
        for j->s in rat_svals do
            Append(~M, [Rationals()!(s)^i : i in [0..2*g+2]] cat [Rationals()!rat_y2vals[j]]);
        end for;
        M := Matrix(M);
        B :=  Basis(Kernel(Transpose(M)));
        require #B le #ds - #rat_svals : "We don't have enough constraints imposed by the rational points";
        Append(~kernels, B);
    end for;
    return kernels;
end intrinsic;


function power_trace(trace, norm, j)
    //given trace of alpha^1 get trace alpha^j
    tr_pvs := 2; //trace of 1 is 2
    tr_curr := trace;
    if j eq 0 then
        return tr_pvs;
    elif j eq 1 then
        return tr_curr;
    else
        for i in [2..j] do
            tr_next :=  trace*tr_curr-tr_pvs*norm;
            tr_pvs := tr_curr;
            tr_curr := tr_next;
        end for;
        return tr_next;
    end if;

end function;

intrinsic QuadraticConstraintsOnEquations(schofer_table::SchoferTable, curves::SeqEnum, kernels::List) -> SeqEnum
    {}
    table := schofer_table`Values;
    keys_fs := schofer_table`Keys_fs;

    k_idxs := schofer_table`K_idxs;
    s_idx := Index(keys_fs,-1);

    all_relns := [* *];
    for j->idx in k_idxs do
        B := kernels[j];
        require not IsEmpty(B) : "Error in Schofer table values at rational points - no solution found!";
        numcoeffs := #Eltseq(B[1]);
        R<[x]>:=PolynomialRing(Rationals(),#B);
        inf_idx_y2 := Index(table[idx], Infinity());
        inf_idx_s := Index(table[s_idx], Infinity());
        y2vals := Remove(table[idx], inf_idx_y2);
        svals := Remove(table[s_idx],inf_idx_s);
        quad_svals := [s :  i->s in svals   | Type(s) eq RngUPolElt ];
        quad_sidxs := [i :  i->s in svals   | Type(s) eq RngUPolElt ];
        quad_y2vals := [ y2vals[i]  : i in quad_sidxs];
        coeff_list := [ &+[x[i]*B[i][k]: i in [1..#B]] : k in [1 .. numcoeffs]]; //this is a list of a_k, and also y^2
        

        relns := [];
        for l->val in quad_svals do 
            trace := -Coefficient(val, 1);
            nm := Coefficient(val, 0);
            rel :=0;
            for i in [0 .. numcoeffs-2] do //omit the y2 var
                for k in [0 .. numcoeffs-2] do //same
                    if i eq k then
                        rel +:= coeff_list[i+1]*coeff_list[k+1]*nm^i;
                    else
                        if i lt k then
                            rel +:= coeff_list[i+1]*coeff_list[k+1]*nm^i*power_trace(trace, nm, k-i);
                        else
                            continue;
                        end if;
                    end if;
                end for;
            end for;
            Append(~relns, rel-quad_y2vals[l] *coeff_list[numcoeffs]^2);
        end for;
        Append(~all_relns, relns);
    end for;
    return all_relns;
end intrinsic;



function solve_quadratic_constraints(relns)
    R := Universe(relns);
    P := ProjectiveSpace(R);
    S := Scheme(P, relns);
    assert Dimension(S) eq 0;
    return RationalPoints(S);
end function;


intrinsic EquationsOfCovers(schofer_table::SchoferTable, all_cm_pts::SeqEnum) -> SeqEnum, Assoc, SeqEnum
{Determine the equations of the covers using the values from Schofers formula}
    R<x> := PolynomialRing(Rationals());
    eqn_list := [ ];
    keys_fs := schofer_table`Keys_fs;
    table := schofer_table`Values;
    ds := schofer_table`Discs;
    k_idxs := schofer_table`K_idxs;
    curves := schofer_table`Curves;

    kernels := RationalConstraintsOnEquations(schofer_table, curves);
    relns := QuadraticConstraintsOnEquations(schofer_table, curves, kernels);

    for i->B in kernels do //indexed by k_idxs
        if #relns[i] eq 0 or #B eq 1 then
            require #B eq 1 : "Try adding quadratic points -- not enough constraints from the rational points";
            v :=  Eltseq(B[1]);
            monic_v := [-v[i]/v[#v] : i in [1..#v-1]];
            f := R!monic_v;
            Append(~eqn_list, f);
        else
            coeffs := solve_quadratic_constraints(relns[i]);
            require #coeffs eq 1 : "We do not have enough constraints coming from quadratic and rational points";
            coeffs := Eltseq(coeffs[1]);
            v := &+[B[j]*coeffs[j] : j in [1..#B]];
            v := Eltseq(v);
            monic_v := [-v[i]/v[#v] : i in [1..#v-1]];
            f := R!monic_v;
            Append(~eqn_list, f);
        end if;
    end for;


    crv_list := [HyperellipticCurve(e) : e in eqn_list];
    keys :=  [keys_fs[i] : i in k_idxs];
    Xstar := schofer_table`Xstar;
    all_W := Xstar`W;

    ws := AssociativeArray();
    for i->k in keys do
        my_W := curves[k]`W;
        nontriv := all_W diff my_W;
        ws[k] := AssociativeArray();
        for w in my_W do
            ws[k][w] := IdentityMap(crv_list[i]);
        end for;
        w := Representative(nontriv);
        for w in nontriv do
            ws[k][w] := HyperellipticInvolution(crv_list[i]);
        end for;
    end for;

    return crv_list, ws, keys;
end intrinsic;

intrinsic EquationsOfCovers(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot] : Prec := 100) -> SeqEnum, Assoc, SeqEnum
{Determine the equations of the immediate covers of X.}
    fs := BorcherdsForms(Xstar, curves : Prec := Prec);
    d_divs := &cat[[T[1]: T in DivisorOfBorcherdsForm(f, Xstar)] : f in [fs[-1], fs[-2]]]; //include zero infinity of hauptmoduls
    all_cm_pts := CandidateDiscriminants(Xstar, curves); // !!! This is slow, figure out why !!!
    genus_list := [curves[i]`g: i in Xstar`CoveredBy];
    
    // num_vals := Maximum([2*g+4 : g in genus_list]); // This is what we need for the equation part, but
    num_vals := Maximum([2*g+5 : g in genus_list]); // This is what we need for finding the y2 scales
    // Note that y^2 may vanish at 2*g+2 CM points, and be infinity at another one (2g+3).
    // We would need two other CM pts to determine the correct scaling, based on the fields of definition. 
    abs_schofer_tab, all_cm_pts := AbsoluteValuesAtCMPoints(Xstar, curves, all_cm_pts, fs : 
                                                            MaxNum := num_vals, Prec := Prec, 
                                                            Exclude := {}, Include := Set(d_divs));
    ReduceTable(abs_schofer_tab);
    schofer_tab := ValuesAtCMPoints(abs_schofer_tab, all_cm_pts);
    return EquationsOfCovers(schofer_tab, all_cm_pts);
end intrinsic;

function curves_above_P1_and_conics(crv_eqns, labels, curves)
    P1s := AssociativeArray();
    conics := AssociativeArray();
    for label in labels do
        eqns := crv_eqns[label];
        P1s[label] := {base : base in Keys(eqns) | Degree(crv_eqns[label][base]) eq 1};
        conics[label] := {base : base in Keys(eqns) | Degree(crv_eqns[label][base]) eq 2 and HasRationalPoint(Conic(crv_eqns[label][base]))};
    end for;

    curves_above_P1s := AssociativeArray();
    curves_above_conics := AssociativeArray();

    for P1_label in Keys(P1s) do
        for c in curves[P1_label]`CoveredBy do
            if not IsDefined(curves_above_P1s, c) then
                curves_above_P1s[c] := AssociativeArray();
            end if;
            curves_above_P1s[c][P1_label] := P1s[P1_label];
        end for;
    end for;

    for conic_label in Keys(conics) do
        for c in curves[conic_label]`CoveredBy do
            if not IsDefined(curves_above_conics, c) then
                curves_above_conics[c] := AssociativeArray();
            end if;
            curves_above_conics[c][conic_label] := conics[conic_label];
        end for;
    end for;

    return curves_above_P1s, curves_above_conics;
end function;

function equation_above_P1(covered_gplus1, covered_P1)
    fpoly := HyperellipticPolynomials(covered_P1);
    c0 := Coefficient(fpoly,0);
    c1 := Coefficient(fpoly,1);
    _<x>:=Parent(fpoly);
    eqn := HyperellipticPolynomials(covered_gplus1);
    gcd_poly := GCD(eqn, fpoly);
    eqn div:= gcd_poly;
    eqn := Evaluate(eqn, (x^2 - c0)/c1);
    C := HyperellipticCurve(eqn);
    return C, gcd_poly;
end function;

function ws_above_P1(H, label, P1_label, gplus1_label, curves, common_base, crv_ws)
    ws_label := AssociativeArray();
    C := Domain(crv_ws[P1_label][common_base][1]);
    Cpoly := HyperellipticPolynomials(C);
    c0 := Coefficient(Cpoly,0);
    c1 := Coefficient(Cpoly,1);
    P1<s,t> := ProjectiveSpace(Rationals(), 1);
    _<x,y,z> := Ambient(C);
    C_to_P1 := iso<C -> P1 | [y, z], [s^2-c0*t^2, c1*s*t, c1*t^2]>;
    deg_H := Degree(H);
    for m in Keys(crv_ws[P1_label][common_base]) do
        wm_crv := crv_ws[P1_label][common_base][m];
        assert m in Keys(crv_ws[gplus1_label][common_base]);
        wm_gplus1 := crv_ws[gplus1_label][common_base][m];
        deg_f := Degree(Domain(wm_gplus1));
        wm_P1 := Inverse(C_to_P1)*wm_crv*C_to_P1;
        wm_param := wm_P1(s)/wm_P1(t);
        _<x,y,z> := Ambient(H);
        wm_x_top := Evaluate(wm_param, [x,z]);
        wm_y_top := Evaluate(wm_gplus1(y)/wm_gplus1(z)^((deg_f+1) div 2), [((x/z)^2 - c0)/c1, y/z^((deg_H + 1) div 2), 1]);
        ws_label[m] := iso< H -> H | [wm_x_top, wm_y_top, 1], [wm_x_top, wm_y_top, 1]>;
    end for;
    return ws_label;
end function;

function equation_above_conic(covered_gplus1, covered_conic)
    C := Conic(covered_conic);
    P2<x,y,z> := Ambient(C);
    assert HasRationalPoint(C); // for now not implemented if C does not have a rational point
    P1_to_C := Parametrization(C);
    C_to_P1 := Inverse(P1_to_C);
    s_param := C_to_P1(x) / C_to_P1(z); // conic was constructed such that this is the hauptmodul
    fpoly := HyperellipticPolynomials(covered_gplus1);
    Cpoly := HyperellipticPolynomials(covered_conic);
    gcd_poly := GCD(fpoly, Cpoly);
    _<t> := PolynomialRing(Rationals());
    s_of_t := Evaluate(s_param, [t,1]);
    s_num := Numerator(s_of_t);
    s_denom := Denominator(s_of_t);
    // homogenized gcd poly, after substituion -
    // if s = F(t)/G(t), then hom_gcd = G^deg(gcd_poly)*gcd_poly(F(t)/G(t))
    hom_gcd := s_denom^Degree(gcd_poly)*Evaluate(Evaluate(gcd_poly, s_param), [t, 1]);
    assert Denominator(hom_gcd) eq 1;
    hom_gcd := Numerator(hom_gcd);
    lc := LeadingCoefficient(hom_gcd);
    monic_hom_gcd := hom_gcd / lc;
    is_sqr, root_gcd := IsSquare(monic_hom_gcd);
    error if not is_sqr, "GCD of hyperelliptic polynomials is not a square";
    f_prime := fpoly div gcd_poly;
    target_poly := s_denom^Degree(f_prime)*Evaluate(Evaluate(f_prime, s_param), [t, 1]);
    assert Denominator(target_poly) eq 1;
    target_poly := Numerator(target_poly);
    eps := Degree(fpoly) mod 2;
    target_poly := target_poly * s_denom^eps;
    // We have (y*s_denom^([d/2])/root_gcd)^2 = lc*target_poly,
    // where y is the y variable of covered_gplus1, and d is the degree (g+1) of the polynomial.
    H := HyperellipticCurve(lc*target_poly);
    y_factor := s_denom^((Degree(fpoly) + 1) div 2) / root_gcd;
    return H, C_to_P1, y_factor;
end function;


function ws_above_conic(H, C_to_P1, y_factor, label, conic_label, gplus1_label, curves, common_base, crv_ws)
    ws_label := AssociativeArray();
    deg_H := Degree(H);
    for m in Keys(crv_ws[conic_label][common_base]) do
        wm_conic := crv_ws[conic_label][common_base][m];
        assert m in Keys(crv_ws[gplus1_label][common_base]);
        wm_gplus1 := crv_ws[gplus1_label][common_base][m];
        deg_f := Degree(Domain(wm_gplus1));
        wm_P1 := Inverse(C_to_P1)*wm_conic*C_to_P1;
        _<s,t>  := Codomain(C_to_P1);
        wm_param := wm_P1(s)/wm_P1(t);
        _<x,y,z> := Ambient(H);
        wm_x_top := Evaluate(wm_param, [x,z]);
        wm_y_factor := Evaluate(y_factor, Evaluate(wm_param, [x/z,1]))/Evaluate(y_factor, x/z);
        s_param := C_to_P1(x) / C_to_P1(z);
        wm_y_top := wm_y_factor*Evaluate(wm_gplus1(y)/wm_gplus1(z)^((deg_f+1) div 2), [Evaluate(s_param, [x/z,1]), y/z^((deg_H + 1) div 2), 1]);
        ws_label[m] := iso< H -> H | [wm_x_top, wm_y_top, 1], [wm_x_top, wm_y_top, 1]>;
    end for;
    return ws_label;
end function;

function process_P1_cover(label, curves_above_P1s, curves, crv_eqns, crv_ws)
    vprintf ShimuraQuotients, 3 : "\n\t\tProcessing curve %o covering a P1...", label;
    g := curves[label]`g;
    crv_eqns_label := AssociativeArray();
    crv_ws_label := AssociativeArray();
    covered_curves := Keys(crv_eqns) meet curves[label]`Covers;
    for P1_label in Keys(curves_above_P1s[label]) do
        covered_others := covered_curves diff {P1_label};
        for other_label in covered_others do
            common_bases := Keys(crv_eqns[other_label]) meet curves_above_P1s[label][P1_label];
            found_base := false;
            for base in common_bases do
                if (Degree(HyperellipticPolynomials(crv_eqns[other_label][base])) eq g+1) then
                    common_base := base;
                    found_base := true;
                    break;
                end if;
            end for;
            if found_base then
                gplus1_label := other_label;
                break;
            end if;
        end for;
        if not found_base then continue; end if;
        covered_P1 := crv_eqns[P1_label][common_base];
        covered_gplus1 := crv_eqns[gplus1_label][common_base];
        H := equation_above_P1(covered_gplus1, covered_P1);
        vprintf ShimuraQuotients, 3 : "Found equation.";
        ws_H := ws_above_P1(H, label, P1_label, gplus1_label, curves, common_base, crv_ws);
        crv_eqns_label[P1_label] := H;
        crv_ws_label[P1_label] := ws_H;
    end for;
    return crv_eqns_label, crv_ws_label;
end function;

function process_conic_cover(label, curves_above_conics, curves, crv_eqns, crv_ws)
    vprintf ShimuraQuotients, 3 : "\n\t\tProcessing curve %o covering a conic...", label;
    g := curves[label]`g;
    crv_eqns_label := AssociativeArray();
    crv_ws_label := AssociativeArray();
    covered_curves := Keys(crv_eqns) meet curves[label]`Covers;
    for conic_label in Keys(curves_above_conics[label]) do
        covered_others := covered_curves diff {conic_label};
        for other_label in covered_others do
            common_bases := Keys(crv_eqns[other_label]) meet curves_above_conics[label][conic_label];
            found_base := false;
            for base in common_bases do
                if (Degree(HyperellipticPolynomials(crv_eqns[other_label][base])) eq g+1) 
                    or (g eq 0 and other_label in Keys(curves_above_conics[label])) then
                    common_base := base;
                    found_base := true;
                    break;
                end if;
            end for;
            if found_base then
                gplus1_label := other_label;
                break;
            end if;
        end for;
        if not found_base then continue; end if;
        covered_conic := crv_eqns[conic_label][common_base];
        covered_gplus1 := crv_eqns[gplus1_label][common_base];
        H, C_to_P1, yfactor := equation_above_conic(covered_gplus1, covered_conic);
        vprintf ShimuraQuotients, 3 : "Found equation.";
        ws_H := ws_above_conic(H, C_to_P1, yfactor, label, conic_label, gplus1_label, curves, common_base, crv_ws);
        crv_eqns_label[conic_label] := H;
        crv_ws_label[conic_label] := ws_H;
    end for;
    return crv_eqns_label, crv_ws_label;
end function;


intrinsic EquationsAboveP1s(crv_list::SeqEnum[CrvHyp], ws::Assoc, keys::SeqEnum[RngIntElt], curves::SeqEnum[ShimuraQuot]) -> Assoc, Assoc
{Using Riemann Roch, leverage covered equations to get higher cover equations}

    // initializing data structures to store for each curve all equations and the corresponding ws,
    // depending on the P1s (and conics) that it covers
    crv_eqns := AssociativeArray();
    crv_ws := AssociativeArray();
    for i->key in keys do
        crv_eqns[key] := AssociativeArray();
        crv_ws[key] := AssociativeArray();
        covered_P1s := curves[key]`Covers;
        assert #covered_P1s eq 1;
        star_key := Representative(covered_P1s);
        crv_eqns[key][star_key] := crv_list[i];
        crv_ws[key][star_key] := ws[key];
    end for;

    curves_above_P1s, curves_above_conics := curves_above_P1_and_conics(crv_eqns, keys, curves);

    new_keys := Keys(curves_above_P1s) join Keys(curves_above_conics);

    while not IsEmpty(new_keys) do
        vprintf ShimuraQuotients, 2 : "\n\tRemaining curves above P1s: %o,", Keys(curves_above_P1s);
        vprintf ShimuraQuotients, 2 : "\n\tRemaining curves above conics: %o...", Keys(curves_above_conics);
        for label in Keys(curves_above_P1s) do
            crv_eqns_label, crv_ws_label := process_P1_cover(label, curves_above_P1s, curves, crv_eqns, crv_ws);
            crv_eqns[label] := crv_eqns_label;
            crv_ws[label] := crv_ws_label;
        end for;

        for label in Keys(curves_above_conics) do
            crv_eqns_label, crv_ws_label := process_conic_cover(label, curves_above_conics, curves, crv_eqns, crv_ws);
            for base in Keys(crv_eqns_label) do
                crv_eqns[label][base] := crv_eqns_label[base];
                crv_ws[label][base] := crv_ws_label[base];
            end for;
        end for;

        curves_above_P1s, curves_above_conics := curves_above_P1_and_conics(crv_eqns, new_keys, curves);
        new_keys := Keys(curves_above_P1s) join Keys(curves_above_conics);
    end while;
    vprintf ShimuraQuotients, 2 : "\n";
    return crv_eqns, crv_ws;

end intrinsic;

intrinsic EquationsAbovePointlessConics(all_eqns::Assoc, all_ws::Assoc, curves::SeqEnum : base_label := 0) -> Assoc, Assoc
    {Find equations above pointless conics, as a last step}
    all_keys := Keys(all_eqns);
    not_done := [k : k in all_keys | #Keys(all_eqns[k]) eq 0]; // don't have an equation over anything they cover
    starcurve := Representative(curves[Maximum(all_keys)]`Covers);
    assert IsStarCurve(curves[starcurve]);
    known_conics := [k : k in all_keys | k notin not_done and curves[k]`g eq 0]; 
    //now find all curves lying over a conic without an equation
    curves_to_do := [k : k in not_done | #(curves[k]`Covers meet Set(known_conics)) gt 0];
    for k in curves_to_do do
        g := curves[k]`g;
        assert exists(conic_key){x : x in (curves[k]`Covers meet Set(known_conics))}; //find the conic that it covers
        for other_curve in curves[k]`Covers do
            found_gplus1 := false;
            bases := Keys(all_eqns[other_curve]);
            if IsEmpty(bases) then continue; end if;// all eqns for all bases have the same degree
            if (base_label ne 0) and base_label notin bases then continue; end if;
            base := (base_label eq 0) select Representative(bases) else base_label;
            if (Degree(HyperellipticPolynomials(all_eqns[other_curve][base])) eq g+1) then
                gplus1key := other_curve; //found the gplus1
                found_gplus1 := true;
                break;
            end if;
        end for;
        if not found_gplus1 then continue; end if;

        //combine equations to get the equation for the curve
        covered_gplus1 := all_eqns[gplus1key][base];
        covered_conic:= all_eqns[conic_key][base];
        wt := curves[gplus1key]`g +1;
        P3<x,y,s,z> := WeightedProjectiveSpace(Rationals(),[1,wt,1,1]);
        fconic := HyperellipticPolynomials(covered_conic);
        fgplus1 := HyperellipticPolynomials(covered_gplus1);
        eqn1 := Evaluate(fconic, s);
        eqn2 := Evaluate(fgplus1,s);
        eqn2 := Homogenization(eqn2,z);
        eqn1 := Homogenization(eqn1,z);
        C := Curve(P3, [y^2 - eqn2, x^2 - eqn1]);
        all_eqns[k][conic_key] := C;

        all_ws[k][conic_key] := AssociativeArray(); 
        //now find the ws
        ws_to_do := Keys(all_ws[gplus1key][base]);
        for w in ws_to_do do
            w_mapgplus1 := all_ws[gplus1key][base][w];
            alg_map_gplus1 := AlgebraMap(w_mapgplus1);
            y_var := Domain(alg_map_gplus1).2;
            s_var := Domain(alg_map_gplus1).1;
            z_var := Domain(alg_map_gplus1).3;
            im_y := Evaluate(alg_map_gplus1(y_var), [s,y,z]);
            im_s1 := Evaluate(alg_map_gplus1(s_var), [s,y,z]);
            im_z1 := Evaluate(alg_map_gplus1(z_var), [s,y,z]);
            w_mapconic := all_ws[conic_key][base][w];
            alg_map_conic := AlgebraMap(w_mapconic);
            s_var := Domain(alg_map_conic).1;
            x_var := Domain(alg_map_conic).2;
            z_var := Domain(alg_map_conic).3;
            im_x :=  Evaluate(alg_map_conic(x_var), [s,x,z]);
            im_s2 := Evaluate(alg_map_conic(s_var), [s,y,z]);
            im_z2 := Evaluate(alg_map_conic(z_var), [s,y,z]);
            assert im_s1*im_z2 eq im_s2*im_z1; //acts the same way
            wmap := map<C->C | [im_x, im_y, im_s1, im_z1]>;
            b, inv := IsIsomorphism(wmap);
            assert b;
            all_ws[k][conic_key][w] := inv^(-1);
        end for;
    end for;
    return all_eqns, all_ws;

end intrinsic;

intrinsic AllEquationsAboveCovers(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot] : Prec := 100, base_label := 0)-> Assoc, Assoc
{Get equations of all covers (not just immediate covers)}
    require IsStarCurve(Xstar): "Xstar must be a star curve";
    vprintf ShimuraQuotients, 1 : "Computing Borcherds forms...";
    fs := BorcherdsForms(Xstar, curves : Prec := Prec);
    vprintf ShimuraQuotients, 1 : "Done!\n";
    vprintf ShimuraQuotients, 1 : "Computing divisors of hauptmodules...";
    d_divs := &cat[[T[1]: T in DivisorOfBorcherdsForm(f, Xstar)] : f in [fs[-1], fs[-2]]]; //include zero infinity of hauptmoduls
    vprintf ShimuraQuotients, 4 : "\n";
    vprintf ShimuraQuotients, 1 : "Done!\n";
    vprintf ShimuraQuotients, 1 : "Computing candidate discriminants...";
    all_cm_pts := CandidateDiscriminants(Xstar, curves);
    genus_list := [curves[i]`g: i in Xstar`CoveredBy];
    num_vals := Maximum([2*g+5 : g in genus_list]);
    vprintf ShimuraQuotients, 1 : "Computing absolute values at CM points...";
    abs_schofer_tab, all_cm_pts:= AbsoluteValuesAtCMPoints(Xstar, curves, all_cm_pts, fs : 
                                                           MaxNum := num_vals, Prec := Prec, 
                                                           Exclude := {}, Include := Set(d_divs));
    vprintf ShimuraQuotients, 2 : "\n";
    vprintf ShimuraQuotients, 1 : "Done!\n";
    ReduceTable(abs_schofer_tab);
    vprintf ShimuraQuotients, 1 : "Computing actual values at CM points...";
    schofer_tab := ValuesAtCMPoints(abs_schofer_tab, all_cm_pts : Exclude := {});
    vprintf ShimuraQuotients, 1 : "Done!\n";  
    vprintf ShimuraQuotients, 1 : "Computing equations of covers...";
    crv_list, ws, new_keys := EquationsOfCovers(schofer_tab, all_cm_pts);
    vprintf ShimuraQuotients, 1 : "Done!\n";
    vprintf ShimuraQuotients, 1 : "Computing equations above P1s and conics...";
    all_eqns, all_ws := EquationsAboveP1s(crv_list, ws, new_keys, curves); //still adding ws here in the conic case
    vprintf ShimuraQuotients, 1 : "Done\n";
    vprintf ShimuraQuotients, 1 :"Computing equations above pointless conics...";
    all_eqns, all_ws := EquationsAbovePointlessConics(all_eqns, all_ws, curves : base_label := base_label);
    vprintf ShimuraQuotients, 1 : "Done\n";
    return all_eqns, all_ws;
end intrinsic;

intrinsic AllEquationsAboveCovers(D::RngIntElt, N::RngIntElt, curves::SeqEnum[ShimuraQuot] : Prec := 100, base_label := 0)-> Assoc, Assoc
{Get equations of all covers (not just immediate covers)}
    _ := exists(Xstar){X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)};
    return AllEquationsAboveCovers(Xstar, curves : Prec := Prec, base_label := base_label);
end intrinsic;