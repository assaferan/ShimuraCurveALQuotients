// Compute the order of vanishing of eta(delta tau) at the cusp a/b in Gamma_0(N)
// multiplied by 24
function OrderOfVanishingOfEta(delta, b, N)
    return N*GCD(b,delta)^2 div (GCD(N,b^2)*delta);
end function;

function what_we_wrote()
    Genus(Gamma0(60));
    Cusps(Gamma0(60));
    #Cusps(Gamma0(60));
    M := 60;
    Divisors(M);
    ds := Divisors(M);
    #ds;
    cond_1 := Matrix(Integers(),[[1 : i in [1..#ds]]]);
    cond_1;
    ps := PrimeDivisors(M);
    ps;
    cond_2 := Matrix([[Valuation(d,p) : d in Divisors(M)] : p in ps]);
    cond_2;
    cond_1;
    cond_3 := Matrix(Integers(),ds);
    cond_3 := Matrix(Integers(),[d : d in ds]);
    cond_3 := Matrix(Integers(),[ds]);
    cond_3;
    cond_4 := Matrix(Integers(),[[M div d : d in ds]]);
    cond_4;
    VerticalJoin([cond_1, cond_2, cond_3, cond_4]);
    VerticalJoin([* cond_1, cond_2, cond_3, cond_4 *]);
    VerticalJoin(cond_1, cond_2);
    VerticalJoin(VerticalJoin(cond_1, cond_2), cond_3);
    VerticalJoin(VerticalJoin(VerticalJoin(cond_1, cond_2), cond_3), cond_4);
    mat := VerticalJoin(VerticalJoin(VerticalJoin(cond_1, cond_2), cond_3), cond_4);
    mat;
    cusps := [1/d : d in ds];
    cusps;
    function OrderOfVanishingOfEta(delta, b, N)
        return N*GCD(b,delta)^2 div (GCD(N,b^2)*delta);
    end function;
    eq_8 := Matrix([[OrderOfVanishingOfEta(delta, b, M) : b in ds] : delta in ds]);
    eq_8;
    return 0;
end function;

function get_D0_M_g(D, N)
    assert IsEven(D) and IsSquarefree(N);
    D0 := (D*N) div 2^Valuation(D,2);
    M := 4*D0; // this is 2*D*N
    g := Genus(Gamma0(M));
    return D0,M,g;
end function;

function lhs_integer_programming(M)
    ds := Divisors(M);
    ps := PrimeDivisors(M);
    cond_1 := Matrix(Integers(),[[1 : i in [1..#ds]]]);
    cond_2 := Matrix([[Valuation(d,p) : d in Divisors(M)] : p in ps]);
    cond_3 := Matrix(Integers(),[ds]);
    cond_4 := Matrix(Integers(),[[M div d : d in ds]]);
    eq_7 := VerticalJoin(VerticalJoin(VerticalJoin(cond_1, cond_2), cond_3), cond_4);
    eq_8 := Matrix([[OrderOfVanishingOfEta(delta, b, M) : delta in ds] : b in ds]);
    delta_epsilon_n := DiagonalMatrix([-2 : p in ps] cat [-24,-24,0]);
    ncols:= Ncols(delta_epsilon_n);
    r_coeffs := VerticalJoin(eq_7, eq_8);
    zerotop := Matrix(Integers(), 1,ncols, [0 : i in [1..ncols]]);
    nrows := Nrows(r_coeffs)-1-Nrows(delta_epsilon_n);
    zerobottom := Matrix(Integers(), nrows ,ncols, [0 : i in [1..nrows*ncols]]);
    add_coeffs := VerticalJoin(zerotop, VerticalJoin(delta_epsilon_n, zerobottom));
    lhs := HorizontalJoin(r_coeffs, add_coeffs);
    lhs := VerticalJoin(lhs, Matrix([lhs[Nrows(lhs)]])); // making sure that it is not holomorphic at infinity
    lhs[Nrows(lhs)-1,Ncols(lhs)] := 24; // bounding the pole order by -n
    lhs := VerticalJoin(lhs, Matrix([[0 : i in [1..Ncols(lhs)-1]] cat [1]])); // bound n from below (by 1)
    return lhs;
end function;

function integer_programming_input(D,N)
    D0,M,g := get_D0_M_g(D, N);
    ds := Divisors(M);
    ps := PrimeDivisors(M);
    lhs := lhs_integer_programming(M);
    n_eq := #ps + 3;
    rels_seq := [0 : i in [1..n_eq]] cat [1 : i in [1..#ds]] cat [-1,1];
    rels := Matrix(Integers(), Nrows(lhs), 1, rels_seq);
    rhs := Matrix(Integers(), Nrows(lhs), 1, [0 : i in [1..n_eq + #ds]] cat [-1,1]);
    objective := Matrix(Integers(), 1, Ncols(lhs), [0 : i in [1..Ncols(lhs)-1]] cat [1]);
    // This is what we want but doesn't work because assumes all variables are nonnegative
    // MinimalIntegerSolution(lhs, rels, rhs, objective);
    LP := LPProcess(Integers(),Ncols(lhs));
    AddConstraints(LP,Matrix(lhs[1..n_eq]), Matrix(rhs[1..n_eq]) : Rel := "eq");
    AddConstraints(LP,Matrix(lhs[n_eq + 1..n_eq + #ds]), Matrix(rhs[n_eq + 1..n_eq + #ds]) : Rel := "ge");
    idx := n_eq + #ds + 1;
    assert idx + 1 eq Nrows(lhs);
    AddConstraints(LP,Matrix(lhs[idx..idx]), Matrix(rhs[idx..idx]) : Rel := "le");
    AddConstraints(LP,Matrix(lhs[idx+1..idx+1]), Matrix(rhs[idx+1..idx+1]) : Rel := "ge");
    for n in [1..Ncols(lhs)] do
        SetLowerBound(LP,n,-100);
    end for;
    SetObjectiveFunction(LP, objective); 
    t := Solution(LP);
    // need to define g and D0
    n0 := Maximum(2*g-2 - &+[d div 4 : d in Divisors(D0)],0);
    k := t[1,Ncols(t)]; // the order of pole for t
    n := n0 + k;
    lhs := Submatrix(lhs, [1..Nrows(lhs)-1], [1..Ncols(lhs)-1]);
    rhs := Submatrix(rhs, [1..Nrows(rhs)-1], [1..1]);
    rhs[Nrows(rhs)-1,1] := -24*n;
    rhs[1,1] := 1; // admissibility condition
    rhs[2,1] := 1; // discriminant is 2 times a square - equation corresponding to valuation at 2
    // eqs in a format for polymake / sage
    eqs := [Eltseq(row) : row in Rows(HorizontalJoin(-Matrix(rhs[1..n_eq]),Matrix(lhs[1..n_eq])))];
    // inequalities in format for polymake / sage
    ieqs := [Eltseq(row) : row in Rows(HorizontalJoin(-Matrix(rhs[n_eq + 1.. n_eq + #ds]),Matrix(lhs[n_eq+1..n_eq+#ds])))];
    return eqs, ieqs;
    /*
    LP := LPProcess(Integers(),Ncols(lhs));
    AddConstraints(LP,Matrix(lhs[1..n_eq]), Matrix(rhs[1..n_eq]) : Rel := "eq");
    AddConstraints(LP,Matrix(lhs[n_eq + 1..n_eq + #ds]), Matrix(rhs[n_eq + 1..n_eq + #ds]) : Rel := "ge");
    idx := n_eq + #ds + 1;
    assert idx eq Nrows(lhs);
    AddConstraints(LP,Matrix(lhs[idx..idx]), Matrix(rhs[idx..idx]) : Rel := "le");
    for j in [1..Ncols(lhs)] do
        SetLowerBound(LP,j,-100);
    end for;
    // using Adam's code to find all solutions
    R := Integers();
    d := Ncols(lhs);
    L := LP;
    rels := [0 : i in [1..n_eq]] cat [1 : i in [1..#ds]] cat [-1];
    intbds := [];
    for i in [1..d] do
        SetObjectiveFunction(L,Matrix([[j eq i select R!1 else R!0: j in [1..d]]]));
        SetMaximiseFunction(L,false);
        sol,val := Solution(L); assert val eq 0;
        lb := sol[1,i];
        SetMaximiseFunction(L,true);
        sol,val := Solution(L); assert val eq 0;
        ub := sol[1,i];
        Append(~intbds,[Ceiling(lb-10^-10),Floor(ub+10^-10)]);
    end for;
    */
    /*
    > intbds;
    [
    [ -5, 3 ],
    [ -3, 12 ],
    [ -4, 6 ],
    [ -5, 3 ],
    [ 2, 15 ],
    [ -10, -4 ],
    [ -4, 0 ],
    [ -1, 0 ],
    [ -11, -3 ],
    [ 0, 7 ]
    ]
    */
    /*
    box_elts := [];
    S := CartesianProduct([{i[1]..i[2]}: i in intbds]);
    j := 0;
    for pt in S do
        if (j mod 100000 eq 0) then
            print j, "/", #S;
        end if;
        v := Vector([pt[i] : i in [1..#pt]]);
        res := v*Transpose(lhs) - Vector(Transpose(rhs));
        if [Sign(res[i]) : i in [1..Degree(res)]] eq rels then
            Append(~box_elts, v);
        end if;
        j +:= 1;
    end for;
    elim := Matrix(lhs[1..n_eq]);
    non_elim := Matrix(lhs[n_eq+1..Nrows(lhs)]);
    H,T := HermiteForm(elim);
    lhs_sub := non_elim - Submatrix(non_elim,[1..Nrows(non_elim)],[1..4])*Submatrix(H,[1..4],[1..Ncols(H)]);
    lhs_sub := Submatrix(lhs_sub, [1..Nrows(lhs_sub)], [5..Ncols(lhs_sub)]);
    elim_rhs := Matrix(rhs[1..n_eq]);
    H_rhs := T*elim_rhs;
    non_elim_rhs := Matrix(rhs[n_eq+1..Nrows(rhs)]);
    rhs_sub := non_elim_rhs - Submatrix(non_elim,[1..Nrows(non_elim)],[1..4])*Submatrix(H_rhs,[1..4],[1..Ncols(H_rhs)]);
    LP := LPProcess(Integers(),Ncols(lhs_sub));
    AddConstraints(LP,Matrix(lhs_sub[1..Nrows(lhs_sub)-1]), Matrix(rhs_sub[1..Nrows(lhs_sub)-1]) : Rel := "ge");
    idx := Nrows(lhs_sub);
    AddConstraints(LP,Matrix(lhs_sub[idx..idx]), Matrix(rhs_sub[idx..idx]) : Rel := "le");
    for j in [1..Ncols(lhs_sub)] do
        SetLowerBound(LP,j,-100);
    end for;
    // using Adam's code to find all solutions
    R := Integers();
    d := Ncols(lhs_sub);
    L := LP;
    rels := [0 : i in [1..n_eq]] cat [1 : i in [1..#ds]] cat [-1];
    intbds := [];
    for i in [1..d] do
        SetObjectiveFunction(L,Matrix([[j eq i select R!1 else R!0: j in [1..d]]]));
        SetMaximiseFunction(L,false);
        sol,val := Solution(L); assert val eq 0;
        lb := sol[1,i];
        SetMaximiseFunction(L,true);
        sol,val := Solution(L); assert val eq 0;
        ub := sol[1,i];
        Append(~intbds,[Ceiling(lb-10^-10),Floor(ub+10^-10)]);
    end for;
    box_elts := [];
    S := CartesianProduct([{i[1]..i[2]}: i in intbds]);
    j := 0;
    for pt in S do
        if (j mod 100000 eq 0) then
            print j, "/", #S;
        end if;
        v := Vector([pt[i] : i in [1..#pt]]);
        res := v*Transpose(lhs) - Vector(Transpose(rhs));
        if [Sign(res[i]) : i in [1..Degree(res)]] eq rels then
            Append(~box_elts, v);
        end if;
        j +:= 1;
    end for;
    return t;*/
end function;

procedure write_polymake_scriptfile(D,N)
    eqs, ieqs := integer_programming_input(D,N);
    output_lines := [];
    Append(~output_lines, "use application \"polytope\";");
    Append(~output_lines, "use vars '$ieqs', '$eqs', '$p';");
    Append(~output_lines, Sprintf("$ieqs = %o;", ieqs));
    Append(~output_lines, Sprintf("$eqs = %o;", eqs));
    Append(~output_lines, "$p = new Polytope(INEQUALITIES=>$ieqs, EQUATIONS=>$eqs);");
    Append(~output_lines, "print $p->LATTICE_POINTS;");
    output := Join(output_lines, "\n");
    fname := Sprintf("polymake_script_%o_%o", D, N);
    Write(fname, output : Overwrite);
    return;
end procedure;

function get_integer_prog_solutions(D,N)
    write_polymake_scriptfile(D,N);
    fname := Sprintf("polymake_script_%o_%o", D, N);
    polymake := Read(POpen("polymake --script " cat fname, "r"));
    sol_lines := Split(polymake, "\n");
    sol_vecs := [Split(line, " ") : line in sol_lines];
    sols := [[eval(x) : x in vec] : vec in sol_vecs];
    M := 2*D*N;
    rs := [sol[2..1 + #Divisors(M)] : sol in sols];
    return rs;
end function;

function get_weakly_holomorphic_basis(D,N)
    D0,M,g := get_D0_M_g(D,N);
    R<q> := PowerSeriesRing(Rationals());
    eta<q> := DedekindEta(q);  
    nor_eta := eta / q^(1/24);
    rs := get_integer_prog_solutions(D,N);
    M := 2*D*N;
    eta_quotients := [&*[(Evaluate(nor_eta,q^d)*q^(d/24))^r[i] : 
                        i->d in Divisors(M)] : r in rs];
    min_v := Minimum([Valuation(eta_quot) : eta_quot in eta_quotients]);
    min_prec := Minimum([RelativePrecision(eta_quot) - min_v + Valuation(eta_quot) : eta_quot in eta_quotients]);
    coeffs := Matrix([AbsEltseq(q^(-min_v)*(R!eta_quo) : FixedLength)[1..min_prec] : eta_quo in eta_quotients]);
    E, T := EchelonForm(coeffs);
    // sanity checks
    n := -min_v;
    dim := n + &+[d div 4 : d in Divisors(D0)] + 1 - g;
    n_gaps := g - &+[d div 4 : d in Divisors(D0)];
    assert Rank(E) eq dim;
    pole_orders := [PivotColumn(E,i) - n - 1 : i in [1..Rank(E)]];
    assert (n + 1 - #pole_orders) eq n_gaps;
    return E, n;
end function;