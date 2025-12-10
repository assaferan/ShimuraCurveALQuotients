// Compute the order of vanishing of eta(delta tau) at the cusp a/b in Gamma_0(N)
// multiplied by 24

function order_of_vanishing_of_eta(delta, b, N)
    return N*GCD(b,delta)^2 div (GCD(N,b^2)*delta);
end function;

function get_D0_M_g(D, N)
    // assert IsEven(D) and IsSquarefree(N);
    assert IsSquarefree(N);
    D0 := (D*N) div 2^Valuation(D,2);
    M := 4*D0;
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
    eq_8 := Matrix([[order_of_vanishing_of_eta(delta, b, M) : delta in ds] : b in ds]);
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

function create_lhs_rhs(M)
    ds := Divisors(M);
    ps := PrimeDivisors(M);
    lhs := lhs_integer_programming(M);
    n_eq := #ps + 3;
    // rels_seq := [0 : i in [1..n_eq]] cat [1 : i in [1..#ds]] cat [-1,1];
    // rels := Matrix(Integers(), Nrows(lhs), 1, rels_seq);
    rhs := Matrix(Integers(), Nrows(lhs), 1, [0 : i in [1..n_eq + #ds]] cat [-1,1]);
    return lhs, rhs, n_eq, #ds;
end function;

function find_t(M)
    lhs, rhs, n_eq, n_ds := create_lhs_rhs(M);
    objective := Matrix(Integers(), 1, Ncols(lhs), [0 : i in [1..Ncols(lhs)-1]] cat [1]);
    // This is what we want but doesn't work because assumes all variables are nonnegative
    // MinimalIntegerSolution(lhs, rels, rhs, objective);
    LP := LPProcess(Integers(),Ncols(lhs));
    AddConstraints(LP,Matrix(lhs[1..n_eq]), Matrix(rhs[1..n_eq]) : Rel := "eq");
    AddConstraints(LP,Matrix(lhs[n_eq + 1..n_eq + n_ds]), Matrix(rhs[n_eq + 1..n_eq + n_ds]) : Rel := "ge");
    idx := n_eq + n_ds + 1;
    assert idx + 1 eq Nrows(lhs);
    AddConstraints(LP,Matrix(lhs[idx..idx]), Matrix(rhs[idx..idx]) : Rel := "le");
    AddConstraints(LP,Matrix(lhs[idx+1..idx+1]), Matrix(rhs[idx+1..idx+1]) : Rel := "ge");
    for n in [1..Ncols(lhs)] do
        SetLowerBound(LP,n,-100);
    end for;
    SetObjectiveFunction(LP, objective); 
    t := Solution(LP);
    return t, lhs, rhs, n_eq, n_ds;
end function;

// n is the order of the pole at infty,
// m the order of the pole at 0
function integer_programming_input(lhs, rhs, n_eq, n_ds, n, m : k := 1/2, sq_disc := false, cuspidal := false)
    lhs := Submatrix(lhs, [1..Nrows(lhs)-1], [1..Ncols(lhs)-1]);
    rhs := Submatrix(rhs, [1..Nrows(rhs)-1], [1..1]);
    rhs[n_eq+1,1] := -24*m;
    rhs[Nrows(rhs)-1,1] := -24*n;
    rhs[1,1] := 2*k; // admissibility condition
    rhs[2,1] := (sq_disc select 0 else 1); // discriminant is 2 times a square - equation corresponding to valuation at 2
    if cuspidal then
        for j in [1..n_ds] do
            rhs[n_eq+j,1] := 1;
        end for;
    end if;
    // eqs in a format for polymake / sage
    eqs := [Eltseq(row) : row in Rows(HorizontalJoin(-Matrix(rhs[1..n_eq]),Matrix(lhs[1..n_eq])))];
    // inequalities in format for polymake / sage
    ieqs := [Eltseq(row) : row in Rows(HorizontalJoin(-Matrix(rhs[n_eq + 1.. n_eq + n_ds]),Matrix(lhs[n_eq+1..n_eq+n_ds])))];
    // return eqs, ieqs, t;
    return eqs, ieqs;
end function;

procedure write_polymake_scriptfile(M, lhs, rhs, n_eq, n_ds, n, m : k := 1/2, sq_disc := false, cuspidal := false)
    eqs, ieqs := integer_programming_input(lhs, rhs, n_eq, n_ds, n, m : k := k, sq_disc := sq_disc, cuspidal := cuspidal);
    output_lines := [];
    Append(~output_lines, "use application \"polytope\";");
    Append(~output_lines, "use vars '$ieqs', '$eqs', '$p';");
    Append(~output_lines, Sprintf("$ieqs = %o;", ieqs));
    Append(~output_lines, Sprintf("$eqs = %o;", eqs));
    Append(~output_lines, "$p = new Polytope(INEQUALITIES=>$ieqs, EQUATIONS=>$eqs);");
    Append(~output_lines, "print $p->LATTICE_POINTS;");
    output := Join(output_lines, "\n");
    fname := Sprintf("polymake/polymake_script_%o_%o_%o", M, n, m);
    Write(fname, output : Overwrite);
    return;
end procedure;

function get_integer_prog_solutions(M, lhs, rhs, n_eq, n_ds, n, m : k := 1/2, sq_disc := false, cuspidal := false)
    vprintf ShimuraQuotients, 3 : "\n\t\tMaking polymake file for (%o, %o, %o)...", M, n, m;
    if FileExists(Sprintf("polymake/polymake_solution_%o_%o_%o", M, n, m)) then
        vprintf ShimuraQuotients, 3 : "File found.";
        return eval Read(Sprintf("polymake/polymake_solution_%o_%o_%o", M, n, m));
    end if;
    vprintf ShimuraQuotients, 3 : "File not found, computing...";
    write_polymake_scriptfile(M, lhs, rhs, n_eq, n_ds, n, m : k := k, sq_disc := sq_disc, cuspidal := cuspidal);
    fname := Sprintf("polymake/polymake_script_%o_%o_%o", M, n, m);
    polymake := Read(POpen("polymake --script " cat fname cat " 2>/dev/null", "r"));
    if IsEof(polymake) then return []; end if;

    sol_lines := Split(polymake, "\n");
    sol_vecs := [Split(line, " ") : line in sol_lines];
    sols := [[eval(x) : x in vec] : vec in sol_vecs];
    rs := [sol[2..1 + #Divisors(M)] : sol in sols];

    Write(Sprintf("polymake/polymake_solution_%o_%o_%o", M, n, m), Sprint(rs, "Magma"));
    return rs;
end function;

intrinsic HolomorphicEtaQuotients(N::RngIntElt, k::RngIntElt : Prec := 100) -> SeqEnum[EtaQuot]
{Returns a basis of M_k(N) using eta quotients.}
    R := EtaQuotientsRing(N, 1);
    lhs, rhs, n_eq, n_ds := create_lhs_rhs(N);
    rs := get_integer_prog_solutions(N, lhs, rhs, n_eq, n_ds, 0, 0 : k := k, sq_disc := true); // , cuspidal := true);
    return [EtaQuotient(R,r) : r in rs];
end intrinsic;

function coeff_height(vec)
    return &+[Log(AbsoluteValue(Numerator(x)))+Log(AbsoluteValue(Denominator(x))) : x in Eltseq(vec) | x ne 0];
end function;

function find_minimal_expression_of_length(etas, l, mat, v)
    min_height := Infinity();
    found := false;
    for S in Subsets({1..#etas},l) do
        S_seq := [j : j in S];
        sub_mat := [mat[i] : i in S_seq];
        V := sub< RowSpace(mat) | sub_mat>;
        if v notin V then continue; end if;
        found := true;
        sol := Solution(Matrix(sub_mat), v);
        height := coeff_height(sol);
        if height lt min_height then
            eta := &+[sol[i]*etas[S_seq[i]] : i in [1..#S_seq]];
            min_height := height;
        end if;
    end for;
    if not found then return false, _; end if;
    return true, eta;
end function;

function find_short_expression_using_LP(etas, mat, v : M := 10, R := Rationals(), D := 1)

    // We will create a system in variables x_i, y_i
    LP := LPProcess(R, 2*#etas);
    
    // adding equations Ax = D*b
    lhs_eq := Transpose(ChangeRing(mat,R));
    lhs_eq := HorizontalJoin(lhs_eq, Parent(lhs_eq)!0);
    rhs_eq := Transpose(Matrix([D*ChangeRing(v,R)]));
    AddConstraints(LP, lhs_eq, rhs_eq : Rel := "eq");

    // assing constraints x_i <= D M y_i
    one := MatrixAlgebra(R, #etas)!1;
    lhs_le := HorizontalJoin(one, -D*M*one); 
    rhs_le := Matrix(R, [[0] : eta in etas]);
    AddConstraints(LP, lhs_le, rhs_le : Rel := "le");

    // assing constraints x_i >= -D M y_i
    lhs_ge := HorizontalJoin(one, D*M*one); 
    rhs_ge := Matrix(R, [[0] : eta in etas]);
    AddConstraints(LP, lhs_ge, rhs_ge : Rel := "ge");

    // set bounds -D M <= x_i <= D M, 0 <= y_i <= 1
    for n in [1..#etas] do
        SetLowerBound(LP, n, -R!(D*M));
        SetUpperBound(LP, n, R!(D*M));
        SetLowerBound(LP, n + #etas, R!0);
        SetUpperBound(LP, n + #etas, R!1);
    end for;

    if Type(R) ne RngInt then
        // set the y_i to be integral, so that they are either 0 or 1
        for n in [(#etas+1)..(2*#etas)] do
            SetIntegerSolutionVariables(LP, [n], true);      
        end for;
    end if;

    // Set the objective to be sum(y_i) (sparsity)
    objective := Matrix(R, 1, 2*#etas, [0 : i in [1..#etas]] cat [1 : i in [1..#etas]]);
    SetObjectiveFunction(LP, objective);

    // have LP minimize the number of non-zero entries in a solution    
    SetMaximiseFunction(LP, false);

    sol := Solution(LP);
    sol_Q := (1/D)*ChangeRing(sol[1], Rationals());
    eta := &+[sol_Q[i]*etas[i] : i in [1..#etas]];
    return eta;
end function;

intrinsic FindAsShortEtaQuotientLP(f::ModFrmElt, N::RngIntElt, k::RngIntElt) -> EtaQuot
{Uses LP to try to find a short expression of f as linear combination of eta quotients.}
    vprintf ShimuraQuotients, 1 : "Trying to find form as a linear combination of holomorphic eta quotients of weight %o and level ", k;
    eta, etas, mat, v := FindAsEtaQuotient(f, N, k);
    sol := Solution(mat, v);
    return find_short_expression_using_LP(etas, mat, v : R := Integers(), D := Denominator(sol));
end intrinsic;

import "EtaQuotient.m" : valuation_at_oo_lb;

intrinsic FindAsEtaQuotient(f::ModFrmElt, N::RngIntElt, k::RngIntElt) -> EtaQuot, SeqEnum[EtaQuots], Mtrx, Mtrx
{Returns f as an eta quotient.}
    found := false;
    M := N;
    vprintf ShimuraQuotients, 1 : "Trying to find form as a linear combination of holomorphic eta quotients of weight %o and level ", k;
    while not found do
        vprintf ShimuraQuotients, 1 : "%o ", M;
        etas := HolomorphicEtaQuotients(M, k);
        M *:= 2;
        if IsEmpty(etas) then continue; end if;
        vals := [valuation_at_oo_lb(eta) : eta in etas];
        prec := Maximum(Maximum(vals), #etas) + 2;
        mat := Matrix([AbsEltseq(qExpansionAtoo(eta,prec) : FixedLength) : eta in etas]);
        v := Vector(Integers(), AbsEltseq(qExpansion(f,prec) : FixedLength));
        mat_Q := ChangeRing(mat, Rationals());
        v_Q := ChangeRing(v, Rationals());
        found := v_Q in RowSpace(mat_Q);
    end while;
    sol := Solution(mat_Q, v_Q);
    eta := &+[sol[i]*etas[i] : i in [1..#etas]];
    return eta, etas, mat_Q, v_Q;
end intrinsic;

intrinsic FindMinimalEtaQuotient(f::ModFrmElt, N::RngIntElt, k::RngIntElt) -> EtaQuot
{Returns f as an eta quotient.}
    vprintf ShimuraQuotients, 1 : "Trying to find form as a linear combination of holomorphic eta quotients of weight %o and level ", k;
    eta, etas, mat_Q, v_Q := FindAsEtaQuotient(f, N, k);
    
    v := ChangeRing(v_Q, Integers());
    mat := ChangeRing(mat_Q, Integers());

    sol := Solution(mat_Q, v_Q);
    length := #[x : x in Eltseq(sol) | x ne 0];
    height := coeff_height(sol);
    vprintf ShimuraQuotients, 1 : "\nFound a solution with length %o and height %o\n", length, height;

    // Try to find best integral solution first
    if v in RowSpace(mat) then
        vprintf ShimuraQuotients, 1 : "Found an integral solution,"; 
        sol_Z := Solution(mat, v);
        ker := Kernel(mat);
        // The 126 < 200 < 262 was set according to performance of ClosestVector
        if (Dimension(ker) gt 0) and (Dimension(ker) lt 200) then
            vprintf ShimuraQuotients, 1 : "looking for closest vector to lattice of dimension %o...", Dimension(ker);
            sol_Z := sol_Z - ClosestVector(Lattice(ker),sol_Z);
            vprintf ShimuraQuotients, 1 : "Done!\n";
        end if;
       
        length_Z := #[x : x in Eltseq(sol_Z) | x ne 0];
        height_Z := coeff_height(sol_Z);
        if (length_Z lt length) or ((length_Z eq length) and (height_Z lt height)) then
            sol := sol_Z;
            length := length_Z;
            height := height_Z;
            vprintf ShimuraQuotients, 1 : "Improved to length %o and height %o\n", length, height;
        end if;
    end if;

    vprintf ShimuraQuotients, 1 : "Trying to improve the upper bound...";

    length_new := length;
    height_new := height;
    count_stable := 0;

    // 10000 is arbitrary here
    while (count_stable lt 10000) do
        count_stable +:= 1;
        eta_idxs := [i : i in RandomSubset({1..#etas}, length-1)];
        basis := [mat_Q[i] : i in eta_idxs];
        V := sub<RowSpace(mat_Q) | basis>;
        if v_Q notin V then continue; end if;
        sol_new := Solution(Matrix(basis), v_Q);
        length_new := #[x : x in Eltseq(sol_new) | x ne 0];
        height_new := coeff_height(sol_new);
        
        if (length_new lt length) or ((length_new eq length) and (height_new lt height)) then
            sol := sol_new;
            length := length_new;
            height := height_new;
            vprintf ShimuraQuotients, 1 : "Improved to length %o and height %o\n", length, height;
            count_stable := 0;
        end if;
    end while;

    vprintf ShimuraQuotients, 1 : "Verifying there are no expressions of length ";
    for l in [1..length-1] do
        vprintf ShimuraQuotients, 1 : "%o ", l;
        found, eta := find_minimal_expression_of_length(etas, l, mat_Q, v_Q);
        if found then return eta; end if;
    end for;
    
    if height eq 0 then
        // In this case, we know found solution is minimal
        eta := &+[sol[i]*etas[i] : i in [1..#etas]];
        return eta; 
    end if;

    vprintf ShimuraQuotients, 1 : "\nFinding minimal expression of length %o...\n", length;
    // Now this must work
    found, eta := find_minimal_expression_of_length(etas, length, mat_Q, v_Q);
    assert found;
    return eta;

end intrinsic;

intrinsic WeaklyHolomorphicBasis(D::RngIntElt,N::RngIntElt : Prec := 100, Zero := false, n0 := 0) -> .
{Returns a weakly holomorphic basis corresponding to D, N.}
    D0,M,g := get_D0_M_g(D,N);
    Ldata := ShimuraCurveLattice(D,N);
    L := Ldata`L;
    Ldual := Ldata`Ldual;
    disc := #(Ldual/L);
    rk := -1;
    dim := 0;
    t, lhs, rhs, n_eq, n_ds := find_t(M);
    t := Eltseq(t);
    n_gaps := Zero select 0 else g - &+[d div 4 : d in Divisors(D0)];
    k := t[#t]; // the order of pole for t
    n := n_gaps;
    if not Zero then
        // Trying n0 here
        n0 := n_gaps;
    end if;
    // The n0 below is guaranteed to work (Lemma 17 and Lemma 27 in [GY]), but might be too large
    // n0 := Maximum(2*g-2 - &+[d div 4 : d in Divisors(D0)] - (Zero select 0 else k), 0);
    n := n0;
    R := EtaQuotientsRing(M, disc);
    gap_condition := false;
    pole_string := Zero select "{0,oo}" else "{oo}";
    vprintf ShimuraQuotients, 2: "\n\tComputing generators for the ring of %o-weakly holomorphic modular forms of level %o...", pole_string, M;
    while (rk lt dim) or (not gap_condition) do
        
        vprintf ShimuraQuotients, 3: "\n\t\tprec = %o, n = %o, k = %o, rk = %o, dim = %o...", Prec, n, k, rk, dim; 
        
        rs := get_integer_prog_solutions(M, lhs, rhs, n_eq, n_ds, n + (Zero select 0 else k), Zero select k else 0);
        
        eta_quotients := [EtaQuotient(R,r) : r in rs];
        t_eta_quotient := EtaQuotient(R,t[1..n_ds]);

        if Zero then
            eta_quotients_oo := eta_quotients;
            
            eta_quotients := [SAction(eta) : eta in eta_quotients];
            t_eta_quotient := SAction(t_eta_quotient : Admissible := false);
        end if;

        qexps := [qExpansionAtoo(eta, Prec) : eta in eta_quotients];
        _<q> := Universe(qexps);
        min_v := Minimum([Valuation(f) : f in qexps]);
        coeffs := Matrix(Rationals(), [AbsEltseq(q^(-min_v)*f : FixedLength) : f in qexps]);
        
        E, T := EchelonForm(coeffs);
        E := Submatrix(E, [1..Rank(E)], [1..Ncols(E)]);
        if Zero then
            eta_quotients_oo := [&+[T[i][j]*eta_quotients_oo[j] : j in [1..#eta_quotients_oo]] : i in [1..Nrows(E)] ];
        end if;
        dim := n + k + &+[d div 4 : d in Divisors(D0)] + 1 - g;
        rk := Rank(E);
        pole_orders := [PivotColumn(E,i) + min_v - 1 : i in [1..Rank(E)]];
        gaps := &cat[[pole_orders[i]+1..pole_orders[i+1]-1] : i in [1..Minimum(#pole_orders-1,1-min_v)] 
                                                            | pole_orders[i+1] - pole_orders[i] gt 1];
        max_pole := (#gaps eq 0) select 0 else -gaps[1];
        
        if (rk lt dim) then
            if (dim gt Prec) then
                Prec := dim;
            else
                Prec +:= k;
                // update n
                n +:= k;
            end if;
        else
            n0 := (n0 le max_pole) select max_pole else n + max_pole;
            gap_condition := (#gaps eq n_gaps) and (n ge max_pole);
            if (#gaps ne n_gaps) then
                n := n0;
            end if;
            if (n lt max_pole) then
                n := max_pole + k;
            end if;
        end if;
    end while;
    vprintf ShimuraQuotients, 3 : "\n";
    vprintf ShimuraQuotients, 2 : "Done!";
    // sanity checks
    assert rk eq dim;
    
    n := -min_v;
    if Zero then
        E := Submatrix(E, [1..n], [1..Ncols(E)]);
    else
        n0 := max_pole;
    end if;
    
    eta_quotients := [&+[T[i][j]*eta_quotients[j] : j in [1..#eta_quotients]] : i in [1..Nrows(E)] ];
    
    if Zero then 
        return E, n, n0, eta_quotients_oo, eta_quotients; 
    end if;
    return E, n, n0, t_eta_quotient, eta_quotients; 
end intrinsic;

function fourth_power_free(a)
    ps := PrimeDivisors(Numerator(a)) cat PrimeDivisors(Denominator(a));
    vals := [Valuation(a,p) : p in ps];
    v_rad := [v div 4 : v in vals];
    v_free := [v mod 4 : v in vals];
    rad := &*[Rationals() | p^v_rad[i] : i->p in ps];
    free := &*[Rationals() | p^v_free[i] : i->p in ps];
    assert rad^4 * free eq a;
    return free, rad;
end function;

// Returns the q-expansion of eta|[g] where g is a 2x2 matrix
// with integral entries
// eta|[g] = eta(gz)/phi_g(z),
// where phi_g is the multiplier system such that phi_S(z) = sqrt(-iz),
// phi_T(z) = zeta_{24}, phi_{a*I}(z) = 1 and phi_{Diagonal(d,1)} = d^{-1/4}
// returns c*eta|[g], c^4 for some constant c.
function eta_action(g)
    B, T := HermiteForm(g);
    assert g eq T^(-1)*B; // so eta|[g] = eta|[B]
    assert B[2,1] eq 0; // making sure this matrix is Borel
    a := B[1,1];
    b := B[1,2];
    d := B[2,2];
    n := d div GCD(b,d); // Order of root of unity
    N := LCM(n, 24);
    K<zeta> := CyclotomicField(N);
    zeta_24 := zeta^(N div 24);
    zeta_n := zeta^(N div n);
    _<q> := PuiseuxSeriesRing(K);
    eta<q> := DedekindEta(q);
    // We use the fact that if q(z) = e^{2 pi i z} then 
    // q((az+b)/d) = q^{a/d} zeta_d^b
    nor_eta := eta / q^(1/24);
    etaB := zeta_n^(b div GCD(b,d)) * q^(a/(24*d)) * Evaluate(nor_eta, q^(a/d));
    // etaB is eta(Bz)
    // we need to divide by the multiplier which is (d/a)^(1/4) zeta_24^b
    free, rad := fourth_power_free(d/a);
    return rad^(-1)*zeta_24^(-b)*etaB, free;
end function;

// Returns the q-expansion of f|[g],
// where f is the product of eta(dz)^(rz) where (d,r)
// run over ds and rs repsectively
function eta_quotient_action(rs, ds, g)
    ret_ceta := 1;
    ret_c4 := 1;
    for i->d in ds do
        alpha_d := DiagonalMatrix([d,1]); // eta(dz) = d^(-1/4) eta|[alpha_d]
        c_eta, c4 := eta_action(alpha_d*g);
        c_eta := c_eta^rs[i];
        c4 := (c4*d)^rs[i];
        ret_ceta *:= c_eta;
        ret_c4 *:= c4;
    end for;
    /*
    is_fourth, c := IsPower(ret_c4, 4);
    assert is_fourth; // For now we hope this is true for our functions
    */
    free, rad := fourth_power_free(ret_c4);
    return rad^(-1)*ret_ceta, free;
    // return ret_ceta / c;
end function;

// Returns the q-expansion of f|[g],
// for the functions f described by alphas.
// alphas describes a linear combination
// of the eta quotients given by rs
function linear_comb_eta_quotients_action(alphas, rs, ds, g)
    L := Rationals();
    ret_eta := PuiseuxSeriesRing(L)!0;
    for j->r in rs do
        ceta, c4 := eta_quotient_action(r, ds, g);
        ceta := alphas[j]*ceta;
        Feta := BaseRing(Parent(ceta));
        Fetax<x> := PolynomialRing(Feta);
        fac := Factorization(x^4 - c4);
        min_deg, min_place := Minimum([Degree(fa[1]) : fa in fac]);
        if min_deg eq 1 then 
            Leta := Feta;
            c := -Coefficient(fac[min_place][1],0);
            assert c^4 eq c4;
        else
            Leta<c> := ext<Feta | fac[min_place][1]>;
        end if;
        eta := ceta/c;
        L := CompositeFields(L, Leta)[1];
        eta := ChangeRing(Parent(eta),L)!eta;
        ret_eta := ChangeRing(Parent(ret_eta),L)!ret_eta + eta;
    end for;
    // return &+[alphas[j]*eta_quotient_action(r, ds, g) : j->r in rs];
    // free, rad := fourth_power_free(ret_c4);
    // return rad^(-1)*ret_ceta, free;
    return ret_eta;
end function;

/*
// assuming v_i is the coefficient of eta_i in Ldual / L
function WeilRepresentation(gamma, v, Ldual, discL, Qdisc, to_disc)
    PSL2Z := PSL2(Integers());
    gens := Generators(PSL2Z);
    T := gens[1];
    S := gens[2];
    assert Eltseq(T) eq [1,1,0,1];
    assert Eltseq(S) eq [0,1,-1,0];
    word := FindWord(PSL2Z, gamma);
    // assert &*[gens[Abs(word[i])]^Sign(word[i]) : i in [1..#word]] eq gamma;
    w := v;
    deltas := [delta : delta in discL];
    B := BasisMatrix(Ldual); // represent in the basis of Ldual, for which the gram matrix is Qdisc
    delta_lifts := [ChangeRing(Solution(B,delta@@to_disc), Rationals()) : delta in deltas];
    norms := [-(delta*Qdisc, delta) / 2 : delta in delta_lifts]; // -<delta,delta>/2
    cycl_order := LCM([Denominator(n) : n in norms] cat [8]);
    K<zeta> := CyclotomicField(cycl_order);
    is_sqr, sqrt_disc := IsSquare(K!(#discL));
    assert is_sqr;
    for i in [1..#word] do
        if Abs(word[i]) eq 1 then
            // rho_L(T)
            w := [zeta^(Sign(word[i])*Integers()!(norms[j]*cycl_order)) * x : j->x in w];
        else // Abs(word[i]) eq 2 
            scalar := zeta^(-cycl_order div 8) / sqrt_disc;
            w := [ scalar*&+[zeta^(Integers()!(cycl_order*(delta_lifts[k]*Qdisc, delta_lifts[j]))) * w[k] : k in [1..#w]] : j->x in w];
        end if;
    end for;
    return w;
end function;

// We would like to be able to start with a f,
// given as a linear combination of Eta quotients
// and create the coefficients of the Borcherds form F
// obtained from it. 
// !! TODO : Not working yet - the constants are off
function SpreadBorcherds(alphas, rs, ds, Ldual, discL, Qdisc, to_disc, M)

function SpreadBorcherds(alphas, rs, ds, Ldual, discL, Qdisc, to_disc)
    M := #discL;
    G := PSL2(Integers());
    gens := Generators(G);
    gammas := CosetRepresentatives(Gamma0(M));
    T := gens[1];
    S := gens[2];
    assert Eltseq(T) eq [1,1,0,1];
    assert Eltseq(S) eq [0,1,-1,0];
    F := [* PuiseuxSeriesRing(Rationals())!0 : i in [1..#discL] *];
    e0 := [1] cat [0 : i in [1..#discL-1]];
    for gamma in gammas do
        // compute f|gamma rho_L(gamma^(-1)) e_0 and add to the sum
        cf_gamma := linear_comb_eta_quotients_action(alphas, rs, ds, Matrix(gamma));
        rho_L := WeilRepresentation(gamma^(-1), e0, Ldual, discL, Qdisc, to_disc);
        for i in [1..M] do
            coord := rho_L[i]*cf_gamma;
            L := BaseRing(Parent(coord));
            F[i] := ChangeRing(Parent(F[i]),L)!F[i] + coord;
        end for;
    end for;
    return F;
end function;*/

intrinsic FindLambda(Q::AlgMatElt, d::RngIntElt, Order::AlgQuatOrd, basis_L::SeqEnum : bound := 10)-> BoolElt, ModTupRngElt
{.}
    require d gt 0: "d must be positive";

    Q := ChangeRing(Q, Integers());
    n := Nrows(Q);
    idxs := CartesianPower([-bound..bound], n);
    for idx in idxs do
        v := Vector([idx[j] : j in [1..n]]);
        v := ChangeRing(v, BaseRing(Q));
        if (v*Q,v) eq 2*d then
            // checking whether this is an optimal embedding of the order of discriminant d
            elt := &+[v[i]*basis_L[i] : i in [1..#basis_L]];
            if d mod 4 ne 3 then
                assert d mod 4 eq 0;
                if elt/2 in Order then
                    return true, v;
                end if;
            end if;
            // d mod 4 eq 3
            if (1+elt)/2 in Order then
                return true, v;
            end if;
        end if;
    end for;
    return false, _;
end intrinsic;

intrinsic FindLambdas(Q::AlgMatElt, ds::SeqEnum[RngIntElt], Order::AlgQuatOrd, basis_L::SeqEnum : bound := 10, lambda_array := AssociativeArray())-> BoolElt, ModTupRngElt
{.}
    require &and[d gt 0 : d in ds]: "All ds must be positive";
    lambdas := lambda_array;
    Q := ChangeRing(Q, Integers());
    n := Nrows(Q);
    idxs := CartesianPower([-bound..bound], n);
    twice_ds := [2*d : d in ds];
    for idx in idxs do
        v := Vector([idx[j] : j in [1..n]]);
        v := ChangeRing(v, BaseRing(Q));
        if (v*Q,v) in twice_ds then
            d := (v*Q,v) div 2;
            if d in Keys(lambdas) then continue; end if; //already found
            // checking whether this is an optimal embedding of the order of discriminant d
            if d mod 4 ne 3 then
                assert d mod 4 eq 0;
                elt := &+[v[i]/2*basis_L[i] : i in [1..#basis_L]]; //checking elt/2 embedding
            // d mod 4 eq 3
            elif d mod 4 eq 3 then
            //(1+elt)/2 now
                elt := &+[v[i]/2*basis_L[i] : i in [1..#basis_L]]+Parent(basis_L[1])!1/2;
            end if;
            //check that this gives an optimal embedding
            B_O := Basis(Order);
            Mat_O := Matrix([Eltseq(B_O[i]) : i in [1..#B_O]]);
            CoB := Matrix([Eltseq(Solution(Mat_O, Vector(Rationals(),[1,0,0,0]))), Eltseq(Solution(Mat_O, Vector(Rationals(), Eltseq(elt))))]);
            den := Denominator(CoB);
            CoBZZ := ChangeRing(den*CoB, Integers());
            S, _, _ := SmithForm(CoBZZ);
            Sprime := HorizontalJoin(IdentityMatrix(Integers(),2), ZeroMatrix(Integers(),2));
            if S eq Sprime then
                lambdas[d] := v;
                if Keys(lambdas) eq Set(ds) then
                    //then we found all the lambdas
                        return true, lambdas;
                end if;
            end if;
        end if;
    end for;
    if Keys(lambdas) eq Set(ds) then
        return true, lambdas;
    else
        vprint ShimuraQuotients, 2 : "Could not find all lambdas";
        vprint ShimuraQuotients, 2 : Set(ds) diff Keys(lambdas);
        return false, lambdas; //return the partial progress
    end if;
end intrinsic;

intrinsic ElementOfNorm(Q::AlgMatElt, d::RngIntElt, Order::AlgQuatOrd, basis_L::SeqEnum) -> ModTupRngElt
{Return element of norm d in the quadratic space with Gram matrix Q.
Warning - does it in a silly way via enumeration. }
    require d gt 0: "d must be a positive norm";
    bd := 10;
    found_lambda := false;
    while not found_lambda do
        bd *:= 2;
        found_lambda, lambda := FindLambda(Q, d, Order, basis_L : bound := bd);
    end while;
    assert found_lambda;
    return lambda;
end intrinsic;

intrinsic ElementsOfNorm(Q::AlgMatElt, ds::SeqEnum[RngIntElt], Order::AlgQuatOrd, basis_L::SeqEnum) -> ModTupRngElt
{Return elements of norm in the list of ds in the quadratic space with Gram matrix Q.
Warning - does it in a silly way via enumeration. }
    require &and[d gt 0 : d in ds]: "All ds must be positive";
    max_d := Maximum(ds);
    bd := max_d div 2;
    lambdas := AssociativeArray();
    found_lambdas := false;
    vprintf ShimuraQuotients, 2 : "\n\tFinding lambdas for norms in %o...", ds;
    while not found_lambdas do
        found_lambdas, lambdas := FindLambdas(Q, ds, Order, basis_L : bound := bd, lambda_array := lambdas);
        if not found_lambdas then
            bd *:=2;
            vprintf ShimuraQuotients, 2 : "Increasing lambda bound to %o\n", bd;
        end if;
    end while;
    vprintf ShimuraQuotients, 2 : "Found lambdas.";
    assert found_lambdas;
    return lambdas;
end intrinsic;

intrinsic VerticalJoinList(mats::List)->.
    {}
    m := mats[1];
    for i in [2..#mats] do
        m := VerticalJoin(m, mats[i]);
    end for;
    return m;
end intrinsic;

function my_legendre_symbol(alpha, p)
    return LegendreSymbol(Integers()!(GF(p)!alpha),p);
end function;

// W_{m,p}
// L should be Lminus
function Wpoly(m,p,mu,L,K,Q)
    _<sqrtp> := K;
    F := QNF();
    S := BasisMatrix(L)*Q*Transpose(BasisMatrix(L));
    n := Rank(L);
    Lnf := NumberFieldLatticeWithGram(ChangeRing(S,F));
    // l is the sequence of exponents
    assert p ne 2; // take care of p = 2 later
    bases, Jblocks, exps := JordanDecomposition(Lnf,p*Integers(F)); 
    l := &cat[[e : j in [1..Nrows(Jblocks[i])]] : i->e in exps];
    eps := &cat[[Rationals() | 1/2 * x / p^exps[i] : x in Diagonal(b)]  : i->b in Jblocks]; // so that S is equivalent to (2 eps_1 p^{l_1},..., 2 eps_n p^{l_n})
    assert &and[Valuation(e,p) eq 0 : e in eps];
    B := ChangeRing(VerticalJoinList(bases), Rationals());
    mu_wrt_L := Solution(ChangeRing(BasisMatrix(L),Rationals()), ChangeRing(mu, Rationals()));
    Q_mu := 1/2*(mu_wrt_L * ChangeRing(S,Rationals()), mu_wrt_L);
    R<x> := PolynomialRing(K);
    if not IsIntegral(m - Q_mu) then
        return R!0;
    end if;
    mu_wrt_B := mu_wrt_L*B^(-1);
    H_mu := {i : i in [1..n] | Valuation(mu_wrt_B[i],p) ge 0};
    vals := [l[i] + Valuation(mu_wrt_B[i], p) : i in [1..n] | i notin H_mu];
    if IsEmpty(vals) then
        K_0 := Infinity();
    else
        K_0 := Minimum(vals);
    end if;
    L_mu := func<k | {i : i in H_mu | IsOdd(l[i] - k) and (l[i] - k lt 0)}>;
    l_mu := func<k | #L_mu(k)>;
    // we compute twice d_mu for technical reasons
    d2_mu := func<k | 2*k + &+[Minimum(l[i]-k, 0) : i in H_mu]>;
    eps_mu := func<k | LegendreSymbol(-1,p)^(l_mu(k) div 2) * &*[Integers() | my_legendre_symbol(eps[i],p) : i in L_mu(k)]>;
    f_1 := function(x)
        a, alpha := Valuation(x,p);
        return IsEven(l_mu(a+1)) select -1/p else my_legendre_symbol(alpha,p) / sqrtp;
    end function;
    t_mu := m - &+[Rationals() | eps[i]*p^l[i]*mu_wrt_B[i]^2 : i in [1..n] | i notin H_mu];
    a := Valuation(t_mu, p);
    
    if a lt K_0 then
        ret := 1;
        ret +:= (1 - 1/p)*&+[R | eps_mu(k)*sqrtp^d2_mu(k)*x^k : k in [1..a] | IsEven(l_mu(k))];
        ret +:= eps_mu(a+1)*f_1(t_mu)*sqrtp^d2_mu(a+1)*x^(a+1);
    else
        ret := 1;
        ret +:= (1 - 1/p)*&+[R | eps_mu(k)*sqrtp^d2_mu(k)*x^k : k in [1..K_0] | IsEven(l_mu(k))];
    end if;
    return ret;
end function;

// This is for p = 2
function Wpoly2(m,mu,L,K,Q)
    p := 2;
    Zp := pAdicRing(p);
    _<sqrtp> := K;
    F := QNF();
    S := BasisMatrix(L)*Q*Transpose(BasisMatrix(L));
    n := Rank(L);
    Lnf := NumberFieldLatticeWithGram(ChangeRing(S,F));
    // l is the sequence of exponents
    assert p eq 2;
    bases, Jblocks, exps := JordanDecomposition(Lnf,p*Integers(F)); 
    bases := [* ChangeRing(B, Rationals()) : B in bases *];
    Jblocks := [* ChangeRing(J, Rationals()) : J in Jblocks *];
    bases := [* ChangeRing(B, Zp) : B in bases *];
    Jblocks := [* ChangeRing(J, Zp) : J in Jblocks *];
    l_list := [];
    m_list := [];
    n_list := [];
    eps := [];
    mu_indices := [];
    // For these we record the first index of two, so that (i, i+1) are the indices
    mu_prime_indices := [];
    mu_prime_prime_indices := [];
    // It seems like eps_prime and eps_prime_prime can always be taken to be 1
    // eps_prime := [];
    // eps_prime_prime := [];
    row_ind := 0;
    for i->Jblock in Jblocks do
        if Nrows(Jblock) eq 1 then
            row_ind +:= 1;
            Append(~l_list, exps[i]);
            Append(~eps, Jblock[1,1] / p^exps[i]);
            Append(~mu_indices, row_ind);
        end if;
        for j in [2..Nrows(Jblock)] do
            row_ind +:= 1;
            b := Jblock[j-1,j] / p^(exps[i]);
            if IsWeaklyZero(b) then
                Append(~l_list, exps[i]);
                Append(~eps, Jblock[j-1,j-1] / p^exps[i]);
                Append(~mu_indices, row_ind);
                if (j eq Nrows(Jblock)) then
                    row_ind +:= 1;
                    Append(~l_list, exps[i]);
                    Append(~eps, Jblock[j,j] / p^exps[i]);
                    Append(~mu_indices, row_ind);
                end if;
                continue; 
            end if;
            a := Jblock[j-1,j-1] / p^(exps[i]);
            d := Jblock[j,j] / p^(exps[i]);
            disc := b^2 - a*d;
            if (Integers(8)!(Integers()!disc) eq 5) then 
                disc +:= 2*a; 
                Append(~n_list, exps[i]);
                Append(~mu_prime_prime_indices, row_ind);
                aniso := true;
            else   
                Append(~m_list, exps[i]);
                Append(~mu_prime_indices, row_ind);
                aniso := false;
            end if;
            is_sqr, sqrt_disc := IsSquare(disc);
            assert is_sqr;
            if IsWeaklyZero(a) then
                B := Matrix(Zp, [[1,0],[0,1]]);
            else
                // solving the quadratic
                x1 := (-b + sqrt_disc) div a;
                x2 := (-b - sqrt_disc) div a;
                if aniso then
                    inner_product := 2-(2*disc div a);
                    x := Zp!inner_product;
                    beta := Sqrt((4-x^2)/3);
                    alpha := (-beta + x)/2;
                    // if v_1, v2 are a basis for [[2,1],[1,2]], then v_1, alpha v1 + beta v2 are a basis for [[2,x],[x,2]]
                    get_to_x := Matrix([[1,0],[alpha, beta]]);
                    assert &and[IsWeaklyZero(e) : e in Eltseq(get_to_x * Matrix([[2,1],[1,2]]) * Transpose(get_to_x) - Matrix([[2,x],[x,2]]))];
                    B := Matrix([[x1, 1], [x2, 1]]);
                    B := get_to_x^(-1)*B;
                else
                    z2 := (-a)/(2*disc); // constant to get scalar product equal to 1 
                    // Change of basis matrix
                    B := Matrix([[x1, 1], [z2*x2, z2]]);
                end if;
            end if;
            cans := [SymmetricMatrix([0,1,0]), SymmetricMatrix([2,1,2])];
            assert B*Matrix([[a,b],[b,d]])*Transpose(B) in cans;
            // !!!! This stopped working in V2.29 !!!!
            // B_big := Parent(Jblock)!1;
            B_big := IdentityMatrix(BaseRing(Jblock),Nrows(Jblock));
            B_big[j-1,j-1] := B[1,1];
            B_big[j-1,j] := B[1,2];
            B_big[j,j-1] := B[2,1];
            B_big[j,j] := B[2,2];
            bases[i] := B_big * bases[i];
            Jblocks[i] := B_big * Jblocks[i] * Transpose(B_big);
            row_ind +:= 1;
        end for;
    end for;
    eps_prime := [1 : m in m_list];
    eps_prime_prime := [1 : n in n_list];
    H := #l_list;
    M := #m_list;
    N := #n_list;
    assert n eq H + 2*M + 2*N;
    assert &and[Valuation(e) eq 0 : e in eps];
    B := ChangeRing(VerticalJoinList(bases), Rationals());
    mu_wrt_L := Solution(ChangeRing(BasisMatrix(L),Rationals()), ChangeRing(mu, Rationals()));
    Q_mu := 1/2*(mu_wrt_L * ChangeRing(S,Rationals()), mu_wrt_L);
    R<x> := PolynomialRing(K);
    if not IsIntegral(m - Q_mu) then
        return R!0;
    end if;
    mu_wrt_B := mu_wrt_L*B^(-1);
    mu_list := [mu_wrt_B[i] : i in mu_indices];
    mu_prime_list := [[mu_wrt_B[i], mu_wrt_B[i+1]] : i in mu_prime_indices];
    mu_prime_prime_list := [[mu_wrt_B[i], mu_wrt_B[i+1]] : i in mu_prime_prime_indices];
    M_mu := {i : i in [1..M] | Valuation(mu_prime_list[i][1],p) ge 0 and Valuation(mu_prime_list[i][2],p) ge 0};
    N_mu := {i : i in [1..N] | Valuation(mu_prime_prime_list[i][1],p) ge 0 and Valuation(mu_prime_prime_list[i][2],p) ge 0};
    H_mu := {i : i in [1..H] | Valuation(mu_list[i],p) ge 0};

    L_mu := func<k | {i : i in H_mu | IsOdd(l_list[i] - k) and (l_list[i] - k lt 0)}>;
    l_mu := func<k | #L_mu(k)>;
    // we compute twice d_mu for technical reasons
    d2_mu := func<k | 2*k + &+[Integers()|Minimum(l_list[i]-k, 0) : i in H_mu] + 2*&+[Integers()|Minimum(m_list[i]-k,0) : i in M_mu] + 2*&+[Integers()|Minimum(n_list[i]-k,0) : i in N_mu]>;

    vals := [l_list[i] + Valuation(mu_list[i], p) : i in [1..H] | i notin H_mu and Valuation(mu_list[i],p) lt -1];
    vals cat:= [l_list[i] : i in [1..H] | i notin H_mu and Valuation(mu_list[i],p) eq -1];
    vals cat:= [m_list[i] + Minimum(Valuation(mu_prime_list[i][1], p), Valuation(mu_prime_list[i][2], p)) : i in [1..M] | i notin M_mu];
    vals cat:= [n_list[i] + Minimum(Valuation(mu_prime_prime_list[i][1], p), Valuation(mu_prime_prime_list[i][2], p)) : i in [1..N] | i notin N_mu];
    if IsEmpty(vals) then
        K_0 := Infinity();
    else
        K_0 := Minimum(vals);
    end if;
    
    p_mu := func<k | (-1)^(&+[Integers()|Minimum(n_list[j] - k, 0) : j in N_mu])>;
    eps_mu := func<k | &*[Zp|eps[h] : h in L_mu(k)]>;

    delta_mu := func<k | exists(h){h : h in H_mu | l_list[h] eq k} select 0 else 1>;

    two_block := func< x | x[1]^2 + x[1]*x[2] + x[2]^2>;

    Q_prime_mu := &+[Rationals() | eps[i]*p^(l_list[i]-1)*mu_list[i]^2 : i in [1..H] | i notin H_mu];
    Q_prime_mu +:= &+[Rationals() | eps_prime[i]*p^m_list[i]*(&* mu_prime_list[i]) : i in [1..M] | i notin M_mu];
    Q_prime_mu +:= &+[Rationals() | eps_prime_prime[i]*p^n_list[i]*two_block(mu_prime_prime_list[i]) : i in [1..N] | i notin N_mu];

    t_mu := m - Q_prime_mu;
    a := Valuation(t_mu, p);

    nu := func< k | t_mu*p^(3-k) - &+[Zp|eps[h] : h in H_mu | l_list[h] lt k]>;

    psi_char := func<k | (Valuation(nu(k)) ge 2) select (-1)^(Integers()!(GF(2))!(nu(k)/ 4)) else 0>;

    K_0 := Minimum(K_0, a+3);

    
    ret := 1;
    ret +:= &+[R | delta_mu(k)*p_mu(k)*sqrtp^(d2_mu(k)-3)*KroneckerSymbol(2,Integers()!(eps_mu(k)*nu(k)))*x^k : k in [1..K_0] | IsOdd(l_mu(k))];
    ret +:= &+[R | delta_mu(k)*p_mu(k)*sqrtp^(d2_mu(k)-2)*KroneckerSymbol(2,Integers()!(eps_mu(k)))*psi_char(k)*x^k : k in [1..K_0] | IsEven(l_mu(k))];

    return ret;
end function;

function Wpoly_scaled(m,p,mu,L,Q : scaled := true)
    S := BasisMatrix(L)*Q*Transpose(BasisMatrix(L));
    D := Determinant(S);
    vpD := Valuation(D,p);
    K<sqrtp> := QuadraticField(p);
    scale := sqrtp^(-vpD);
    euler := p eq 2 select Wpoly2(m,mu,L,K,Q) else Wpoly(m,p,mu,L,K,Q);
    assert CanChangeRing(euler, Rationals());
    return (scaled) select scale*euler else ChangeRing(euler, Rationals()), p^(-vpD);
end function;

function W(m,p,mu,L,Q)
    Wpoly := Wpoly_scaled(m,p,mu,L,Q);
    _<sqrtp> := BaseRing(Wpoly);
    n := Rank(L);
    // s0 := n/2 - 1;
    // s := -s0;
    s2 := 2 - n; // s2 = 2*s always integral
    return Evaluate(Wpoly, sqrtp^(-s2));
end function;

function get_Wpolys(m,mu,Lminus,Q, Sm_mu : scaled := true)
    if scaled then 
        Wpolys := [* Wpoly_scaled(m,p,mu,Lminus,Q) : p in Sm_mu *];
        wpolyseval := [* Evaluate(Wpolys[i],1) : i in [1..#Wpolys] *];
        return Wpolys, wpolyseval;
    end if;

    Wpolys := [];
    scales_sqr := [Rationals() | ];
    wpolyseval := [];
    for p in Sm_mu do
        Wpol, scale_sqr :=  Wpoly_scaled(m,p,mu,Lminus,Q : scaled := false);
        Append(~Wpolys, Wpol);
        Append(~scales_sqr, scale_sqr);
        Append(~wpolyseval, Evaluate(Wpol, 1));
    end for;
    
    return Wpolys, wpolyseval, scales_sqr;
end function;

function find_v(F, places, sqrtps)
    // assert exists(v){v : v in places | &and[Evaluate(F!sqrtp, v) gt 0 : sqrtp in sqrtps]};
    pos := false;
    for v in places do
        pos := true;
        for sqrtp in sqrtps do
            if Evaluate(F!sqrtp, v) lt 0 then pos := false; break; end if;
        end for;
        if pos then v_pos := v; break; end if;
    end for;
    assert pos;
    return v_pos;
end function;

function get_kappa_minus_squared(d, Wpolys, Wpol, Sm_mu, i, scales_sqr)
    W_prod := &*[ Rationals() | Evaluate(Wpolys[j],1) : j in [1..#Sm_mu] | j ne i];
    W_prod *:= Evaluate(Derivative(Wpol),1); // this should be multiplied by -log(p_prime)

    kron_prod := &*[Rationals() | 1 - Evaluate(KroneckerCharacter(d),p)/p : p in Sm_mu];

    h := ClassNumber(d);
    w := #UnitGroup(QuadraticField(d));

    scale_sqr := &*scales_sqr;
    
    W_kron := W_prod / kron_prod;
    // km_sqr := -d*scale_sqr*(w*W_kron / h)^2;
    // Using Yang's code to try to work the non-maximal case
    // !!! Not sure why this works !!!
    _, f := SquarefreeFactorization(d div FundamentalDiscriminant(d));
    km_sqr := -d*scale_sqr*(w*W_kron / h)^2 / f^2;
    km_sign := -Sign(W_kron);

    return km_sqr, km_sign;
   
end function;

// returns x,y such that the answer is x logy
function kappaminus(mu, m, Lminus, Q, d)
    // error if m eq 0, "Not implemented for m eq 0 at CM point!\n", d;  
    if m eq 0 then
        printf "\nWarning: Adding 0 for m eq 0 at CM point!\n";
        return 0, 1;
    end if;
    Bminus := BasisMatrix(Lminus);
    Delta := Determinant(Bminus*Q*Transpose(Bminus));
    
    Sm_mu := {p : p in PrimeDivisors(Delta)} join {p : p in PrimeDivisors(Numerator(m))};
    Sm_mu := [p : p in Sm_mu];
    
    vprintf ShimuraQuotients, 5: "\t\t\tSm_mu = %o\n", Sm_mu;

    Wpolys, wpolyseval, scales_sqr := get_Wpolys(m,mu,Lminus,Q, Sm_mu : scaled := false);
   
    assert exists(i){i : i in [1..#Sm_mu] | wpolyseval[i] eq 0};
    p_prime := Sm_mu[i];
    if exists(j){j : j in [1..#Sm_mu] | wpolyseval[j] eq 0 and j ne i} then
        return 0, p_prime;
    end if;
    Wpol := Wpolys[i];

    // F, sqrtps := get_field(Wpolys);

    // ret := get_kappa_minus(F, d, Wpolys, Wpol, Sm_mu, i, sqrtps);
    // ret_squared, ret_sign := get_kappa_minus_squared(F, d, Wpolys, Wpol, Sm_mu, i, sqrtps);
    ret_squared, ret_sign := get_kappa_minus_squared(d, Wpolys, Wpol, Sm_mu, i, scales_sqr);
    // assert ret_squared eq ret^2 and ret_sign eq Sign(ret);
    is_sqr, ret := IsSquare(ret_squared);
    assert is_sqr;
    ret := ret_sign*AbsoluteValue(ret);

    vprintf ShimuraQuotients, 5 : "\t\t\tadding %o log %o\n", -ret, p_prime;
    return -ret, p_prime; // to get x logy instead of -xlogy
    // return p_prime^(-ret);
end function;

function kappaminuszero(D,N,d)
    log_coeffs := AssociativeArray();
    for p in PrimeDivisors(D div GCD(d,D)) do
        log_coeffs[p] := (p-1)/(p+1);
    end for;
    for p in PrimeDivisors(N div GCD(d,N)) do
        log_coeffs[p] := 1;
    end for;
    // !! TODO := think about precision
    RR := RealField();
    pi := Pi(RR);
    gamma :=  EulerGamma(RR);
    chi := KroneckerCharacter(d);
    mu := #UnitGroup(QuadraticField(d));
    h := ClassNumber(d);
    chowla_selberg := &+[chi(a)*mu*Log(Gamma(a/AbsoluteValue(d)))/h : a in [1..AbsoluteValue(d)-1]];
    chowla_selberg +:= Log(4*pi) - 3*Log(AbsoluteValue(d)) + gamma;
    return log_coeffs, chowla_selberg;
end function;

// Computes kappa0(m) in Schofer's formula
intrinsic Kappa0(m::RngIntElt, d::RngIntElt, Q::AlgMatElt, lambda_v::ModTupRngElt) -> LogSm
{Computing coefficients Kappa0(m) in Schofers formula}
    return Kappa(Parent(lambda_v)!0,Rationals()!m,d,Q,lambda_v);
end intrinsic;

intrinsic Kappa(gamma::ModTupRngElt, m::FldRatElt, d::RngIntElt, Q::AlgMatElt, lambda_v::ModTupRngElt) -> LogSm
{Computing coefficients Kappa(gamma, m) in Schofers formula.}
    vprintf ShimuraQuotients, 4: "\n\t\t\tKappa_%o of %o", gamma, m;
    Qrat := ChangeRing(Q, Rationals());
    Q := ChangeRing(Q, Integers());
    
    c_Lplus := Content(lambda_v);
    Lplus := RSpaceWithBasis(Matrix(lambda_v div c_Lplus));
    Lminus := Kernel(Transpose(Matrix(lambda_v*Q)));
    L := RSpaceWithBasis(IdentityMatrix(Integers(),3));
    L_quo, L_quo_map := L / (Lplus + Lminus);

    lambda_rat := ChangeRing(lambda_v, Rationals());
    gamma_rat := ChangeRing(gamma, Rationals());
    c_gamma_plus := ((gamma_rat*Qrat, lambda_rat)/(lambda_rat*Qrat,lambda_rat));
    gamma_plus:= c_gamma_plus*lambda_rat;
    gamma_minus := gamma_rat - gamma_plus;
    log_coeffs := LogSum();
    // This is the condition from Yang code, if we have a vector with Q(x) = m
    Yang_tt := false;
    for mu_bar in L_quo do
        mu := mu_bar@@L_quo_map;
        c_mu_plus := ((mu*Q, lambda_v)/(lambda_v*Q,lambda_v));
        mu_plus:= c_mu_plus*ChangeRing(lambda_v, Rationals());
        mu_minus := mu - mu_plus;
        // finding the possible range of x in mu_plus + L_plus
        // use that mu_plus = c_mu_plus * lambda, L_plus = Z * c_Lplus^(-1) * lambda
        // gamma_plus = c_gamma_plus*lambda
        // that <lambda,lambda> = -2d, and we only need x with <x,x> <= 2m
        // so if x = c_gamma_plus + c_mu_plus + c_Lplus^(-1)*k, we need only those with
        // (c_gamma_plus + c_mu_plus + c_Lplus^(-1)*k)^2 le m/(-d)
        // thus k is between the following bounds
        sqr_bd := m/(-d);
        
        lb := Ceiling((-Sqrt(sqr_bd) - c_mu_plus - c_gamma_plus)*c_Lplus);
        ub := Floor((Sqrt(sqr_bd) - c_mu_plus - c_gamma_plus)*c_Lplus);
        for k in [lb..ub] do
            x := (c_gamma_plus + c_mu_plus + k * c_Lplus^(-1)) * lambda_rat;
            assert (m - (x*Qrat,x)/2) ge 0;
            // if (m - (x*Qrat,x)/2) lt 0 then printf "skipping...\n"; continue; end if;
            vprintf ShimuraQuotients, 5: "\n\t\t\tmu_minus = %o, m - Q(x) = %o\n", gamma_minus + mu_minus, m - (x*ChangeRing(Q,Rationals()),x)/2;
            norm_mu_minus := ((gamma_minus + mu_minus)*Qrat, gamma_minus + mu_minus)/2;
            vprintf ShimuraQuotients, 5: "\t\t\tQ(mu_minus) = %o, Q(mu_minus) - m + Q(x) = %o\n", norm_mu_minus, norm_mu_minus - m + (x*ChangeRing(Q,Rationals()),x)/2;
            if (m - (x*Qrat,x)/2 eq 0) then // and (gamma_minus + mu_minus ne 0) then // This condition is for Chowla-Selberg constant
                if (gamma ne 0) then
                    Yang_tt := true;
                else
                    m0, m_cond := SquareFreeFactorization(Integers()!m);
                    fac := Factorization(m_cond);
                    for pe in fac do
                        p,e := Explode(pe);
                        log_coeffs -:= LogSum(Rationals()!2*e,p);
                    end for;
                end if;
            else
                a, p := kappaminus(gamma_minus + mu_minus, m - (x*Qrat,x)/2, Lminus, Q, d);
                log_coeffs +:= LogSum(Rationals()!a,p);
            end if;
        end for;
    end for;

    // trying to imitate Yang's code
    // !!! Don't know why this is working !!!
    if Yang_tt then
        d0 := FundamentalDiscriminant(d);
        f2 := d div d0;
        is_pp, p, e2 := IsPrimePower(f2);
        if is_pp then
            e := e2 div 2;
            log_coeffs +:= LogSum(-4*p^(1-e)/(p-KroneckerSymbol(d0,p)),p);
        end if;
    end if;

    return log_coeffs;
end intrinsic;

intrinsic SchoferFormula(fs::SeqEnum[RngSerLaurElt], d::RngIntElt, Q::AlgMatElt, lambda::ModTupRngElt, scale::FldRatElt) -> SeqEnum[LogSm]
{Assuming that fs are q-expansions of oo-weakly holomorphic modular forms at oo, 
 returns the log of the absolute value of Psi_F_f for every f in fs at the CM point with CM d.
 Here Q is the Gram matrix of the lattice L and lambda is a vecotr of norm -d.}
    ns := [-Valuation(f) : f in fs];
    n := Maximum(ns);
    log_coeffs := [LogSum() : f in fs];
    for m in [1..n] do
        if &and[Coefficient(f, -m) eq 0 : f in fs] then continue; end if;
        log_coeffs_m := Kappa0(m,d,Q,lambda);
        vprintf ShimuraQuotients, 5 : "\t\t";
        vprintf ShimuraQuotients, 4 : " is %o", log_coeffs_m;
        for i->f in fs do
            log_coeffs[i] +:= Coefficient(f,-m)*log_coeffs_m;
        end for;
    end for;

    // rescaling

    for i in [1..#fs] do
        log_coeffs[i] := scale * log_coeffs[i];
    end for;

    return log_coeffs;
end intrinsic;

function SchoferFormula0(fs_0, d, Q, lambda_v, scale, M, disc_grp, to_disc)

    log_coeffs := [LogSum() : f in fs_0];

    ns := [-Valuation(f) : f in fs_0];
    n := Maximum(ns);
    
    // computing norms of elements in the discriminant group
    mod_M_to_vecs := AssociativeArray([0..M-1]);
    for j in [0..M-1] do
        mod_M_to_vecs[j] := [];
    end for;
    for eta in disc_grp do
        v := ChangeRing(eta@@to_disc,Rationals());
        norm_v := (v*Q,v)/(2*M);
        if not IsIntegral(norm_v) then continue; end if;
        norm_mod_M := Integers()!norm_v mod M;
        Append(~mod_M_to_vecs[norm_mod_M], eta);
    end for;

    for mM in [1..n] do
        if &and[Coefficient(f, -mM) eq 0 : f in fs_0] then continue; end if;
        gammas:= [1/M*ChangeRing(gammaM@@to_disc, Rationals()) : gammaM in mod_M_to_vecs[mM mod M]];
        log_coeffs_m := &+([Kappa(gamma,mM/M,d,Q,lambda_v) : gamma in gammas] cat [LogSum()]);
        vprintf ShimuraQuotients, 4 : " is %o", log_coeffs_m;
        for i->f in fs_0 do
            log_coeffs[i] +:= Coefficient(f,-mM)*log_coeffs_m;
        end for;
    end for;

    // rescaling

    for i in [1..#fs_0] do
        log_coeffs[i] := scale * log_coeffs[i];
    end for;

    return log_coeffs;
end function;

intrinsic SchoferFormula(f::RngSerLaurElt, d::RngIntElt, Q::AlgMatElt, lambda::ModTupRngElt, scale::FldRatElt) -> LogSm
{Assuming that f is the q-expansions of a oo-weakly holomorphic modular form at oo, 
 returns the log of the absolute value of Psi_F_f at the CM point with CM d.
 Here Q is the Gram matrix of the lattice L and lambda is a vecotr of norm -d.}
    return SchoferFormula([f], d, Q, lambda, scale)[1];
end intrinsic;

intrinsic ScaleForSchofer(d::RngIntElt, D::RngIntElt, N::RngIntElt) -> FldRatElt
{Return the scaling factor in Schofer formula for CM(d) on X*(D,N).}
    D0 := D div 2^Valuation(D,2);
    M := 4*D0;
    
    OK := MaximalOrder(QuadraticField(d));
    is_sqr, cond := IsSquare(d div Discriminant(OK));
    assert is_sqr;
    // require cond eq 1 : "Not implemented for non-maximal orders!";
    O := sub<OK | cond>;
    n_d := NumberOfOptimalEmbeddings(O, D, N);
    require n_d gt 0 : "Curve does not have a CM point of discirminant d!";
    W_size := 2^#PrimeDivisors(D*N);

    // Not sure?? Think this what happens to the number of CM points on the full quotient
    /*
    sqfree, sq := SquarefreeFactorization(d);
    Ogg_condition := (cond eq 1) or (((sqfree mod 4 eq 1) and (cond eq 2)));
    if ((D*N) mod sqfree eq 0) and Ogg_condition then
        W_size div:= 2;
    end if;
    */
    // This follows from Ogg's description of the fixed points 
    // of Atkin-Lehner w_m
    Ogg_condition := ((d eq -4) and IsEven(D*N)) or
                     ((d mod 4 eq 0) and ((D*N mod (d div 4)) eq 0)) or
                     ((d mod 4 eq 1) and (D*N mod d eq 0));
    if Ogg_condition then
        W_size div:= 2;
    end if;

    scale := -n_d / (4*W_size);

    return scale;
end intrinsic;

intrinsic SchoferFormula(f::RngSerLaurElt, d::RngIntElt, D::RngIntElt, N::RngIntElt, Ldata::QuaternionLatticeData : Lambda := false) -> LogSm
{Assuming that f is the q-expansions of a oo-weakly holomorphic modular form at oo, 
 returns the log of the absolute value of Psi_F_f at the CM point with CM d.}
    // _,_,_,_,_,Q,O,basis_L := ShimuraCurveLattice(D,N);
    Q := Ldata`Q;
    O := Ldata`O;
    basis_L := Ldata`basis_L;

    scale := ScaleForSchofer(d,D,N);

     if Type(Lambda) eq BoolElt then 
        lambda := ElementOfNorm(Q,-d, O, basis_L);
    else
        lambda := Lambda;
    end if;

    return SchoferFormula(f, d, Q, lambda, scale);
end intrinsic;


// !! TODO - cache the Kappa0's or do it for a bunch of fs simultaneously
// We use a combination of the two versions of Schofer's formula from [GY] and [Err]
// We write sum log|psi|^2 = -|CM(d)|/4 * sum c_m kappa(-m)
// Note that in [GY] there is no square on the lhs, and 
// in [Err] there is no division by 4 on the rhs,
// but this seems to match with the examples in [Err] !?
intrinsic SchoferFormula(etas::SeqEnum[EtaQuot], d::RngIntElt, D::RngIntElt, N::RngIntElt, Ldata::QuaternionLatticeData : Lambda := false) -> SeqEnum[LogSm]
{Return the log of the absolute value of Psi_F_f for every f in fs at the CM point with CM d.}
    // _,_,disc_grp,to_disc,_, Q, O, basis_L := ShimuraCurveLattice(D,N);
    Q := Ldata`Q;
    O := Ldata`O;
    basis_L := Ldata`basis_L;
    disc_grp := Ldata`disc_grp;
    to_disc := Ldata`to_disc;

    scale := ScaleForSchofer(d,D,N);
    if Type(Lambda) eq BoolElt then 
        lambda := ElementOfNorm(Q, -d,  O, basis_L);
    else
        lambda := Lambda;
    end if;

    fs := [qExpansionAtoo(eta,1) : eta in etas];
    fs_0 := [qExpansionAt0(eta,1) : eta in etas];

    // Taking care of the principal part at infinity
    log_coeffs := SchoferFormula(fs, d, Q, lambda, scale);

    // Taking care of the principal part at zero
    M := IsOdd(D*N) select 4*D*N else 2*D*N;
    log_coeffs_0 := SchoferFormula0(fs_0, d, Q, lambda, scale, M, disc_grp, to_disc);

    // summing up
    for i->s in log_coeffs do
        log_coeffs[i] +:= log_coeffs_0[i];
    end for;

    return log_coeffs;
end intrinsic;

intrinsic SchoferFormula(eta::EtaQuot, d::RngIntElt, D::RngIntElt, N::RngIntElt, Ldata::QuaternionLatticeData : Lambda := false) -> LogSm
{Return the log of the absolute value of Psi_F_f at the CM point with CM d.}
    return SchoferFormula([eta], d, D, N, Ldata : Lambda := Lambda)[1];
end intrinsic;

intrinsic AbsoluteValuesAtRationalCMPoint(fs::SeqEnum[EtaQuot], d::RngIntElt, Xstar::ShimuraQuot, Ldata::QuaternionLatticeData : Lambda := false) -> SeqEnum[LogSm]
{Returns the absolute value of f for every f in fs at the rational CM point with CM d.}
    vals := [LogSum() : f in fs];
    for i->f in fs do
        div_f := DivisorOfBorcherdsForm(f, Xstar);
        in_support := exists(pt){pt : pt in div_f | pt[1] eq d};
        if in_support then
            if pt[2] lt 0 then vals[i] := LogSum(Infinity()); end if;
            if pt[2] gt 0 then vals[i] := LogSum(0); end if;
        end if;
    end for;
    rest_idxs := [i : i in [1..#fs] | vals[i] eq LogSum()];
    if IsEmpty(rest_idxs) then return vals; end if;
    rest_fs := [fs[i] : i in rest_idxs];
    log_coeffs := SchoferFormula(rest_fs, d, Xstar`D, Xstar`N, Ldata : Lambda := Lambda);
    for i->log_coeff in log_coeffs do
        vals[rest_idxs[i]] := log_coeff;
    end for;
    return vals;
end intrinsic;

function coeffs_to_divisor_matrix(min_m, D, N, num_coeffs : Zero := false, const_coeff := true)
    disc_ms := AssociativeArray();
    scale := Zero select N^2 else 4;
    for d in [1..-min_m] do
        discs := [scale*d div r2: r2 in Divisors(scale*d) | IsSquare(r2)];
        discs := [disc : disc in discs | disc mod 4 in [0,3]];
        for disc in discs do
            S := QuadraticOrder(BinaryQuadraticForms(-disc));
            n_d := NumberOfOptimalEmbeddings(S,D,N);
            if (n_d ne 0) then
                if not IsDefined(disc_ms, disc) then disc_ms[disc] := []; end if;
                Append(~disc_ms[disc], d);
            end if;
        end for;
    end for;
    relevant_ds := Sort([d : d in Keys(disc_ms)]);
    mat := ZeroMatrix(Integers(), num_coeffs, #relevant_ds + (const_coeff select 1 else 0)); // also the constant coefficient
    for j->d in relevant_ds do
        for m in disc_ms[d] do
            // collecting conributions for all m
            mat[1 - min_m - m,j] := 1;
        end for;
    end for;

    // consant coefficient
    if const_coeff then mat[Nrows(mat),Ncols(mat)] := 1; end if;

    return mat, relevant_ds;
end function;

function basis_of_weakly_holomorphic_forms(pole_order, fs_E, n0, n, t : Zero := false)
   
    k := -Valuation(Zero select qExpansionAt0(t,1 : Admissible := false) else qExpansionAtoo(t,1));
   
    r := (pole_order - n0) div k;
    s := pole_order - r*k;

    assert n + 2 gt n0 + k; // making sure we have enough forms to complete to a basis

    assert n + 1 - n0 lt #fs_E; // Making sure the value of n makes sense

    basis_n0 := fs_E[n+2-n0..#fs_E]; // basis for M_{n0-1}^!
    init_basis := fs_E[n+2-n0-k..n+1-n0]; // completing to a basis for M_{n_0+k-1}^!
   
    full_basis := [t^r*f : f in init_basis[n0+k-s..#init_basis]];
    full_basis cat:= &cat[[t^(r-1-j)*f : f in init_basis] : j in [0..r-1]];
    full_basis cat:= basis_n0;
    
    minval := pole_order;
   
    if Zero then
        qexps := [qExpansionAt0(eta, 1) : eta in full_basis];
    else
        qexps := [qExpansionAtoo(eta, 1) : eta in full_basis];
    end if;
    Rq<q> := Universe(qexps);
    R := BaseRing(Rq);
    assert minval eq -Minimum([Valuation(f) : f in qexps]);
   
    coeffs := Matrix(R, [AbsEltseq(q^minval*f : FixedLength) : f in qexps]);
    
    ech_basis, T := EchelonForm(coeffs);
    ech_etas := [&+[T[i][j]*full_basis[j] : j in [1..Ncols(T)] | T[i][j] ne 0] : i in [1..Nrows(T)]];
   
    return ech_basis, ech_etas, T;
end function;

function sum_divisors(div1, div2)
    // sum the divisors
    sum_divs := div1;
    for pt in div2 do
        if exists(j){j : j->x in sum_divs | x[1] eq pt[1]} then
            sum_divs[j][2] +:= pt[2];
        else
            Append(~sum_divs, pt);
        end if;
    end for;
    nonzero_idxs := [j : j->pt in sum_divs | pt[2] ne 0];
    sum_divs := [sum_divs[j] : j in nonzero_idxs];
    return sum_divs;
end function;

intrinsic BorcherdsForms(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot] : Prec := 100, Exclude := {}) -> Assoc
{Returns weakly holomorphic modular forms with divisors that are the ramification divisors of each of the double covers in curves,
along with two different hauptmoduls.}
    rams := RamficationPointsOfCovers(Xstar, curves);
    D0,M,g := get_D0_M_g(Xstar`D,Xstar`N);
    
    // !!! For D = 35, this takes about 13 minutes on lava....
    E, n, n0, t, eta_quotients := WeaklyHolomorphicBasis(Xstar`D, Xstar`N : Prec := Prec);
    k := -Valuation(qExpansionAtoo(t,1));
   
    if IsOdd(Xstar`D*Xstar`N) then
        E0, nE0, _, eta_quotients_oo, eta_quotients_0 := WeaklyHolomorphicBasis(Xstar`D, Xstar`N : Prec := Prec, Zero, n0 := n0);
    end if;
    // we do this twice -- we should remember this
    pts := RationalCMPoints(Xstar); // pts <-> infty, 0, rational
    pts := [p : p in pts | p[1] notin Exclude and GCD(p[1], Xstar`N) eq 1];

    max_pole_order_oo := 0;
    ech_basis_all_oo :=  MatrixAlgebra(Rationals(),0)!0; // zero matrix
    ech_etas_all_oo := [];
    T_all_oo := MatrixAlgebra(Rationals(),0)!0; // zero matrix

    max_pole_order_0 := 0;
    ech_basis_all_0 :=  MatrixAlgebra(Rationals(),0)!0; // zero matrix
    ech_etas_all_0 := [];
    T_all_0 := MatrixAlgebra(Rationals(),0)!0; // zero matrix

    all_ms := [];
    m_idx := 1;
    all_ms := &cat[[(d[1] mod 4 eq 0) select d[1] div 4 else d[1] : d in pts] : pt in pts];
    all_ms := Reverse(Sort([m : m in Set(all_ms)]));

    found_all := false;
    
    while (not found_all) do
        if IsOdd(Xstar`D*Xstar`N) then
            vprintf ShimuraQuotients, 2 : "\n\tAttempting to find Borcherds forms with m = %o...", all_ms[m_idx];
        end if;
        for infty in pts do
            vprintf ShimuraQuotients, 2 : "\n\tTrying infinity = %o...", infty;
            non_infty := [pt : pt in pts | pt ne infty];
            for other_pts in CartesianPower(non_infty,2) do
                if other_pts[1] eq other_pts[2] then continue; end if;
                vprintf ShimuraQuotients, 3 : "\n\t\tTrying other points = %o...", other_pts;
                rams[-1] := [other_pts[1]];
                rams[-2] := [other_pts[2]];
                
                etas := AssociativeArray();
                
                found_all := true;
                for i in Keys(rams) do
                    ram := rams[i];
                    // adding the part at infinity
                    if exists(j){j : j->pt in ram | pt[1] eq infty[1]} then
                        assert ram[j] eq infty;
                        Remove(~ram, j);
                    end if;
                    deg := &+[pt[3] : pt in ram];
                    div_coeffs := [1 : pt in ram] cat [-deg]; // divisor coefficients
                    Append(~ram, infty);

                    vprintf ShimuraQuotients, 4 : "\n\t\t\tWorking on ramification divisor %o...", [<pt[1], div_coeffs[j]> : j->pt in ram];

                    ms := [(d[1] mod 4 eq 0) select d[1] div 4 else d[1] : d in ram];
                    min_m := Minimum(ms);
                    min_m := Minimum(min_m, -(n0 + k - 1));
                    
                    if (max_pole_order_oo lt -min_m) then
                        max_pole_order_oo := -min_m;
                        vprintf ShimuraQuotients, 5 : "\n\t\t\t\tComputing basis of {oo}-weakly holomorphic forms with pole order %o...", -min_m;
                        ech_basis_all_oo, ech_etas_all_oo, T_all_oo := basis_of_weakly_holomorphic_forms(-min_m, eta_quotients, n0+1, n, t);
                        vprintf ShimuraQuotients, 5 : "Done!";
                    end if;
                    
                    first_idx := min_m+max_pole_order_oo+1;
                    ech_basis := SubmatrixRange(ech_basis_all_oo, first_idx, first_idx, Nrows(ech_basis_all_oo), Ncols(ech_basis_all_oo));
                    ech_etas := ech_etas_all_oo[first_idx..#ech_etas_all_oo];
                    assert SubmatrixRange(T_all_oo, first_idx, 1, Nrows(T_all_oo), first_idx-1) eq 0;
                    T := SubmatrixRange(T_all_oo, first_idx, first_idx, Nrows(T_all_oo), Ncols(T_all_oo));

                    if IsOdd(Xstar`D*Xstar`N) then
                        assert m_idx le #all_ms;
                        m_choice := all_ms[m_idx];
                        vprintf ShimuraQuotients, 5 : "\n\t\t\t\tWorking on m = %o for q-expansion at 0", m_choice;
                        pole_order := -D0*m_choice;
                    
                        if (max_pole_order_0 lt pole_order) then
                            max_pole_order_0 := pole_order;
                            t0 := SAction(t : Admissible := false);
                            vprintf ShimuraQuotients, 5 : "\n\t\t\t\tComputing basis of {0,oo}-weakly holomorphic forms with pole orders (%o, %o)...", pole_order/(4*D0), nE0;
                            ech_basis_all_0, ech_etas_all_0, T_all_0 := basis_of_weakly_holomorphic_forms(pole_order, eta_quotients_oo, 1, nE0, t0 : Zero);
                            vprintf ShimuraQuotients, 5 : "Done!";
                        end if;

                        first_idx := -pole_order+max_pole_order_0+1;
                        ech_basis_0 := SubmatrixRange(ech_basis_all_0, first_idx, first_idx, Nrows(ech_basis_all_0), Ncols(ech_basis_all_0));
                        ech_etas_0 := ech_etas_all_0[first_idx..#ech_etas_all_0];
                        assert SubmatrixRange(T_all_0, first_idx, 1, Nrows(T_all_0), first_idx-1) eq 0;
                        T0 := SubmatrixRange(T_all_0, first_idx, first_idx, Nrows(T_all_0), Ncols(T_all_0));

                        vprintf ShimuraQuotients, 5 : "\n\t\t\t\tBuilding q-expansions at oo...";
                        ech_fs_oo := [qExpansionAtoo(eta,1) : eta in ech_etas_0];
                        vprintf ShimuraQuotients, 5 : "Done!";

                        Rq<q> := Universe(ech_fs_oo);
                        R := BaseRing(Rq);
                
                        ech_basis_oo := Matrix(R, [AbsEltseq(q^n0*f : FixedLength) : f in ech_fs_oo]);
                        
                        non_div_idxs := [i : i in [1..Ncols(ech_basis_0)] | (i-1-pole_order) mod D0 ne 0];
                        div_idxs := [i : i in [1..Ncols(ech_basis_0)] | (i-1-pole_order) mod D0 eq 0];
                        // good_forms_0 := BasisMatrix(Kernel(Submatrix(ech_basis_0, [1..Nrows(ech_basis_0)], non_div_idxs)));
                        T := BasisMatrix(Kernel(Submatrix(ech_basis_0, [1..Nrows(ech_basis_0)], non_div_idxs)));
                        good_forms_0 := T*ech_basis_0;
                        assert Submatrix(good_forms_0,[1..Nrows(good_forms_0)], non_div_idxs) eq 0;
                        // T := Solution(ech_basis_0, good_forms_0);
                        // assert T*ech_basis_0 eq good_forms_0;
                        good_forms_oo := ChangeRing(T,Rationals())*ech_basis_oo;
                        // Passingt o q-expansions with q^(1/4) instead of q^(1/4D0)
                        good_forms_0 := Submatrix(good_forms_0,[1..Nrows(good_forms_0)], div_idxs);
                        // This was now verified to give the q-expansion of h in [GY] Example 31, p. 20 
                        mat_0, relevant_ds_0 := coeffs_to_divisor_matrix(m_choice, Xstar`D, Xstar`N, Ncols(good_forms_0) : Zero, const_coeff := false);
                        mat_oo, relevant_ds_oo := coeffs_to_divisor_matrix(-n0, Xstar`D, Xstar`N, Ncols(good_forms_oo) : const_coeff := false);
                        coeffs_0 := good_forms_0*ChangeRing(mat_0, Rationals());
                        coeffs_oo := good_forms_oo*ChangeRing(mat_oo,Rationals());

                        ech_etas_0 := [&+[T[i][j]*ech_etas_0[j] : j in [1..#ech_etas_0]] : i in [1..Nrows(T)]];

                        // collecting contributions from 0 and oo
                        relevant_ds_0_oo := Sort([x : x in Set(relevant_ds_0) join Set(relevant_ds_oo)]);

                        ds_0_to_ds := ZeroMatrix(Integers(), #relevant_ds_0, #relevant_ds_0_oo);
                        for i->d in relevant_ds_0 do
                            ds_0_to_ds[i, Index(relevant_ds_0_oo, d)] := 1;
                        end for;
                        
                        ds_oo_to_ds := ZeroMatrix(Rationals(), #relevant_ds_oo, #relevant_ds_0_oo);
                        for i->d in relevant_ds_oo do
                            ds_oo_to_ds[i, Index(relevant_ds_0_oo, d)] := 1;
                        end for;
                        
                        mat_0_oo := coeffs_0*ChangeRing(ds_0_to_ds,Rationals()) + coeffs_oo*ds_oo_to_ds;
                    end if;

                    mat, relevant_ds := coeffs_to_divisor_matrix(min_m, Xstar`D, Xstar`N, Ncols(ech_basis));
                    coeffs_trunc := ech_basis * ChangeRing(mat, BaseRing(ech_basis));

                    if IsOdd(Xstar`D*Xstar`N) then
                        ds_0_oo_to_ds := ZeroMatrix(Rationals(), #relevant_ds_0_oo, #relevant_ds + 1);
                        for i->d in relevant_ds_0_oo do
                            ds_0_oo_to_ds[i, Index(relevant_ds, d)] := 1;
                        end for;
                        coeffs_0_oo := mat_0_oo*ds_0_oo_to_ds;
                        coeffs_trunc := VerticalJoin(ChangeRing(coeffs_trunc,Rationals()), coeffs_0_oo);
                    end if;

                    V := RSpace(BaseRing(coeffs_trunc), Ncols(mat));
                    target_v := &+[div_coeffs[j]*pt[2]*V.(Index(relevant_ds,-pt[1])) : j->pt in ram];
                    
                    found_v := target_v in Image(coeffs_trunc);
                   
                    if not found_v then found_all := false; break; end if;
                    sol := Solution(coeffs_trunc, target_v);
                    
                    etas[i] := &+[sol[i]*ech_etas[i] : i in [1..#ech_etas]];
                    if IsOdd(Xstar`D*Xstar`N) then
                        etas[i] +:= &+[sol[#ech_etas + i]*ech_etas_0[i] : i in [1..#ech_etas_0]];
                    end if;
                    // check divisor
                    div_f := DivisorOfBorcherdsForm(etas[i], Xstar);
                    
                    assert Set(div_f) eq {<pt[1], div_coeffs[j]> : j->pt in ram};
                end for;
                if found_all then break; end if;
            end for;
            if found_all then break; end if;
        end for;
        m_idx +:= 1;
        if (m_idx gt #all_ms) then break; end if;
    end while;  
    vprintf ShimuraQuotients, 2 : "\n";
    if not found_all then
        error "Failed to find all Borcherds forms";
    end if;
    return etas;
end intrinsic;

intrinsic DivisorOfBorcherdsForm(f::RngSerLaurElt, Xstar::ShimuraQuot : Zero := false) -> SeqEnum
{Return the divisor of the Borcherds form associated to the (oo)-Weakly holomorphic modular form f.
If Zero is set to true, returns the divisor associated to the (0)-Weakly holomorphic modular form with q-expansion f(q^(1/4))}
    if IsWeaklyZero(f) then return []; end if;
    N := Valuation(f);
    ells := NumberOfEllipticPointsByCMOrder(Xstar);
    ell_lut := AssociativeArray();
    for q in Keys(ells) do
        for d in Keys(ells[q]) do
            ell_lut[d] := q;
        end for;
    end for;
    divisor := [];
    scale := Zero select (Xstar`N)^2 else 4;
    for m in [N..-1] do
        c_m := Coefficient(f,m);
        if c_m eq 0 then continue; end if;
        sq_divs := [r2 : r2 in Divisors(-scale*m) | IsSquare(r2)];
        for r2 in sq_divs do
            d := (scale*m) div r2;
            if (d mod 4 notin [0,1]) then continue; end if;
            S := QuadraticOrder(BinaryQuadraticForms(d));
            n_d := NumberOfOptimalEmbeddings(S,Xstar`D,Xstar`N);
            if n_d eq 0 then continue; end if;
            e := 1;
            if IsDefined(ell_lut, d) then
                e := ell_lut[d];
            end if;
            Append(~divisor, <d, c_m/e>);
        end for;
    end for;

    simple_divisor := [];
    discs := {pair[1] : pair in divisor};
    for d in discs do
        mult_d := &+[pair[2] : pair in divisor | pair[1] eq d];
        if mult_d ne 0 then
            Append(~simple_divisor, <d, mult_d>);
        end if;
    end for;
    return simple_divisor;
end intrinsic;

// This is [GY, Lemma 25]
intrinsic DivisorOfBorcherdsForm(foo::RngSerLaurElt, f0::RngSerLaurElt, Xstar::ShimuraQuot) -> SeqEnum
{Return the divisor of the Borcherds form associated to the (0,oo)-Weakly holomorphic modular form f with q-expansion foo(q) at oo and f0(q^(1/4)) at 0.}
    div_oo := DivisorOfBorcherdsForm(foo, Xstar);
    div_0 := DivisorOfBorcherdsForm(f0, Xstar : Zero);
   
    return sum_divisors(div_oo, div_0);
end intrinsic;

intrinsic DivisorOfBorcherdsForm(f::EtaQuot, Xstar::ShimuraQuot) -> SeqEnum
{Return the divisor of the Borcherds form associated to the (0,oo)-Weakly holomorphic modular form f.}
    foo := qExpansionAtoo(f,0);
    f0_D0 := qExpansionAt0(f,0);
    D0 := Xstar`D div 2^Valuation(Xstar`D,2);
    v := Valuation(f0_D0) div D0;
    coeffs_f0_D0 := [Coefficient(f0_D0,m*D0) : m in [v..-1]];
    R<q> := Parent(f0_D0);
    f0 := q^v*&+[R | coeffs_f0_D0[i]*q^(i-1) : i in [1..#coeffs_f0_D0]];

    vprintf ShimuraQuotients, 4 : "\n\t\t\tComputing divisor of %o,", f;
    vprintf ShimuraQuotients, 4 : " qExpansion at 0 is %o...", f0;

    ret := DivisorOfBorcherdsForm(foo, f0, Xstar);
    vprintf ShimuraQuotients, 4 : "Done!";
    return ret;
end intrinsic;


intrinsic CandidateDiscriminants(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot] : Exclude := {}, bd := 4) -> SeqEnum
{Returns list of candidate discriminats for Schofer's formula} //'
    cm_pts := RationalCMPoints(Xstar : Exclude := Exclude);
    cm_pts := Reverse(Sort(cm_pts));

    quad_cm := QuadraticCMPoints(Xstar : Exclude := Exclude, bd := bd);
    quad_cm := Reverse(Sort(quad_cm));

    return [cm_pts, quad_cm];

end intrinsic;

intrinsic AbsoluteValuesAtCMPoints(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot], all_cm_pts ::SeqEnum, fs::Assoc : MaxNum := 7, Prec := 100,  Exclude := {}, Include := {}) -> SchoferTable, SeqEnum
{Returns the absolute values of y^2 for all degree 2 covers and two hauptmodules at CM points.}
    
    cm_pts := [a : a in all_cm_pts[1]| a[1] notin Exclude];
    quad_cm := [a : a in  all_cm_pts[2] | a[1] notin Exclude];
    all_cm_pts := [cm_pts, quad_cm];

    keys_fs := [k : k in Keys(fs)];
    all_fs := [fs[k] : k in keys_fs];
   
    if #Include gt 0 then
        bd := Maximum([ ClassNumber(d) : d in Include]);
    else
        bd:= 4;
    end if;

    cm_pts_must := [p : p in cm_pts | p[1] in Include];
    other_cm := [p : p in cm_pts | p[1] notin Include];
    need := MaxNum - #Include;
    // if #cm_pts gt need then  
    if #other_cm ge need then
    //need to make space for include points, but otherwise fill up with rational points as much as possible
        // pt_list_rat := cm_pts_must cat other_cm[1..(need - #cm_pts_must)];
        // The above does not make sense - we need to complete to MaxNum points
        pt_list_rat := cm_pts_must cat other_cm[1..need];
    else
        pt_list_rat := cm_pts_must cat other_cm;
    end if;

    //remove points that we already included
    Include := SetToSequence(Include);
    for p in cm_pts_must do
        pidx := Index(Include, p[1]);
        Include := Remove(Include,Index(Include, p[1]));
    end for;

    if #Include gt 0 then
        bd := Maximum([ ClassNumber(d) : d in Include]);
    else
        bd:= 2;
    end if;
    need := MaxNum - #pt_list_rat;
    //Whatever is left, we need

    cm_pts_must := [p : p in quad_cm | p[1] in Include];
    other_cm := [p : p in quad_cm | p[1] notin Include];

    pt_list_quad := cm_pts_must cat other_cm[1..(need -#Include)];

    for p in cm_pts_must do
        pidx := Index(Include, p[1]);
        Include := Remove(Include,Index(Include, p[1]));
    end for;

    assert #Include eq 0;

    table := [[] : f in all_fs];
    // _,_,_,_,_,Q,O,basis_L := ShimuraCurveLattice(Xstar`D,Xstar`N);
    Ldata := ShimuraCurveLattice(Xstar`D,Xstar`N);
    Q := Ldata`Q;
    O := Ldata`O;
    basis_L := Ldata`basis_L;

    lambdas := ElementsOfNorm(Q, [-pt[1] : pt in pt_list_rat cat pt_list_quad], O, basis_L);
    for pt in pt_list_rat do
        d := pt[1];
        vals := AbsoluteValuesAtRationalCMPoint(all_fs, d, Xstar, Ldata : Lambda := lambdas[-d]);
        for i->v in vals do
            Append(~table[i], vals[i]);
        end for;
    end for;

    for pt in pt_list_quad do
        d := pt[1];
        norm_val := AbsoluteValuesAtRationalCMPoint(all_fs, d, Xstar, Ldata : Lambda := lambdas[-d]);
        for i->v in norm_val do
            Append(~table[i], norm_val[i]);
        end for;
    end for;

    ds := [[pt[1] : pt in pt_list_rat], [ pt[1] : pt in pt_list_quad ]];

    schofer_tab := CreateSchoferTable(table, keys_fs, ds, curves,Xstar);
    schofer_tab`BorcherdsForms := fs;

    return schofer_tab, all_cm_pts;
end intrinsic;

function find_signs(s, stilde, ds)
    ratds := ds[1];
    quadds := ds[2];
    inf_zero_indices := [Index(s,0), Index(stilde,0), Index(s,Infinity())];
    assert stilde[inf_zero_indices[3]] eq Infinity();
    scale_tilde := stilde[Index(s,0)];
    scale := s[Index(stilde,0)];
    idxs := [i : i in [1..#s] | i notin inf_zero_indices and i le #ratds];
    signs := &cat[[[eps1, eps2] : eps1,eps2 in [-1,1] | eps1*s[i]/scale + eps2*stilde[i]/scale_tilde eq 1] : i in idxs];
    degs := [1 : i in ds[1] ] cat [ 2 : i in ds[2]];
    s_new := [* ss/scale^degs[i] : i->ss in s *];
    stilde_new := [* sstilde/scale_tilde^degs[i] : i->sstilde in stilde *];
    for j->idx in idxs do
        s_new[idx] := signs[j][1]*s_new[idx];
        stilde_new[idx] := signs[j][2]*stilde_new[idx];
    end for;
    return s_new, stilde_new, scale, scale_tilde;
end function;


intrinsic RationalConstraintsOnEquations(schofer_table::SchoferTable, curves::SeqEnum[ShimuraQuot]) -> SeqEnum
{Impose linear constraints on equations given by rational CM points}

    keys_fs := schofer_table`Keys_fs;
    table := schofer_table`Values;
    sep_ds := schofer_table`Discs;
    ds := sep_ds[1] cat sep_ds[2];
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


procedure add_new_column(schofer_tab, dnew, deg)
    //add column associated to dnew of degree deg
    table := schofer_tab`Values;
    ds := schofer_tab`Discs;
    Xstar := schofer_tab`Xstar;
    fs := schofer_tab`BorcherdsForms;
    keys_fs := schofer_tab`Keys_fs;
    row_scales := schofer_tab`RowScales;
    curves := schofer_tab`Curves;
    all_fs := [fs[k] : k in keys_fs];
    Ldata := ShimuraCurveLattice(Xstar`D,Xstar`N);
    norm_val := AbsoluteValuesAtRationalCMPoint(all_fs, dnew, Xstar, Ldata);
    for i->v in norm_val do
        Append(~table[i], norm_val[i]/row_scales[i]^deg);
    end for;
    schofer_tab`Values := table;
    if deg eq 1 then
        Append(~ds[1],dnew);
    else
        Append(~ds[2],dnew);
    end if;
    schofer_tab`Discs := [ds[1], ds[2]];
    return;
end procedure;

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
    ds := ds[1] cat ds[2];
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

intrinsic EquationsAbovePointlessConics(all_eqns::Assoc, all_ws::Assoc, curves::SeqEnum) -> Assoc, Assoc
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
            if exists(base){b : b in Keys(all_eqns[other_curve])} then // all eqns for all bases have the same degree
                if (Degree(HyperellipticPolynomials(all_eqns[other_curve][base])) eq g+1) then
                    gplus1key := other_curve; //found the gplus1
                    found_gplus1 := true;
                    break;
                end if;
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
            im_y := Evaluate(alg_map_gplus1(y_var), [s,y,z]);
            im_s1 := Evaluate(alg_map_gplus1(s_var), [s,y,z]);
            w_mapconic := all_ws[conic_key][base][w];
            alg_map_conic := AlgebraMap(w_mapconic);
            s_var := Domain(alg_map_conic).1;
            x_var := Domain(alg_map_conic).2;
            im_x :=  Evaluate(alg_map_conic(x_var), [s,x,z]);
            im_s2 := Evaluate(alg_map_conic(s_var), [s,y,z]);
            assert im_s1 eq im_s2; //acts the same way
            wmap := map<C->C | [im_x, im_y, im_s1, z]>;
            b, inv := IsIsomorphism(wmap);
            assert b;
            all_ws[k][conic_key][w] := inv^(-1);
        end for;
    end for;
    return all_eqns, all_ws;

end intrinsic;

intrinsic AllEquationsAboveCovers(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot] : Prec := 100)-> Assoc, Assoc
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
    vprintf ShimuraQuotients, 1 : "Done!\n";
    genus_list := [curves[i]`g: i in Xstar`CoveredBy];
    num_vals := Maximum([2*g+5 : g in genus_list]);
    vprintf ShimuraQuotients, 1 : "Computing absolute values at CM points...";
    Exclude := {pt[1] : pt in all_cm_pts[1] cat all_cm_pts[2] | GCD(pt[1], Xstar`N) ne 1};
    abs_schofer_tab, all_cm_pts:= AbsoluteValuesAtCMPoints(Xstar, curves, all_cm_pts, fs : 
                                                           MaxNum := num_vals, Prec := Prec, 
                                                           Exclude := Exclude, Include := Set(d_divs));
    vprintf ShimuraQuotients, 2 : "\n";
    vprintf ShimuraQuotients, 1 : "Done!\n";
    ReduceTable(abs_schofer_tab);
    vprintf ShimuraQuotients, 1 : "Computing actual values at CM points...";
    schofer_tab := ValuesAtCMPoints(abs_schofer_tab, all_cm_pts);
    vprintf ShimuraQuotients, 1 : "Done!\n";  
    vprintf ShimuraQuotients, 1 : "Computing equations of covers...";
    crv_list, ws, new_keys := EquationsOfCovers(schofer_tab, all_cm_pts);
    vprintf ShimuraQuotients, 1 : "Done!\n";
    vprintf ShimuraQuotients, 1 : "Computing equations above P1s and conics...";
    all_eqns, all_ws := EquationsAboveP1s(crv_list, ws, new_keys, curves); //still adding ws here in the conic case
    vprintf ShimuraQuotients, 1 : "Done\n";
    vprintf ShimuraQuotients, 1 :"Computing equations above pointless conics...";
    all_eqns, all_ws := EquationsAbovePointlessConics(all_eqns, all_ws, curves);
    vprintf ShimuraQuotients, 1 : "Done\n";
    return all_eqns, all_ws;
end intrinsic;

intrinsic AllEquationsAboveCovers(D::RngIntElt, N::RngIntElt, curves::SeqEnum[ShimuraQuot] : Prec := 100)-> Assoc, Assoc
{Get equations of all covers (not just immediate covers)}
    _ := exists(Xstar){X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)};
    return AllEquationsAboveCovers(Xstar, curves : Prec := Prec);
end intrinsic;

// This is following [GR, Section 5]
intrinsic FieldsOfDefinitionOfCMPoint(X::ShimuraQuot, d::RngIntElt) -> List
{Return possible fields of definition for CM point with CM by d on X.}
    R := QuadraticOrder(BinaryQuadraticForms(d));
    K := NumberField(R);
    f := Conductor(R);
    H_R := RingClassField(R); // maybe want NumberField(H_R)
    D := X`D;
    N := X`N;
    D_R := &*[Integers()| p : p in PrimeDivisors(D) | KroneckerCharacter(d)(p) eq -1];
    N_R := &*[Integers()| p : p in PrimeDivisors(N) | KroneckerCharacter(d)(p) eq 1];   
    N_star_R := &*[Integers()| p : p in PrimeDivisors(N) | (KroneckerCharacter(d)(p) eq 1) and (f mod p ne 0)];
    W_R := [m : m in X`W | D_R*N_R mod m eq 0];
    assert GCD(D_R * N_star_R, Discriminant(R)) eq 1;
    assert GCD(D_R*N_R, Discriminant(R)) eq GCD(N,f);

    // Proposition 5.6
    if (Discriminant(R) mod ((D*N) div (D_R*N_star_R))) ne 0 then
        return [* *];
    end if;

    rec := ArtinMap(H_R);
    // rec_abs := Components(rec)[1];
    // gal_to_aut := &*Components(rec)[2..#Components(rec)];
    // assert rec_abs*gal_to_aut eq rec;

    // also number of points is 2^PrimeDivisors(D_R*N_R) * ClassNumber(R)

    // Theorem 5.8 - Shimura reciprocity
    // fields := [* F[1] : F in Subfields(AbsoluteField(NumberField(H_R))) *];
    // Q_P_ext := H_R;

    // setting up number fields
    H_R_NF := NumberField(H_R);
    abs_H_R := AbsoluteField(H_R_NF);
    _, H_R_to_abs := IsIsomorphic(H_R_NF, abs_H_R);

    // setting up Picard groups 
    PicR, mPicR := PicardGroup(R);
    A, PicR_to_A := PicR / (2*PicR);
    B := QuaternionAlgebra(D);

    al_action := AssociativeArray();

    // Theorem 5.12 (1) and Lemma 5.10 for complex conjugation
    m := D_R*N_star_R;
    //if m in X`W then
    for a in A do
        alift := a@@PicR_to_A;
        // circumventing a bug in Magma in mPicR
        if (alift eq PicR!0) then
            fraka := 1*R;
        else
            fraka := mPicR(alift);
        end if;
        B_fraka := QuaternionAlgebra(Rationals(), d, m*Norm(fraka));
        if IsIsomorphic(B_fraka, B) then
            break;
        end if;
    end for;
    assert assigned fraka;
    sigma_a := rec(fraka);
    abs_sig_a := H_R_to_abs^(-1)*sigma_a*H_R_to_abs;
    // _, K_to_abs := IsSubfield(K, abs_H_R);
    has_cc, cc := HasComplexConjugate(abs_H_R);
    if not has_cc then
        gal, auts, gal_to_auts := AutomorphismGroup(abs_H_R);
        // elements that restrict to the complex conjugation on K
        cc_candidates := [g : g in gal | Order(g) eq 2 and gal_to_auts(g)(K.1) eq ComplexConjugate(K.1)];  
        cc := Representative(cc_candidates);
        cc := gal_to_auts(cc);
    end if;
    sigma := hom<abs_H_R -> abs_H_R | cc(abs_sig_a(abs_H_R.1))>;
    if (m ne 1) then
        al_action[m] := sigma;
    else
        sigma_for_later := sigma;
    end if;

    // Lemma 5.9
    fixed_sub_gens := [];
    unknown_quotients := 0;
    // for m in ALsToGens(X`W, D*N) do
    als_DN := [Q : Q in Divisors(D*N) | GCD(Q, (D*N) div Q) eq 1];
    for m in als_DN do
        al_is_gal := ((D*N) div (D_R*N_R)) mod m eq 0;
        if al_is_gal then
            frakb := &*[Parent(1*Integers(K)) | pa[1]^(pa[2] div 2) : pa in Factorization(m*Integers(K))];
            assert Norm(frakb) eq m;
            // sigma := rec_abs(frakb);
            al_action[m] := H_R_to_abs^(-1)*rec(frakb)*H_R_to_abs;
            // Append(~fixed_sub_gens, sigma);
        // else
            // unknown_quotients +:= 1;
        end if;
        // !! TODO : figure out what to do if it is not Galois
    end for;

    known_al := Keys(al_action);
    // allws := {Integers()|};
    S := Subsets(known_al);
    for s in S do
        if #s eq 0 then
            al_action[1] := hom<abs_H_R->abs_H_R | abs_H_R.1>;
            // Include(~allws, 1);
        else
            prod := 1;
            for w in s do
                prev_prod := prod;
                prod := AtkinLehnerMul(w,prod,D*N);
                if not IsDefined(al_action, prod) then
                    al_action[prod] := al_action[prev_prod]*al_action[w];
                end if;
            end for;
            // Include(~allws, prod);
        end if;
    end for;

    if (m eq 1) then
        al_action[1] := sigma_for_later;
    end if;

    fixed_by := [al_action[m] : m in X`W meet Keys(al_action)];

    sanity_check, n_fixed := IsPowerOf(#fixed_by, 2);
    sanity_check_trivial, n_W := IsPowerOf(#X`W, 2);

    assert sanity_check;
    assert sanity_check_trivial;

    unknown_quotients := n_W - n_fixed;

    // fixed_sub := sub<Universe(fixed_sub_gens) | fixed_sub_gens>;
    // Q_P_ext := FixedField(H_R, fixed_sub);

    // fixed_by cat:= [H_R_to_abs^(-1)*gal_to_aut(s)*H_R_to_abs : s in fixed_sub_gens];
    Q_P_ext := FixedField(abs_H_R, fixed_by);

    // if Type(Q_P_ext) eq FldRat then return [* Q_P_ext *]; end if;

    // Q_Ps := [* F[1] : F in Subfields(Q_P_ext) | Degree(Q_P_ext) le 2^unknown_quotients * Degree(F[1]) *];

    return [* Q_P_ext *];

    // if (Degree(Q_P_ext) le 2^unknown_quotients) then
    //     Append(~Q_Ps, Rationals());
    // end if;

    /*
    if top_is_H then 
        Q_Ps := [* Q_P_ext *];
    else
        Q_Ps := [* Rationals() *] cat [* F[1] : F in Subfields(AbsoluteField(NumberField(Q_P_ext))) *];
    end if;
    */

    // return Q_Ps;
end intrinsic;

procedure replace_column(schofer_tab, d, dnew, is_log)
    //add column associated to dnew remove column associated to d
    table := schofer_tab`Values;
    ds := schofer_tab`Discs;
    Xstar := schofer_tab`Xstar;
    curveid := Xstar`CurveID;
    fs := schofer_tab`BorcherdsForms;
    keys_fs := schofer_tab`Keys_fs;
    row_scales := schofer_tab`RowScales;
    all_fs := [fs[k] : k in keys_fs];
    UpdateFieldsOfDefn(schofer_tab, dnew);
    flds := (schofer_tab`FldsOfDefn)[curveid][dnew];
    assert #flds eq 1;
    deg := Degree(flds[1]);
    if d in ds[1] then
        d_idx := Index(ds[1],d);
        ds[1][d_idx] := dnew;
    else
        qidx := Index(ds[2], d);
        d_idx := qidx+#ds[1];
        ds[2][qidx] := dnew;
    end if;
    Ldata := ShimuraCurveLattice(Xstar`D,Xstar`N);
    norm_val := AbsoluteValuesAtRationalCMPoint(all_fs, dnew, Xstar, Ldata);
    for i->v in norm_val do
        // table[i][d_idx] := norm_val[i]/row_scales[i]^deg;
        table[i][d_idx] := norm_val[i]-deg*row_scales[i];
        if not is_log then
            table[i][d_idx] := RationalNumber(table[i][d_idx]);
        end if;
    end for;
    schofer_tab`Values := table;
    schofer_tab`Discs := [ds[1], ds[2]];
    curves := schofer_tab`Curves;
    UpdateFieldsOfDefn(schofer_tab, dnew);  
    return;
end procedure;


function find_y2_scales(schofer_table)
    ds := schofer_table`Discs;
    ratds := ds[1];
    table := schofer_table`Values;
    keys_fs := schofer_table`Keys_fs;
    k_idxs := schofer_table`K_idxs;
    fldsofdef := schofer_table`FldsOfDefn;

    //Scale the y2 rows of the table

    scale_factors :=[];
    for i in k_idxs do
        if exists(j1){j : j->d1 in ratds  | #fldsofdef[keys_fs[i]][d1] eq 1 and Degree(fldsofdef[keys_fs[i]][d1][1]) eq 1 and table[i][j] ne LogSum(Infinity()) and table[i][j] ne LogSum(0)} then
            //then we have a rational point on X
            // d1 := ratds[j1];
            v1 := table[i][j1];
            // scale, _ := SquareFree(v1);
            log_scale := SquareFree(v1);
            // Append(~scale_factors, AbsoluteValue(scale));
            Append(~scale_factors, log_scale);
        else 
            assert exists(j1){j : j->d1 in ratds  | #fldsofdef[keys_fs[i]][d1] le 2 and {Degree(fldsofdef[keys_fs[i]][d1][k]) : k in [1..#fldsofdef[keys_fs[i]][d1]]} subset {1,2} and table[i][j] ne LogSum(Infinity()) and table[i][j] ne LogSum(0)};
            assert exists(j2){j : j->d2 in ratds  | #fldsofdef[keys_fs[i]][d2] le 2 and {Degree(fldsofdef[keys_fs[i]][d2][k]) : k in [1..#fldsofdef[keys_fs[i]][d2]]} subset {1,2}  and table[i][j] ne LogSum(Infinity()) and ratds[j1] ne d2 and table[i][j] ne LogSum(0)};
            //otherwise we find two points that are potentially over quadratic fields
            v1 := table[i][j1];
            v2 := table[i][j2];
            d1 := ratds[j1];
            d2 := ratds[j2];
            quad_idx_1 := 1;
            quad_idx_2 := 1;
            if fldsofdef[keys_fs[i]][d1][1] eq Rationals() then
                quad_idx_1 := 2;
            end if;
            if fldsofdef[keys_fs[i]][d2][1] eq Rationals() then
                quad_idx_2 := 2;
            end if;
            d1 := Discriminant(MaximalOrder(fldsofdef[keys_fs[i]][d1][quad_idx_1]));
            d2 := Discriminant(MaximalOrder(fldsofdef[keys_fs[i]][d2][quad_idx_2]));
            log_scale1 := SquareFree(v1 - LogSum(AbsoluteValue(d1)));
            log_scale2 := SquareFree(v1);
            // scale1, _ := SquareFree(v1/d1); //two possibilities
            // scale2, _ := SquareFree(v1);
            // if IsSquare(scale1*v2/d2) then
            if IsSquare(log_scale1 + v2 - LogSum(AbsoluteValue(d2))) then
                // Append(~scale_factors, AbsoluteValue(scale1));
                Append(~scale_factors, log_scale1);
            else
                // assert IsSquare(scale2*v2);
                assert IsSquare(log_scale2 + v2);
                // Append(~scale_factors, AbsoluteValue(scale2));
                Append(~scale_factors, log_scale2);
            end if;
        end if;
    end for;
    return scale_factors;

end function;

function find_y2_signs(table, keys_fs, curves, d, j, flds, k_idxs)
    //find signs of y^2 for rational CM point d on each y^2
    //in keys_fs, where j is index of column of d in table
    for k->i in k_idxs do
        if table[i][j] eq Infinity() then continue; end if;
        if table[i][j] eq 0 then continue; end if;
        fields := flds[keys_fs[i]][d];
        Fs_eps := [* <F, eps> : F in flds, eps in [-1,1] *];
        possible_answers := [* *];
        for eps in [-1,1] do
            y2 := eps*table[i][j];
            for F in fields do
                is_sqr, y := IsSquare(F!y2);
                if (is_sqr) then
                    if (Type(F) eq FldRat) or (Degree(F) eq Degree(sub<F|y>)) then
                        Append(~possible_answers, <F,eps,y>);
                    end if;
                end if;
            end for;
        end for;
        assert #possible_answers eq 1;
        eps := possible_answers[1][2];
        table[i][j] :=  eps * table[i][j];
    end for;
    return table;
end function;

intrinsic ValuesAtCMPoints(abs_schofer_tab::SchoferTable, all_cm_pts::SeqEnum) -> SchoferTable
    {}
    ds := abs_schofer_tab`Discs;
    table := abs_schofer_tab`Values;
    keys_fs := abs_schofer_tab`Keys_fs;
    k_idxs := abs_schofer_tab`K_idxs;
    Xstar := abs_schofer_tab`Xstar;
    row_scales := abs_schofer_tab`RowScales;
    curves := abs_schofer_tab`Curves;
    cid := Xstar`CurveID; 
    ratds := ds[1];
    quadds := ds[2];
    allds := ratds cat quadds;

    s_idx := abs_schofer_tab`sIndex;
    stilde_idx := abs_schofer_tab`sTildeIndex;
    s := table[s_idx];
    stilde := table[stilde_idx];


    table :=[* [* x : i->x in t *] : t in table *];
    scale_tilde := stilde[Index(s,LogSum(0))];
    scale := s[Index(stilde,LogSum(0))];
   
    row_scales[s_idx] +:= scale;
    row_scales[stilde_idx] +:= scale_tilde;
   
    k_idxs := abs_schofer_tab`K_idxs;
    abs_schofer_tab`Values := table;

    //Scale the y2 rows of the table
    scale_factors := find_y2_scales(abs_schofer_tab);

    // make table values into rational numbers

    degs := [1 : i in ds[1] ] cat [ 2 : i in ds[2]];
    for i->k in k_idxs do
        for j in [1 .. #table[k]] do
            // table[k][j] := table[k][j]/scale_factors[i]^degs[j];
            table[k][j] -:= degs[j]*scale_factors[i];
        end for;
        // row_scales[k] := row_scales[k]*scale_factors[i];
        row_scales[k] +:= scale_factors[i];
    end for;
    // abs_schofer_tab`Values := table;
    abs_schofer_tab`Values := [[RationalNumber(x) : x in y] : y in table];
    abs_schofer_tab`RowScales := row_scales;
    
    table := abs_schofer_tab`Values;
    table :=[* [* x : i->x in t *] : t in table *];

    s := table[s_idx];
    stilde := table[stilde_idx];
    s, stilde := find_signs(s, stilde, ds);

    table[s_idx] := s;
    table[stilde_idx] := stilde;

    //Next need to go from norms to values on the hauptmoduls for the quad_cm points

    all_flds := abs_schofer_tab`FldsOfDefn;

    bad_ds := {};
    for i->d in quadds do
        d_idx := #ratds+i;
        good_inds := [];
        currd := d;
        while #good_inds ne 1 do
            all_flds := abs_schofer_tab`FldsOfDefn;
            norm_s := table[s_idx][d_idx];
            norm_stilde := table[stilde_idx][d_idx];
            flds := all_flds[cid][currd];
            assert #flds eq 1;
            assert Degree(flds[1]) eq 2;
            K := flds[1]; //assume fields of definition are exactly quadratic on Xstar
            _<x> := PolynomialRing(Rationals());
            signs := [[1,1], [1,-1],[-1,1],[-1,-1]];
            minpolys := [];
            for eps in signs do
                    trace := 1 - eps[1]*norm_stilde +  eps[2]*norm_s;
                    Append(~minpolys, x^2 - trace*x + eps[2]*norm_s);
            end for;
            roots := [Roots(p,K) : p in minpolys];
            good_inds := [i : i->r in roots | #r ne 0 and not(&and[rt[1] in Rationals() : rt in r])];
            if #good_inds ne 1 then
                vprintf ShimuraQuotients, 1: "We need that there is a unique minpoly left after filtering by roots so we are replacing %o.\n", currd;
                Include(~bad_ds, currd);
                candidates := Set([pt[1] : pt in all_cm_pts[2]]) diff Set(quadds) diff bad_ds;
                require #candidates ge 1: "No possible choices of CM points left which we can pin down the correct minpoly";
                newd := Reverse(Sort(SetToSequence(candidates)))[1];
                replace_column(abs_schofer_tab, currd, newd, false);
                currd := newd;
                table := abs_schofer_tab`Values;
                table :=[* [* x : i->x in t *] : t in table *];
            end if;
        end while;
        table[s_idx][d_idx] := minpolys[good_inds[1]];
        norm_s := Coefficient(minpolys[good_inds[1]], 0);
        trace_s := - Coefficient(minpolys[good_inds[1]], 1);
        table[stilde_idx][d_idx] := x^2 - (2- trace_s)*x + (1- trace_s + norm_s);
        abs_schofer_tab`Values := table;
    end for;
       
    //Find signs on the y2 rows

    //Note that when we are at a quadratic CM point, K the quadratic field the norm is always positive
    //Write v = | y^2(tau)y^2(taubar)|. Then y(tau)y(taubar) = sqrt(eps*norm). Since the fields Q(y(tau)) = Q(y(tau))^sigma
    // where sigma is the unique nontrivial element of Gal(K/Q)
    // y(tau)y(taubar) lies in the fixed field by sigma, i.e. in Q
    // so eps* v is a square in Q, and is positive

    for j->d in ratds do
        table := find_y2_signs(table, keys_fs, curves, d, j, all_flds, k_idxs);
    end for;

    schofer_table := CreateSchoferTable(table, keys_fs, abs_schofer_tab`Discs, curves, Xstar);
    return schofer_table;
end intrinsic;

intrinsic ReduceTable(schofer_tab::SchoferTable)
    {}
    table := schofer_tab`Values;
    sep_ds := schofer_tab`Discs;
    scales := [];
    num_rat_ds := #sep_ds[1];
    for t in table do
        // xs := [x : i->x in t | i le num_rat_ds and x notin [0, Infinity()] ];
        xs := [x : i->x in t | i le num_rat_ds and x notin [LogSum(0), LogSum(Infinity())] ];
        //xs := [x : x in t | x notin [0, Infinity()]];
        /*
        ps := &join[Set(PrimeDivisors(Numerator(x))) : x in xs];
        ps join:= &join[Set(PrimeDivisors(Denominator(x))) : x in xs];
        ps := [p : p in ps];
        */
        ps := [p : p in &join[Keys(x`log_coeffs) : x in xs]];
        // vals := [[Valuation(x,p) : x in xs ] : p in ps];
        vals := [[IsDefined(x`log_coeffs, p) select x`log_coeffs[p] else 0 : x in xs ] : p in ps];
        mins := [Minimum([<AbsoluteValue(v),v> : v in valp]) : valp in vals];
        // scale := &*[Rationals() | p^mins[i][2] : i->p in ps];
        scale := &+([LogSum()] cat [LogSum(mins[i][2], p) : i->p in ps]);
        Append(~scales, scale);
    end for;
    degs := [1 : i in sep_ds[1] ] cat [ 2 : i in sep_ds[2]];
    // schofer_tab`Values :=  [[x/scales[i]^degs[j] : j->x in t] : i->t in table ];
    schofer_tab`Values :=  [[x - degs[j]*scales[i] : j->x in t] : i->t in table ];
    schofer_tab`RowScales := scales;
    return;
end intrinsic;

intrinsic ValuesAtCMPoints(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot] : MaxNum := 7, Prec := 100, Exclude := {}, Include := {}) -> SeqEnum, SeqEnum, SeqEnum
{Returns the values of y^2 for all degree 2 covers and two hauptmodules at CM points.}
    fs := BorcherdsForms(Xstar, curves : Prec := Prec);
    d_divs := &cat[[T[1]: T in  DivisorOfBorcherdsForm(f, Xstar)] : f in [fs[-1], fs[-2]]]; //include zero infinity of hauptmoduls
    all_cm_pts := CandidateDiscriminants(Xstar, curves);
    abs_schofer_tab, all_cm_pts := AbsoluteValuesAtCMPoints(Xstar, curves, all_cm_pts, fs : MaxNum := MaxNum, Prec := Prec, Exclude := {}, Include := Set(d_divs));
    ReduceTable(abs_schofer_tab);
    schofer_tab := ValuesAtCMPoints(abs_schofer_tab, all_cm_pts);
    return schofer_tab;
end intrinsic;

/* Bibliography
[GR] - Gonzales, Rotger - non-elliptic shimura curves of genus one
*/