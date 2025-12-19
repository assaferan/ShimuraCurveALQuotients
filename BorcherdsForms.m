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
        if not IsEmpty(qexps) then
            _<q> := Universe(qexps);
            min_v := Minimum([Valuation(f) : f in qexps]);
            coeffs := Matrix(Rationals(), [AbsEltseq(q^(-min_v)*f : FixedLength) : f in qexps]);
        else
            min_v := 0;
            coeffs := MatrixAlgebra(Rationals(), 0)!0;
        end if;
        
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
   
    if IsOdd(Xstar`D) then
        E0, nE0, _, eta_quotients_oo, eta_quotients_0 := WeaklyHolomorphicBasis(Xstar`D, Xstar`N : Prec := Prec, Zero, n0 := n0);
    end if;
    // we do this twice -- we should remember this
    pts, _ := RationalandQuadraticCMPoints(Xstar); // pts <-> infty, 0, rational
    pts := [p : p in pts | p[1] notin Exclude and GCD(p[1], Xstar`N) eq 1];
    require #pts ge 3 : "Could not find enough rational CM points!";


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
        if IsOdd(Xstar`D) then
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

                    if IsOdd(Xstar`D) then
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

                    if IsOdd(Xstar`D) then
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
                    if IsOdd(Xstar`D) then
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

