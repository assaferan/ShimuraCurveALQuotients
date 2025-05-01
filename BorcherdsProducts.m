// Compute the order of vanishing of eta(delta tau) at the cusp a/b in Gamma_0(N)
// multiplied by 24
function OrderOfVanishingOfEta(delta, b, N)
    return N*GCD(b,delta)^2 div (GCD(N,b^2)*delta);
end function;

/*
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
    fm1 := q^(-1)*Parent(t_eta_quotient)!Eltseq(E[1]);
    f0 := q^(-1)*Parent(t_eta_quotient)!Eltseq(E[2]);
    t := t_eta_quotient;
    fm6 := f0*(t^6 - 8*t^5+17*t^4+4*t^3-37*t^2+22*t);
    fm3 := f0*(t^3-5*t^2+5*t);
    fm2 := f0*(t^2-4*t);    
    psi_s := -2*fm6 + 2*fm1 + 6*f0;
    psi_y := -4*fm6 + 2*fm3 + 2*fm1 + 10*f0;
    B<i,j> := QuaternionAlgebra(Rationals(),-1,3);
    IsSquare(B!(-3));
    (i*j)^2;
    (i+i*j)^2;
    (j+i*j)^2;
    (j-i*j)^2;
    (3*ii*j)^2;
    (3*i+i*j)^2;
    (2*i+i*j)^2;
    (2*i+j)^2;
    (i+2*j)^2;
    psi_s;
    psi_s + 6*f0;
    psi_y;
    psi_y + 10*f0;
    B;
    Units(MaximalOrder(B));
    k := (1 + i + j + i*j)/2
    ;
    k^2;
    k^2-k-1;
    (2*k-1)^2;
    (j + 2*k-1)^2;
    (j - 2*k-1)^2;
    (1+i)^2;
    (i+j)^2;
    (i-j)^2;
    (i-i*j)^2;
    (i+i*j)^2;
    [<a,b,c,d> : a,b,c,d in [-2..2] | (a+b*i+c*j+d*k)^2 eq -3];
    -1+2*i+2*k;
    lambda := -1+2*i+2*k;
    lambda^2;
    Trace(lambda);
    Basis(B);
    [Trace(x*lambda) : x in Basis(B)];
    [Trace(x*lambda) : x in Basis(B)[2..4]];
    Matrix([[Trace(x*lambda) : x in Basis(B)[2..4]]]);
    Kernel(Transpose(Matrix([[Trace(x*lambda) : x in Basis(B)[2..4]]])));
    U := Kernel(Transpose(Matrix([[Trace(x*lambda) : x in Basis(B)[2..4]]])));
    lambda;
    O := MaximalOrder(B);
    lambda/2 in O;
    plus_minus := Matrix([[2,0,1],[2,-1,-1],[-1,0,1]]);
    Determinant(plus_minus);
    basis_O := Basis(O);
    basis_O;
    [Trace(x) : x in basis_O];
    Kernel(Transpose(Matrix([[Trace(x) : x in basis_O]])));
    L_space := Kernel(Transpose(Matrix([[Trace(x) : x in basis_O]])));
    Basis(L_space);
    [&+[b[i]*basis_O[i] : i in [1..4]] : b in Basis(L_space)];
    basis_L := [&+[b[i]*basis_O[i] : i in [1..4]] : b in Basis(L_space)];
    lambda;
    BM_L := Matrix([Eltseq(b) : b in basis_L]);
    BM_L;
    Solution(BM_L, Vector(Eltseq(lambda)));
    lambda_L := Solution(BM_L, Vector(Eltseq(lambda)));
    [Trace(x*lambda) : x in basis_L];
    Transpose(Kernel([Trace(x*lambda) : x in basis_L]));
    Kernel(Transpose(Matrix([[Trace(x*lambda) : x in basis_L]])));
    Matrix(Basis(Kernel(Transpose(Matrix([[Trace(x*lambda) : x in basis_L]])))));
    BLplus := Matrix(Basis(Kernel(Transpose(Matrix([[Trace(x*lambda) : x in 
    basis_L]])))));
    BLplus;
    BLplus * BM_L;
    L;
    BM_L;
    Q := Matrix([[Norm(x+y)^2-Norm(x)-Norm(y) : x in basis_L] : y in basis_L]);
    Q;
    Lat := LatticeWithGram(Q);
    LatticeWithGram;
    Lat := LatticeWithGram(Q : CheckPositive := false);
    Lat;
    Dual(Lat);
    Q;
    Q^(-1);
    L;
    basis_L;
    Q_O := Matrix([[Norm(x+y)^2-Norm(x)-Norm(y) : x in basis_O] : y in basis_O]);
    Q_O;
    Determinant(Q_O);
    basis_O;
    Discriminant(O);
    Q_O := Matrix([[Norm(x+y)-Norm(x)-Norm(y) : x in basis_O] : y in basis_O]);
    Q_O;
    Determinant(Q_O);
    Q := Matrix([[Norm(x+y)-Norm(x)-Norm(y) : x in basis_L] : y in basis_L]);
    Determinant(Q);
    Q;
    Q^(-1);
    Q^(-1)*BM_L;
    BM_Ldual := Q^(-1)*BM_L;
    eta := B!Eltseq(BM_Ldual[1]);
    eta;
    BLplus;
    BLplus * BM_L;
    Q;
    BLplus * Q;
    Kernel(BLplus * Q);
    Kernel(Transpose(BLplus * Q));
    lambda;
    Solution(BM_L, Vector(Eltseq(lambda)));
    eta;
    BM_Ldual[1];
    Solution(BM_L, BM_Ldual[1]);
    eta_L := Solution(BM_L, BM_Ldual[1]);
    eta_L*Transpose(BLplus*Q);
    eta_L;
    eta eq 1/6*lambda;
    eta := B!Eltseq(BM_Ldual[2]);
    eta_L := Solution(BM_L, BM_Ldual[2]);
    eta_L*Transpose(BLplus*Q);
    eta_U := eta_L*Transpose(BLplus*Q);
    eta_U*BLplus;
    eta - eta_U*BLplus;
    eta_L - eta_U*BLplus;
    eta_L*Q;
    (eta_L*Q, lambda_L);
    eta_lambda := ((eta_L*Q, lambda_L)/(lambda*Q,lambda))*lambda;
    eta_lambda := ((eta_L*Q, lambda_L)/(lambda_L*Q,lambda_L))*lambda_L;
    eta_lambda;
    eta - eta_lambda;
    eta_L - eta_lambda;
    eta_U := eta_L - eta_lambda;
    eta_plus := eta_lambda;
    eta_minus := eta_U;
    Solution(BM_L, eltseq(i));
    Solution(BM_L, Eltseq(i));
    Solution(BM_L, Vector(Eltseq(i)));
    mu := Solution(BM_L, Vector(Eltseq(i)));
    mu_plus:= ((mu*Q, lambda_L)/(lambda_L*Q,lambda_L))*lambda_L;
    mu_plus;
    eta_plus;
    mu - mu_plus;
    mu_minus := mu - mu_plus;
    eta_plus + mu_plus;
    lambda_L;
    (lambda_L*Q, lambda);
    (lambda_L*Q, lambda_L);
    SmithForm(Q^(-1));
    Q^(-1);
    BM_Ldual;
    Image(BM_Ldual);
    Module(Integers(),BM_Ldual);
    RModule(Integers(),BM_Ldual);
    RSpaceWithBasis;
    RSpaceWithBasis(BM_Ldual);
    tmpdual := RSpaceWithBasis(BM_Ldual);
    tmp := RSpaceWithBasis(BM_L);
    tmp;
    tmpdual / tmp;
    RSpaceWithBasis(Integers(), BM_Ldual);
    RSpaceWithBasis(ChangeRing(6*BM_Ldual,Integers()));
    tmpdual := RSpaceWithBasis(ChangeRing(6*BM_Ldual,Integers()));
    tmp := RSpaceWithBasis(ChangeRing(6*BM_L,Integers()));
    tmpdual / tmp;
    disc_grp, disc_quo := tmpdual / tmp;
    #disc_grp;
    Gamma0(12);
    CosetRepresentatives(Gamma0(12));
    gammas := CosetRepresentatives(Gamma0(12));
    gammas := gammas[1];
    gammas := CosetRepresentatives(Gamma0(12));
    gamma := gammas[1];
    gamma;
    gamma := gammas[2];
    gamma;
    gamma^(-1);
    WordProblem;
    WordProblem(gamma^(-1));
    WordProblem;
    ListSignatures(GrpPSL2Elt);
    FindWord(PSL2(Integers()), gamma^(-1));
    FindWord;
    G := PSL2(Integers());
    gens := Generators(G);
    gens;
    gens[1]*gens[2];
    gamma^(-1);
    word := FindWord(PSL2(Integers()), gamma^(-1));
    &*[gens[Abs(word[i])]^Sign(word[i]) : i in [1..#word]] eq gamma^(-1);
    gens[1];
    gens[2];
    psi_s;
    psi := psi_s;
    return 0;
end function;
*/

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
    return eqs, ieqs, t;
end function;

function write_polymake_scriptfile(D,N)
    eqs, ieqs, t := integer_programming_input(D,N);
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
    return t;
end function;

function get_integer_prog_solutions(D,N)
    t := write_polymake_scriptfile(D,N);
    fname := Sprintf("polymake_script_%o_%o", D, N);
    polymake := Read(POpen("polymake --script " cat fname, "r"));
    sol_lines := Split(polymake, "\n");
    sol_vecs := [Split(line, " ") : line in sol_lines];
    sols := [[eval(x) : x in vec] : vec in sol_vecs];
    M := 2*D*N;
    rs := [sol[2..1 + #Divisors(M)] : sol in sols];
    return rs, Eltseq(t)[1..#Divisors(M)];
end function;

function get_weakly_holomorphic_basis(D,N)
    D0,M,g := get_D0_M_g(D,N);
    _<q> := PuiseuxSeriesRing(Rationals());
    eta := DedekindEta(q);  
    nor_eta := eta / q^(1/24);
    rs, t := get_integer_prog_solutions(D,N);
    M := 2*D*N;
    eta_quotients := [&*[(Evaluate(nor_eta,q^d)*q^(d/24))^r[i] : 
                        i->d in Divisors(M)] : r in rs];
    t_eta_quotient := &*[(Evaluate(nor_eta,q^d)*q^(d/24))^t[i] : 
                        i->d in Divisors(M)];
    min_v := Minimum([Valuation(eta_quot) : eta_quot in eta_quotients]);
    min_prec := Minimum([RelativePrecision(eta_quot) - min_v + Valuation(eta_quot) : eta_quot in eta_quotients]);
    R<q> := LaurentSeriesRing(Rationals());
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

function FourthPowerFree(a)
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
    free, rad := FourthPowerFree(d/a);
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
    free, rad := FourthPowerFree(ret_c4);
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
    // free, rad := FourthPowerFree(ret_c4);
    // return rad^(-1)*ret_ceta, free;
    return ret_eta;
end function;

function ShimuraCurveLattice(D,N)
    B := QuaternionAlgebra(D);
    O_max := MaximalOrder(B);
    O := Order(O_max,N);
    basis_O := Basis(O);
    L_space := Kernel(Transpose(Matrix([[Trace(x) : x in basis_O]])));
    basis_L := [&+[b[i]*basis_O[i] : i in [1..4]] : b in Basis(L_space)];
    BM_L := Matrix([Eltseq(b) : b in basis_L]);
    Q := Matrix([[Norm(x+y)-Norm(x)-Norm(y) : y in basis_L] : x in basis_L]);
    BM_Ldual := Q^(-1)*BM_L;
    // L := LatticeWithGram(Q : CheckPositive := false);
    // return L;
    denom := Denominator(BM_Ldual);
    // We are modifying it to be always with respect to the basis of L.
    // Ldual := RSpaceWithBasis(ChangeRing(denom*BM_Ldual,Integers()));
    Ldual := RSpaceWithBasis(ChangeRing(denom*Q^(-1), Integers()));
    // L := RSpaceWithBasis(ChangeRing(denom*BM_L,Integers()));
    L := RSpaceWithBasis(ScalarMatrix(3,denom));
    disc_grp, to_disc := Ldual / L;
    return L, Ldual, disc_grp, to_disc, Q^(-1);
end function;

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
function SpreadBorcherds(alphas, rs, ds, Ldual, discL, Qdisc, to_disc)
    M := #discL;
    G := PSL2(Integers());
    gens := Generators(G);
    gammas := CosetRepresentatives(Gamma0(M));
    T := gens[1];
    S := gens[2];
    assert Eltseq(T) eq [1,1,0,1];
    assert Eltseq(S) eq [0,1,-1,0];
    F := [* PuiseuxSeriesRing(Rationals())!0 : i in [1..M] *];
    e0 := [1] cat [0 : i in [1..M-1]];
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
end function;

function FindLambda(Q, d : bound := 10)
    Q := ChangeRing(Q, Integers());
    n := Nrows(Q);
    idxs := CartesianPower([-bound..bound], n);
    for idx in idxs do
        v := Vector([idx[j] : j in [1..n]]);
        v := ChangeRing(v, BaseRing(Q));
        if (v*Q,v) eq 2*d then
            return true, v;
        end if;
    end for;
    return false, _;
end function;

function VerticalJoinList(mats)
    m := mats[1];
    for i in [2..#mats] do
        m := VerticalJoin(m, mats[i]);
    end for;
    return m;
end function;

function MyLegendreSymbol(alpha, p)
    return LegendreSymbol(Integers()!(GF(p)!alpha),p);
end function;

// functions for testing the Wpolys from Kudla Yang paper
function bp_Kudla_Yang_poly(p, kappam, D)
    d, c := SquarefreeFactorization(4*kappam);
    k := Valuation(c,p);
    F := Rationals();
    R<x> := PolynomialRing(F);
    vp := LegendreSymbol(d, p);
    if k lt 0 then return R!1; end if;
    // k >= 0
    if (D mod p ne 0) then
        return (1 - vp*x + p^k*vp*x^(1+2*k)-p^(k+1)*x^(2*k+2))/(1-p*x^2);
    end if;
    // p divides D
    return ((1-vp*x)*(1-p^2*x^2)-vp*p^(k+1)*x^(2*k+1)+p^(k+2)*x^(2*k+2)+vp*p^(k+1)*x^(2*k+3)-p^(2*k+2)*x^(2*k+4))/(1-p*x^2);
end function;

function sigmasp_Kudla_Yang_poly(p, m, kappa, is_even)
    F := Rationals();
    R<x> := FunctionField(F);
    assert IsSquarefree(kappa);
    chi_p := is_even select KroneckerCharacter(kappa)(p) else KroneckerCharacter(2*kappa)(p);
    if m eq 0 then
        return (1 - chi_p*x)^(-1);
    end if;
    return &+[(chi_p*x)^r : r in [0..Valuation(m,p)]];
end function;

procedure test_kronecker_sigma(B)
    kappas := [kappa : kappa in [1..B] | IsSquarefree(kappa)];
    for p in PrimesUpTo(B) do
        for kappa in kappas do
            assert sigmasp_Kudla_Yang_poly(p,0,kappa,true)*EulerFactor(KroneckerCharacter(kappa),p) eq 1;
        end for;
    end for;
end procedure;

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
    if not IsIntegral(m - Q_mu) then
        return 0;
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
    eps_mu := func<k | LegendreSymbol(-1,p)^(l_mu(k) div 2) * &*[Integers() | MyLegendreSymbol(eps[i],p) : i in L_mu(k)]>;
    f_1 := function(x)
        a, alpha := Valuation(x,p);
        return IsEven(l_mu(a+1)) select -1/p else MyLegendreSymbol(alpha,p) / sqrtp;
    end function;
    t_mu := m - &+[Rationals() | eps[i]*p^l[i]*mu_wrt_B[i]^2 : i in [1..n] | i notin H_mu];
    a := Valuation(t_mu, p);
    R<x> := PolynomialRing(K);
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
        for j in [2..Nrows(Jblock)] do
            row_ind +:= 1;
            b := Jblock[j-1,j];
            if b eq 0 then
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
            a := Jblock[j-1,j-1];
            d := Jblock[j,j];
            disc := b^2 - a*d;
            if (Integers(8)!disc eq 5) then 
                disc -:= 2*a; 
                Append(~n_list, exps[i]);
                Append(~mu_prime_prime_indices, row_ind);
            else   
                Append(~m_list, exps[i]);
                Append(~mu_prime_indices, row_ind);
            end if;
            is_sqr, sqrt_disc := IsSquare(Zp!disc);
            assert is_sqr;
            // solving the quadratic
            x1 := (-b + sqrt_disc)/a;
            x2 := (-b - sqrt_disc)/a;
            z2 := (-a)/(2*disc); // constant to get scalar product equal to 1
            // Change of basis matrix
            B := Matrix([[x1, 1], [z2*x2, z2]]);
            cans := [SymmetricMatrix([0,1,0]), SymmetricMatrix([2,1,2])];
            assert B*Matrix([[a,b],[b,d]])*Transpose(B) in cans;
            B_big := ChangeRing(Parent(Jblock)!1, Zp);
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
    assert &and[Valuation(e,p) eq 0 : e in eps];
    B := ChangeRing(VerticalJoinList(bases), Rationals());
    mu_wrt_L := Solution(ChangeRing(BasisMatrix(L),Rationals()), ChangeRing(mu, Rationals()));
    Q_mu := 1/2*(mu_wrt_L * ChangeRing(S,Rationals()), mu_wrt_L);
    if not IsIntegral(m - Q_mu) then
        return 0;
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
    d2_mu := func<k | 2*k + &+[Minimum(l_list[i]-k, 0) : i in H_mu + 2*&+[Minimum(m_list[i]-k,0) : i in M_mu] + 2*&+[Minimum(n_list[i]-k,0) : i in N_mu]]>;

    vals := [l_list[i] + Valuation(mu_list[i], p) : i in [1..H] | i notin H_mu and Valuation(mu_list[i],p) lt -1];
    vals cat:= [l_list[i] : i in [1..H] | i notin H_mu and Valuation(mu_list[i],p) eq -1];
    vals cat:= [m_list[i] + Minimum(Valuation(mu_prime_list[i][1], p), Valuation(mu_prime_list[i][2], p)) : i in [1..M] | i notin M_mu];
    vals cat:= [n_list[i] + Minimum(Valuation(mu_prime_prime_list[i][1], p), Valuation(mu_prime_prime_list[i][2], p)) : i in [1..N] | i notin N_mu];
    if IsEmpty(vals) then
        K_0 := Infinity();
    else
        K_0 := Minimum(vals);
    end if;
    
    p_mu := func<k | (-1)^(&+[Minimum(n_list[j] - k, 0) : j in N_mu])>;
    eps_mu := func<k | &*[eps[h] : h in L_mu(k)]>;

    delta_mu := func<k | exists(h){h : h in H_mu | l_list[h] eq k} select 0 else 1>;

    two_block := func< x | x[1]^2 + x[1]*x[2] + x[2]^2>;

    Q_prime_mu := &+[Rationals() | eps[i]*p^l_list[i]*mu_list[i]^2 : i in [1..H] | i notin H_mu];
    Q_prime_mu +:= &+[Rationals() | eps_prime[i]*p^m_list[i]*(&* mu_prime_list[i]) : i in [1..M] | i notin M_mu];
    Q_prime_mu +:= &+[Rationals() | eps_prime_prime[i]*p^n_list[i]*two_block(mu_prime_prime_list[i]) : i in [1..N] | i notin N_mu];

    t_mu := m - Q_prime_mu;
    a := Valuation(t_mu, p);

    nu := func< k | t_mu*p^(3-k) - &+[eps[h] : h in H_mu | l_list[h] lt k]>;

    psi_char := func<k | (nu(k) mod 4 eq 0) select (-1)^(nu / 4) else 0>;

    K_0 := Minimum(K_0, a+3);

    R<x> := PolynomialRing(K);
    ret := 1;
    ret +:= &+[R | delta_mu(k)*p_mu(k)*sqrtp^(d2_mu(k)-3)*KroneckerSymbol(2,eps_mu(k)*nu(k))*x^k : k in [1..K_0] | IsOdd(l_mu(k))];
    ret +:= &+[R | delta_mu(k)*p_mu(k)*sqrtp^(d2_mu(k)-2)*KroneckerSymbol(2,eps_mu(k))*psi_char(k)*x^k : k in [1..K_0] | IsEven(l_mu(k))];

    return ret;
end function;

// Prop. 2.1 in [KY] says that if chi_p is unramified and p is odd in the odd dimensional case
// we have Wmp(s) = sigma_{-s,p}(m,chi)/Lp(s+1,chi) in the even case
// and Lp(s+1/2,chi_{kappam})/zetap(2s+1)*bp(kappam,s+1/2) in the odd case
// In particular, when m = 0 this should yield
// Lp(s,chi) / Lp(s+1,chi) in the even case, and 
// zeta_p(2s) / zeta_p(2s+1) in the odd case

function W(m,p,mu,L,Q)
    n := Rank(L);
    s0 := n/2 - 1;
    s := -s0;
    s2 := 2 - n; // s2 = 2*s always integral
    S := BasisMatrix(L)*Q*Transpose(BasisMatrix(L));
    D := Determinant(S);
    vpD := Valuation(D,p);
    K<sqrtp> := QuadraticField(p);
    scale := sqrtp^(-vpD);
    euler := p eq 2 select Wpoly2(m,mu,L,K,Q) else Wpoly(m,p,mu,L,K,Q);
    return scale*Evaluate(euler, sqrtp^(-s2));
end function;

function kappaminus(mu, m, Lminus, Q, d)
    ret := 0;
    if not IsIntegral(m) then
        return 0;
    end if;
    Bminus := BasisMatrix(Lminus);
    Delta := Determinant(Bminus*Q*Transpose(Bminus));
    Sm_mu := {p : p in PrimeDivisors(Delta)} join {p : p in PrimeDivisors(m)};
    h := ClassNumber(d);
    w := #UnitGroup(QuadraticField(d));
    return ret;
end function;

// Computes kappa0(m) in Schofer's formula
function kappa0(m, d, Q)
    Q := ChangeRing(Q, Integers());
    _, lambda_v := FindLambda(Q,-d);
    // lambda := &+[lambda_v[i]*basis_L[i] : i in [1..#basis_L]];
    c_Lplus := Content(lambda_v);
    Lplus := RSpaceWithBasis(Matrix(lambda_v div c_Lplus));
    Lminus := Kernel(Transpose(Matrix(lambda_v*Q)));
    L := RSpaceWithBasis(IdentityMatrix(Integers(),3));
    L_quo, L_quo_map := L / (Lplus + Lminus);
    s := 0;
    for mu_bar in L_quo do
        mu := mu_bar@@L_quo_map;
        c_mu_plus := ((mu*Q, lambda_v)/(lambda_v*Q,lambda_v));
        mu_plus:= c_mu_plus*ChangeRing(lambda_v, Rationals());
        mu_minus := mu - mu_plus;
        // finding the possible range of x in mu_plus + L_plus
        // use that mu_plus = c_mu_plus * lambda, L_plus = Z * c_Lplus^(-1) * lambda
        // that <lambda,lambda> = -2d, and we only need x with <x,x> <= 2m
        // so if x = c_mu_plus + c_Lplus^(-1)*k, we need only those with
        // (c_mu_plus + c_Lplus^(-1)*k)^2 le m/(-d)
        // thus k is between the following bounds
        lb := Ceiling((m/d - c_mu_plus)*c_Lplus);
        ub := Floor((m/(-d) - c_mu_plus)*c_Lplus);
        for k in [lb..ub] do
            x := (c_mu_plus + k * c_Lplus^(-1)) * ChangeRing(lambda_v, Rationals());
            // s +:= kappaminus(mu_minus, m - (x*ChangeRing(Q,Rationals()),x)/2);
            print mu_plus, m - (x*ChangeRing(Q,Rationals()),x)/2;
        end for;
    end for;
    return s;
end function;