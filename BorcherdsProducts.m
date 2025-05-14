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

intrinsic GetWeaklyHolomorphicBasis(D::RngIntElt,N::RngIntElt) -> .
{returns a weakly holomorphic basis corresponding to D, N.}
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
    return E, n, t_eta_quotient;
end intrinsic;

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

intrinsic ShimuraCurveLattice(D::RngIntElt,N::RngIntElt) -> .
{return the lattice correpsonding to the Eichler order of level N in the Quaternion algebra of discriminant D.}
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
end intrinsic;

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
    d, c := SquarefreeFactorization(Integers()!(4*kappam));
    k := Valuation(c,p);
    F := Rationals();
    R<x> := PolynomialRing(F);
    vp := KroneckerSymbol(d, p);
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
            B_big := Parent(Jblock)!1;
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

    R<x> := PolynomialRing(K);
    ret := 1;
    ret +:= &+[R | delta_mu(k)*p_mu(k)*sqrtp^(d2_mu(k)-3)*KroneckerSymbol(2,Integers()!(eps_mu(k)*nu(k)))*x^k : k in [1..K_0] | IsOdd(l_mu(k))];
    ret +:= &+[R | delta_mu(k)*p_mu(k)*sqrtp^(d2_mu(k)-2)*KroneckerSymbol(2,Integers()!(eps_mu(k)))*psi_char(k)*x^k : k in [1..K_0] | IsEven(l_mu(k))];

    return ret;
end function;

function Wpoly_scaled(m,p,mu,L,Q)
    S := BasisMatrix(L)*Q*Transpose(BasisMatrix(L));
    D := Determinant(S);
    vpD := Valuation(D,p);
    K<sqrtp> := QuadraticField(p);
    scale := sqrtp^(-vpD);
    euler := p eq 2 select Wpoly2(m,mu,L,K,Q) else Wpoly(m,p,mu,L,K,Q);
    return scale*euler;
end function;


// Testing Proposition 5.1 in [KY]
// Should have Wp(s-1/2,m,mu) = Lp(s,chi_{kappa m})/zeta_p(2s) bp(kappa m, s) * (m - Q(mu) in Zp)
// Note that our Wpoly_scaled is evaluated at (s+s_0) and when n = 1, s0 = n/2 - 1 = -1/2
// We also have:
// zeta_p(2s)^(-1) = 1 - X^2
// Lp(s, chi_{kappa m}) = (1 - chi_{kappa m}(p) X)^(-1) 
// There is a sqrtp factor that I am missing

// procedure test_bp_KY(B)
function test_bp_KY(B)
    _<x> := PolynomialRing(Rationals());
    kappas := [kappa : kappa in [1..B] | IsSquarefree(kappa)];
    L := RSpaceWithBasis(IdentityMatrix(Integers(),1));
    failures := [* *];
    for m0 in [1..B] do
        for kappa in kappas do
            Q := Matrix([[2*kappa]]);
            Q_rat := ChangeRing(Q, Rationals());
            for P in PrimesUpTo(B, Rationals() : coprime_to := kappa) do
                mus := [Vector(Rationals(), [0])];
                p := Norm(P);
                if (p eq 2) then Append(~mus, Vector([1/2])); end if;
                for mu in mus do
                    m := m0 + 1/2*(mu*Q_rat, mu);
                    // kappa is NOT replaced by 2*kappa inthe character because the rank is odd - see [KY] (2.9)
                    rhs := ((1-x^2)/EulerFactor(KroneckerCharacter(Integers()!(kappa*Numerator(m)*Denominator(m))),p))*bp_Kudla_Yang_poly(p, kappa*m,1);
                    assert Denominator(rhs) eq 1;
                    rhs := Numerator(rhs);
                    K<sqrtp> := QuadraticField(p);
                    lhs := (p eq 2) select Wpoly2(m,mu,L,K,Q) else Wpoly(m,p,mu,L,K,Q);
                    // assert lhs eq ChangeRing(rhs, BaseRing(lhs));
                    if lhs ne ChangeRing(rhs, BaseRing(lhs)) then
                        Append(~failures, [* m0, kappa, p, mu *]);
                    end if;
                end for;
            end for;
        end for;
    end for;
    // return;
    return failures;
// end procedure;
end function;

// Prop. 2.1 in [KY] says that if chi_p is unramified and p is odd in the odd dimensional case
// we have Wmp(s) = sigma_{-s,p}(m,chi)/Lp(s+1,chi) in the even case
// and Lp(s+1/2,chi_{kappam})/zetap(2s+1)*bp(kappam,s+1/2) in the odd case
// In particular, when m = 0 this should yield
// Lp(s,chi) / Lp(s+1,chi) in the even case, and 
// zeta_p(2s) / zeta_p(2s+1) in the odd case


function W(m,p,mu,L,Q)
    Wpoly := Wpoly_scaled(m,p,mu,L,Q);
    _<sqrtp> := BaseRing(Wpoly);
    n := Rank(L);
    // s0 := n/2 - 1;
    // s := -s0;
    s2 := 2 - n; // s2 = 2*s always integral
    return Evaluate(Wpoly, sqrtp^(-s2));
end function;

procedure test_W()
    // testing the few values we know from Yang
    L, Ldual, disc_grp, to_disc, Qinv := ShimuraCurveLattice(6,1);
    Q := ChangeRing(Qinv^(-1), Integers());
    for d in [-3,-4] do
        _, lambda_v := FindLambda(Q,-d);
        Lminus := Kernel(Transpose(Matrix(lambda_v*Q)));
        mu := Vector([0,0,0]);
        if d eq -4 then
            w32<x> := Wpoly_scaled(3,2,mu,Lminus,Q);
            assert w32 eq 1/2*(1-x^2);
            w33<x> := Wpoly_scaled(3,3,mu,Lminus,Q);
            assert w33 eq 1/3*(1+2*x+x^2);
            w22<x> := Wpoly_scaled(2,2,mu,Lminus,Q);
            assert w22 eq 1/2*(1+x^3);
            w23<x> := Wpoly_scaled(2,3,mu,Lminus,Q);
            assert w23 eq 1/3*(1-x);
        end if;
        if d eq -3 then
            w12<x> := Wpoly_scaled(1,2,mu,Lminus,Q);
            assert w12 eq 1/2*(1-x);
            w13<x> := Wpoly_scaled(1,3,mu,Lminus,Q);
            _<sqrt3> := BaseRing(w13);
            assert w13 eq 1/sqrt3*(1+x);
        end if;
    end for;
    return;
end procedure;

// returns x,y such that the answer is x logy
function kappaminus(mu, m, Lminus, Q, d)
    if (m eq 0) and (mu ne 0) then
        error Error("Not implemented for m = mu = 0!");
        print "special case";
        return 0, 1;
    end if;
    error if m eq 0, "Not implemented for m eq 0!\n";  
    Bminus := BasisMatrix(Lminus);
    Delta := Determinant(Bminus*Q*Transpose(Bminus));
    
    Sm_mu := {p : p in PrimeDivisors(Delta)} join {p : p in PrimeDivisors(Numerator(m))};
    Sm_mu := [p : p in Sm_mu];
    
    vprintf ShimuraQuotients, 2: "Sm_mu := %o\n", Sm_mu;

    Wpolys := [* Wpoly_scaled(m,p,mu,Lminus,Q) : p in Sm_mu *];
    assert exists(i){i : i in [1..#Sm_mu] | Evaluate(Wpolys[i],1) eq 0};
    p_prime := Sm_mu[i];
    Wpol := Wpolys[i];

    F := BaseRing(Wpolys[1]);
    sqrtps := [* F.1 *];
    for j in [2..#Wpolys] do
        Fj := BaseRing(Wpolys[j]);
        Append(~sqrtps, Fj.1);
        composites := CompositeFields(F, Fj);
        // assert exists(K){K : K in composites | IsTotallyPositive(K!(Fj.1) / K!(F1.1))};
        F := composites[1];
    end for;

    W_prod := &*[F | Evaluate(Wpolys[j],1) : j in [1..#Sm_mu] | j ne i];
    W_prod *:= Evaluate(Derivative(Wpol),1); // this should be multiplied by -log(p_prime)

    kron_prod := &*[F | 1 - Evaluate(KroneckerCharacter(d),p)/p : p in Sm_mu];

    h := ClassNumber(d);
    w := #UnitGroup(QuadraticField(d));

    assert IsTotallyReal(F);

    //is_sqr, sqrtd := IsSquare(F!AbsoluteValue(d));
    Fz<z> := PolynomialRing(F);
    sqrtd :=  Roots(z^2 - AbsoluteValue(d))[1][1];
    
    assert exists(v){v : v in RealPlaces(F) | &and[Evaluate(F!sqrtp, v) gt 0 : sqrtp in sqrtps]};
    if Evaluate(sqrtd, v) lt 0 then sqrtd := -sqrtd; end if;
    assert Evaluate(sqrtd, v) gt 0;
    
    ret := -sqrtd*w*W_prod / (h*kron_prod);

    ret := Rationals()!ret;
    vprintf ShimuraQuotients, 2 : "adding %o log %o\n", -ret, p_prime;
    return -ret, p_prime; // to get xlogy instead of -xlogy
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
intrinsic Kappa0(m::RngIntElt, d::RngIntElt, Q::AlgMatElt) -> FldReElt
{Computing coefficients Kappa0(m) in Schofers formula}
    vprintf ShimuraQuotients, 1:"Kappa0 of %o\n", m;
    Q := ChangeRing(Q, Integers());
    bd := 10;
    found_lambda := false;
    while not found_lambda do
        bd *:= 2;
        found_lambda, lambda_v := FindLambda(Q,-d : bound := bd);
    end while;
    assert found_lambda;
    // lambda := &+[lambda_v[i]*basis_L[i] : i in [1..#basis_L]];
    c_Lplus := Content(lambda_v);
    Lplus := RSpaceWithBasis(Matrix(lambda_v div c_Lplus));
    Lminus := Kernel(Transpose(Matrix(lambda_v*Q)));
    L := RSpaceWithBasis(IdentityMatrix(Integers(),3));
    L_quo, L_quo_map := L / (Lplus + Lminus);
    log_coeffs := AssociativeArray();
    // ret := 1;
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
        sqr_bd := m/(-d);
        // lb := Ceiling((m/d - c_mu_plus)*c_Lplus);
        // ub := Floor((m/(-d) - c_mu_plus)*c_Lplus);
        lb := Ceiling((-Sqrt(sqr_bd) - c_mu_plus)*c_Lplus);
        ub := Floor((Sqrt(sqr_bd) - c_mu_plus)*c_Lplus);
        for k in [lb..ub] do
            x := (c_mu_plus + k * c_Lplus^(-1)) * ChangeRing(lambda_v, Rationals());
            a, p := kappaminus(mu_minus, m - (x*ChangeRing(Q,Rationals()),x)/2, Lminus, Q, d);
            assert (m - (x*ChangeRing(Q,Rationals()),x)/2) ge 0;
            // ret *:= kappaminus(mu_minus, m - (x*ChangeRing(Q,Rationals()),x)/2, Lminus, Q, d);
            if not IsDefined(log_coeffs, p) then log_coeffs[p] := 0; end if;
            log_coeffs[p] +:= Rationals()!a;
            vprintf ShimuraQuotients, 2: "mu_minus = %o, m - Q(x) = %o\n", mu_minus, m - (x*ChangeRing(Q,Rationals()),x)/2;
        end for;
    end for;
    return log_coeffs;
    // return ret;
end intrinsic;

function better_code_we_wrote()
    D := 6;
    N := 1;
    E, n, t_eta_quotient := GetWeaklyHolomorphicBasis(D,N);
    _<q> := Parent(t_eta_quotient);
    fm_n := q^(-n)*Parent(t_eta_quotient)!Eltseq(E[1]);


    t := t_eta_quotient;
    f0 := q^(-1)*Parent(t_eta_quotient)!Eltseq(E[2]);
    fm1 := fm_n;
    fm6 := f0*(t^6 - 8*t^5+17*t^4+4*t^3-37*t^2+22*t);
    fm3 := f0*(t^3-5*t^2+5*t);
    fm2 := f0*(t^2-4*t);
    // s is such that s(tau_{-4}) = 0, s(tau_{-3}) is rational and s(tau_{-24}) = infty  
    // P_{-24} is elliptic of order 2, P_{-4} is elliptic of order 4 and P_{-3} is elliptic of order 6.
    // div(s) = P_{-4} - P_{-24}, 
    // We choose psi_s such that ([GY, Lemma 22]) div(psi_s) = (-2)*1/2*P_{-24} + 4*1/4*P_{-4}
    psi_s := -2*fm6 + 4*fm1 + 6*f0;
    // s_tilde is such that s(tau_{-4}) is rational, s(tau_{-24}) is infinity and s(tau_{-3}) = 0
    // div(s_tilde) = P_{-3} - P_{-24}
    // We choose psi_s_tilde such that ([GY, Lemma 22]) div(psi_s) = (-2)*1/2*P_{-24} + 6*1/6*P_{-3}
    psi_s_tilde := -2*fm6 + 6*fm3;
    // We look for the div(y^2), where y^2 = bs(s-s(tau_{-3})) is a model for X/w_6
    // div(y^2) = div s + div (s - s(tau_{-3})) = P_{-3} + P_{-4} - 2P_{-24}
    // We choose psi_y2 such that ([GY, Lemma 22]) div(psi_y2) = (-4)*1/2*P_{-24} + 4*1/4*P_{-4} + 6*1/6*P_{-3}
    psi_y2_w6 := -4*fm6 + 6*fm3 + 4*fm1 + 6*f0;
    // We look for the div(y^2), where y^2 = c(s-s(tau_{-3})) is a model for X/w_2
    // div(y^2) = div (s - s(tau_{-3})) = P_{-3} - P_{-24}
    // We choose psi_y2 such that ([GY, Lemma 22]) div(psi_y2) = (-4)*1/2*P_{-24} + 4*1/4*P_{-4} + 6*1/6*P_{-3}
    L, Ldual, disc_grp, to_disc, Qinv := ShimuraCurveLattice(D,N);
    Q := ChangeRing(Qinv^(-1), Integers());
    // Q(sqrt(d)) has class number one, so there is a single pt P_{d}
    // There are 4 pts at the top curve X(6,1), but they are not fixed by any AL, so become 1 on the quotient
    // There is also a single point P_{-3}
    psis := [psi_s, psi_s_tilde, psi_y2_w6];
    psi_s_vals := AssociativeArray();
    psi_s_tilde_vals := AssociativeArray();
    psi_y2_vals := AssociativeArray();
    psi_s_vals[-4] := 0;
    psi_s_vals[-24] := Infinity();
    psi_s_tilde_vals[-3] := 0;
    psi_s_tilde_vals[-24] := Infinity();
    psi_y2_vals[-3] := 0;
    psi_y2_vals[-4] := 0;
    psi_s_vals[-24] := Infinity();
    psi_y2_vals[-24] := Infinity();
    psi_vals := [psi_s_vals, psi_s_tilde_vals, psi_y2_vals];
    // These should be rational CM points
    for d in [-3,-4,-19,-43,-67] do
        printf "Computing values of s, s_tilde, y^2 at d = %o...\n", d;
        for i->psi_val in psi_vals do
            if not IsDefined(psi_val, d) then
                schofer := SchoferFormula(psis[i],d,6,1);
                psi_vals[i][d] := &*[p^(Integers()!schofer[p]) : p in Keys(schofer)];
            end if;
        end for;
    end for;
    /*
    epsilon := 10^(-10);
    for d in [-3,-19,-43,-67] do
        n_d := NumberOfOptimalEmbeddings(MaximalOrder(QuadraticField(d)), D, N);
        assert n_d gt 0;
        if (d ne -3) then
            n_d := n_d div 4;
        end if;
        psi_s_at_d := (&*[Kappa0(m,d,Q)^(Integers()!Coefficient(psi_s,-m)) : m in [1..6]])^(-n_d/4);
        psi_y2_at_d := (&*[Kappa0(m,d,Q)^(Integers()!Coefficient(psi_y2,-m)) : m in [1..6]])^(-n_d/4);
        assert AbsoluteValue(psi_s_at_d - Round(psi_s_at_d)) lt epsilon;
        assert AbsoluteValue(psi_y2_at_d - Round(psi_y2_at_d)) lt epsilon;
        Append(~psi_s_vals, Round(psi_s_at_d));
        Append(~psi_y2_vals, Round(psi_y2_at_d));
    end for;
    */
    // psi_s_at_m3 := (&*[kappa0(m,-3,Q)^(Integers()!Coefficient(psi_s,-m)) : m in [1..n]])^(-n_m3/4); // 16
    // Q(sqrt(-19)) has class number one, so there is a single pt P_{-19}
    // There are 4 pts at the top curve X(6,1), but they are not fixed by any AL, so become 1 on the quotient
    // psi_s_at_m19 := (&*[kappa0(m,-19,Q)^(Integers()!Coefficient(psi_s,-m)) : m in [1..n]])^(-n_m19/4); // 1

    D := 26;
    N := 1;
    g1 := 2*q^(-13) - 2*q^(-2) - 4*q^(-1) + 2*q - 2*q^2 - 2*q^3 + O(q^4);
    g2 := q^(-11) + 2*q^(-7)  -2*q^(-2) + 4*q + 4*q^4 + O(q^5);
    g3 := 2*q^(-26) + 6*q^(-7) - 6*q^(-2) + 2*q^(-1) + 10*q - 8*q^2 + O(q^3);
end function;

// !! TODO - cache the Kappa0's or do it for a bunch of fs simultaneously
// We use a combination of the two versions of Schofer's formula from [GY] and [Err]
// We write sum log|psi|^2 = -|CM(d)|/4 * sum c_m kappa(-m)
// Note that in [GY] there is no square on the lhs, and 
// in [Err] there is no division by 4 on the rhs,
// but this seems to match with the examples in [Err] !?
intrinsic SchoferFormula(fs::SeqEnum[RngSerPuisElt], d::RngIntElt, D::RngIntElt, N::RngIntElt) -> SeqEnum[Assoc]
{Return the log of the absolute value of f for every f in fs at the CM point with CM d, as a sequence of associative array}
    L, Ldual, disc_grp, to_disc, Qinv := ShimuraCurveLattice(D,N);
    Q := ChangeRing(Qinv^(-1), Integers());
    OK := MaximalOrder(QuadraticField(d));
    is_sqr, cond := IsSquare(d div Discriminant(OK));
    assert is_sqr;
    require cond eq 1 : "Not implemented for non-maximal orders!";
    O := sub<OK | cond>;
    n_d := NumberOfOptimalEmbeddings(O, D, N);
    W_size := 2^#PrimeDivisors(D*N);
    // Not sure?? Think this what happens to the number of CM points on the full quotient
    sqfree, sq := SquarefreeFactorization(d);
    Ogg_condition := (cond eq 1) or (((sqfree mod 4 eq 1) and (cond eq 2)));
    if ((D*N) mod sqfree eq 0) and Ogg_condition then
        W_size div:= 2;
    end if;
    // n := -Valuation(f);
    ns := [-Valuation(f) : f in fs];
    n := Maximum(ns);
    log_coeffs := [AssociativeArray() : f in fs];
    for m in [1..n] do
        if &and[Coefficient(f, -m) eq 0 : f in fs] then continue; end if;
        log_coeffs_m := Kappa0(m,d,Q);
        for p in Keys(log_coeffs_m) do
            if (log_coeffs_m[p] ne 0) then
                for i->f in fs do
                    if not IsDefined(log_coeffs[i], p) then
                        log_coeffs[i][p] := 0;
                    end if;
                    log_coeffs[i][p] +:= Coefficient(f,-m)*log_coeffs_m[p];
                end for;
            end if;
        end for;
    end for;
    for i in [1..#fs] do
        for p in Keys(log_coeffs[i]) do
            log_coeffs[i][p] *:= -n_d / (4*W_size);
        end for;
    end for;
    return log_coeffs;
end intrinsic;

intrinsic SchoferFormula(f::RngSerPuisElt, d::RngIntElt, D::RngIntElt, N::RngIntElt) -> Assoc
{Return the log of the absolute value of f at the CM point with CM d, as an associative array}
    return SchoferFormula([f], d, D, N)[1];
end intrinsic;

procedure test_Kappa0()
    D := 6;
    N := 1;
    L, Ldual, disc_grp, to_disc, Qinv := ShimuraCurveLattice(D,N);
    Q := ChangeRing(Qinv^(-1), Integers());
    // verifying [Yang, Example 21, p. 24-25]
    // assert Round(1/Kappa0(3,-4,Q)) eq 2^8*3^4;
    log_coeffs := Kappa0(3,-4,Q);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, -8>, <3, -4> };
    // assert Round(1/Kappa0(1,-3,Q)) eq 2^4;
    log_coeffs := Kappa0(1,-3,Q);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, -4> };
    // verifying [Err, p. 850]
    // assert Round(1/Kappa0(1,-24,Q)) eq 2^6;
    log_coeffs := Kappa0(1,-24,Q);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, -6> };
    // assert Round(1/Kappa0(3,-24,Q)) eq 2^8*3^4;
    log_coeffs := Kappa0(3,-24,Q);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, -8>, <3, -4> };
    // assert Round(1/Kappa0(1,-163,Q)) eq 2^4*3^11*7^4*19^4*23^4;
    log_coeffs := Kappa0(1,-163,Q);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, -4>, <3, -11>, <7, -4>, <11,0>, <19,-4>, <23,-4> };
    // assert Round(1/Kappa0(3,-163,Q)^3) eq 2^40*3^12*5^12*11^12*17^12;
    log_coeffs := Kappa0(3,-163,Q);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, -40/3>, <3, -4>, <5, -4>, <11,-4>, <17,-4>, <23,0>, <89, 0> };
    D := 10;
    L, Ldual, disc_grp, to_disc, Qinv := ShimuraCurveLattice(D,N);
    Q := ChangeRing(Qinv^(-1), Integers());
    log_coeffs := Kappa0(3,-68,Q);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -8>,  <5, -14/3>};
    // This does not work. This is probably a typo! Isn't it supposed to be 2. This would be the relevant number.
    //log_coeffs := Kappa0(1,-68,Q);
    log_coeffs := Kappa0(2,-68,Q);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -6>,  <5, -6>};
    return;
end procedure;

procedure test_Schofer_6()
    _<q> := PuiseuxSeriesRing(Rationals());
    // testing [Errthum, p. 850]
    f6 := -6*q^(-3) + 4*q^(-1) + O(q);
    log_coeffs := SchoferFormula(f6, -24, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <3, -6>, <2, -6> };

    // |t6| = 6^6 | psi_f6 |^2
    log_coeffs := SchoferFormula(f6, -163, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <3, 5>, <2, -16>, <7,4>, <5,-6>, <19,4>, <11,-6>, <23,4>, <17,-6> };

    // testing [Errthum, Table 2] (Section 8.1.1)
    log_coeffs := SchoferFormula(f6, -40, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <3, 1>, <2, -6>, <5,-3> };

    log_coeffs := SchoferFormula(f6, -19, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <3, 1>, <2, -16> };

    log_coeffs := SchoferFormula(f6, -52, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -4>, <3, 1>, <5,-6> };

    log_coeffs := SchoferFormula(f6, -84, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -4>, <3, -9>, <7,2> };

    log_coeffs := SchoferFormula(f6, -88, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -6>, <3, 1>, <5,-6>, <7,4>, <11,-3> };

    // This one does not work! - the problem with zeros
    // log_coeffs := SchoferFormula(f6, -100, 6, 1);
    // assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -2>, <3, 1>, <5,1>, <7,4>, <11,-6> };

    log_coeffs := SchoferFormula(f6, -120, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -6>, <3, -9>, <5,-3>, <7,4> };

    log_coeffs := SchoferFormula(f6, -132, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -2>, <3, -6>, <5,-6>, <11,2> };

    log_coeffs := SchoferFormula(f6, -148, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -4>, <3, 1>, <5,-6>, <7,4>, <11,4>, <17,-6> };

    log_coeffs := SchoferFormula(f6, -168, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -6>, <3, -6>, <5,-6>, <7,2>, <11,4> };

    log_coeffs := SchoferFormula(f6, -43, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -16>, <3, 1>, <5,-6>, <7,4>};

    log_coeffs := SchoferFormula(f6, -51, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -16>, <3, -6>, <7,4>};

    log_coeffs := SchoferFormula(f6, -228, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <3, -12>, <5,-6>, <7,4>, <19,2> };

    log_coeffs := SchoferFormula(f6, -232, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,-6>, <3, 1>, <5,-6>, <7,4>, <11,4>, <19,4>, <23,-6>, <29,-3> };

    log_coeffs := SchoferFormula(f6, -67, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,-22>, <3, 1>, <5,-6>, <7,4>, <11,4> };

    // This does not work - hits mu = m = 0 !!!!???
    // log_coeffs := SchoferFormula(f6, -75, 6, 1);
    // assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,-16>, <3, -9>, <5,-1>, <11,4> };

    log_coeffs := SchoferFormula(f6, -312, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,-6>, <3, -6>, <5,-6>, <7,4>, <11,-6>, <23,4> };

    log_coeffs := SchoferFormula(f6, -372, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,-4>, <3, -9>, <5,-6>, <7,4>, <11,-6>, <19,4>, <31,2> };

    log_coeffs := SchoferFormula(f6, -408, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,-6>, <3, -12>, <5,-6>, <7,4>, <11,4>, <17,-3>, <31,4> };

    log_coeffs := SchoferFormula(f6, -123, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,-16>, <3, -6>, <5,-6>, <7,4>, <19,4> };

    // This one does not work! - the problem with zeros
    // Also see p.828 of [Err]
    // log_coeffs := SchoferFormula(f6, -147, 6, 1);
    // assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,-16>, <3, -9>, <5,-6>, <7,-1>, <11,4>, <23,4> };
    
    log_coeffs := SchoferFormula(f6, -163, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,-16>, <3, 5>, <5,-6>, <7,4>, <11,-6>, <17,-6>, <19,4>, <23,4> };

    log_coeffs := SchoferFormula(f6, -708, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,2>, <3, -6>, <5,-6>, <7,4>, <11,4>, <17,-6>, <29,-6>, <47,4>, <59,2> };

    log_coeffs := SchoferFormula(f6, -267, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,-22>, <3, -6>, <5,-6>, <7,4>, <11,-6>, <31,4>, <43,4> };

    //this is very close but off at 2 and 3?
    // log_coeffs := SchoferFormula(f6, -996, 6, 1);
    // assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, 10>, <3, -6>, <41, -6>, <29, -6>, <17,-6>, <7,12>, <83,2>, <71,4> };

    return;




end procedure;


procedure test_Schofer_10()

    _<q> := PuiseuxSeriesRing(Rationals());

    //This works!
    f10 := 3*q^(-3) - 2*q^(-2) + O(q);
    log_coeffs := SchoferFormula(f10, -20, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, 3> };

    //This does still not work at 2! Off by a factor of 4, this is confusing. The kappas are all correct now.
    log_coeffs := SchoferFormula(f10, -68, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, 2>, <5, 1> };

    //t_10 = 2^(-2)*|Psi_f_10|^2

    //Testing Err. Table 4, Section 8.2.1
    //this works!
    log_coeffs := SchoferFormula(f10, -40, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <3,3>, <2, 2> };

    //this works!
    log_coeffs := SchoferFormula(f10, -52, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <3,3>, <2, 3>, <5, -2> };

    //this is wrong - used to have the 0 error
    log_coeffs := SchoferFormula(f10, -72, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <3,1>, <2, 2>, <5, -2>, <7,-2> };

    //this works!
    log_coeffs := SchoferFormula(f10, -120, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <3,3>, <2, 2>,  <7,-2> };

    //this works!
    log_coeffs := SchoferFormula(f10, -88, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <3,3>, <2, 1>, <5,3>, <7,-2> };

    //this is wrong -used to have the 0 error
    log_coeffs := SchoferFormula(f10, -27, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <3,1>, <2, 8>, <5,-2> };

    //this works!
    log_coeffs := SchoferFormula(f10, -35, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {  <2, 8>, <7,-1> };
    //this works!
    log_coeffs := SchoferFormula(f10, -148, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {  <2, 3>, <3,3>, <11,3>, <5,-2>, <7,-2>, <13,-2> };
    //this works!
    log_coeffs := SchoferFormula(f10, -43, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {  <2, 8>, <3,3>,  <5,-2>, <7,-2>};

    //This (-180) doesn't work! Why??? This was proved in Elkies, also
    log_coeffs := SchoferFormula(f10, -180, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {  <2, 3>,  <11,3>,<13,-2>};

    // this works!
    log_coeffs := SchoferFormula(f10, -232, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq {   <3,3>,  <11,3>,<17,3>, <5,-2>, <7,-2>, <23,-2>};

    //this works!
   log_coeffs := SchoferFormula(f10, -67, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {   <2,8>,  <3,3>,  <5,3>, <7,-2>, <13,-2>};

    //this works!
   log_coeffs := SchoferFormula(f10, -280, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {   <2,1>,  <3,3>, <11,3>, <7,-1>, <23,-2>};
    //this works!
   log_coeffs := SchoferFormula(f10, -340, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {   <2,3>,  <3,3>, <23,3>, <7,-2>, <29,-2>};

    //this works!
   log_coeffs := SchoferFormula(f10, -115, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {   <2,11>,  <3,3>, <13,-2>, <23,-1>};
    //this works!
   log_coeffs := SchoferFormula(f10, -520, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {   <2,-1>,  <3,3>, <29,3>, <7,-2>, <13,-1>, <47,-2>};

    //this works but takes a long time!
   log_coeffs := SchoferFormula(f10, -163, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {   <2,11>,  <3,3>, <5,3>, <11,3>, <7,-2>, <13,-2>,<29,-2>, <31,-2> };

    //this works!
   log_coeffs := SchoferFormula(f10, -760, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {   <2,2>,  <3,3>, <17,3>, <47,3>, <7,-2>, <31,-2>,<71,-2> };
    //this works!
   log_coeffs := SchoferFormula(f10, -235, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {   <2,8>,  <3,3>, <17,3>,  <7,-2>, <37,-2>,<47,-1> };


end procedure;


procedure test_Schofer_26()

    //Test as in Guo Yang p.25
    _<q> := PuiseuxSeriesRing(Rationals());

    g1 := 2*q^(-13) - 2*q^(-2) - 4*q^(-1) + O(q);
    g2 := q^(-11) + 2*q^(-7) - 2*q^(-2) + O(q);
    g3 := 2*q^(-26) + 6*q^(-7)- 6*q^(-2) +2*q^(-1) + O(q);

    log_coeffs := SchoferFormula(g1, -11, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq {  };

    log_coeffs := SchoferFormula(g1, -19, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <3,2>  };

    log_coeffs := SchoferFormula(g1, -20, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <5,1>  };

    log_coeffs := SchoferFormula(g1, -24, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <3,1>  };

    log_coeffs := SchoferFormula(g1, -67, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <3,4>, <5,-2>  };


    log_coeffs := SchoferFormula(g2, -19, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,6>  };

    log_coeffs := SchoferFormula(g2, -20, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,5>  };

    log_coeffs := SchoferFormula(g2, -24, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,5>  };

    log_coeffs := SchoferFormula(g2, -52, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,3>  };

    log_coeffs := SchoferFormula(g2, -67, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,6>, <5,-2>, <7,1>  };


    log_coeffs := SchoferFormula(g3, -11, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,10>,<11,1>,<13,3> };

    log_coeffs := SchoferFormula(g3, -19, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,10>,<13,3>,<19,1>  };

    log_coeffs := SchoferFormula(g3, -20, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,12>,<13,3>  };

    log_coeffs := SchoferFormula(g3, -24, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,13>,<13,3>  };

    log_coeffs := SchoferFormula(g3, -52, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,6>,<13,5>  };

    log_coeffs := SchoferFormula(g3, -67, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,10>, <5,-6>, <13,3>,<41,2>,<67,1>  };

end procedure;