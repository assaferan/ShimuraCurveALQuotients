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
    fm1 := q^(-1)*Parent(t_eta_quotient)!Eltseq(E[1]);
    f0 := q^(-1)*Parent(t_eta_quotient)!Eltseq(E[2]);
    f0;
    f0*t;
    t;
    t := t_eta_quotient;
    t;
    f0*t;
    f1;
    fm1;
    f0;
    f0*t - 3*f0;
    f0*t - 3*f0 - fm1;
    f0*t^6;
    [f0*t^n : n in [0..6]];
    f0*(t^6 - 8*t^5);
    f0*(t^6 - 8*t^5+17*t^4);
    f0*(t^6 - 8*t^5+17*t^4-4*t^3);
    f0*(t^6 - 8*t^5+17*t^4-4*t^3-8*t^2);
    f0*(t^6 - 8*t^5+17*t^4+4*t^3);
    f0*(t^6 - 8*t^5+17*t^4+4*t^3-37*t^2);
    f0*(t^6 - 8*t^5+17*t^4+4*t^3-37*t^2+22*t);
    fm6 := f0*(t^6 - 8*t^5+17*t^4+4*t^3-37*t^2+22*t);
    f0*(t^3);
    f0*(t^3-5*t^2);
    f0*(t^3-5*t^2+5*t);
    fm3 := f0*(t^3-5*t^2+5*t);
    fm6;
    fm3;
    f0*(t^2);
    f0*(t^2-4*t);
    fm2 := f0*(t^2-4*t);
    fm6;
    fm3;
    fm2;
    f := fm6 - fm3 + fm2 + fm1;
    f;
    psi_s := -2*fm6 + 2*fm1 + 6*f0;
    psi_y := -4*fm6 + 2*fm3 + 2*fm1 + 10*f0;
    psi+s;
    psi_s;
    psi_y;
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
    plus+minus := Matrix([[2,0,1],[2,-1,-1],[-1,0,1]]);
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
    Kernel(Transpose(Matrix([[Trace(x*lambda) : x in basis_L]]));
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
    _<q> := PowerSeriesRing(Rationals());
    eta<q> := DedekindEta(q);  
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

// We woudl like to be able to start with a f,
// and create the coefficients of the Borcherds form F
// obtained from it. 
function SpreadBorcherds(f, L)
    M := Level(L);
    G := PSL2(Integers());
    gens := Generators(G);
    gammas := CosetRepresentatives(Gamma0(M));
    T := gens[1];
    S := gens[2];
    assert Eltseq(T) eq [1,1,0,1];
    assert Eltseq(S) eq [0,1,-1,0];
    s := 0;
    for gamma in gammas do
        word := FindWord(PSL2(Integers()), gamma^(-1));
        assert &*[gens[Abs(word[i])]^Sign(word[i]) : i in [1..#word]] eq gamma^(-1);
        // compute f|gamma rho_L(gamma^(-1)) e_0 and add to the sum
        // Probably needs to remember how f is a linear combination of eta quotients
    end for;
    return s;
end function;