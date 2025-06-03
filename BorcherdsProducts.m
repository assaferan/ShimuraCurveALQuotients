// Compute the order of vanishing of eta(delta tau) at the cusp a/b in Gamma_0(N)
// multiplied by 24
function order_of_vanishing_of_eta(delta, b, N)
    return N*GCD(b,delta)^2 div (GCD(N,b^2)*delta);
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

function integer_programming_input(D,N : n := -1)
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
    if n eq -1 then
        n := n0 + k;
    end if;
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

function write_polymake_scriptfile(D,N : n := -1)
    eqs, ieqs, t := integer_programming_input(D,N : n := n);
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

function get_integer_prog_solutions(D,N : n := -1)
    t := write_polymake_scriptfile(D,N : n := n);
    fname := Sprintf("polymake_script_%o_%o", D, N);
    polymake := Read(POpen("polymake --script " cat fname, "r"));
    sol_lines := Split(polymake, "\n");
    sol_vecs := [Split(line, " ") : line in sol_lines];
    sols := [[eval(x) : x in vec] : vec in sol_vecs];
    M := 2*D*N;
    rs := [sol[2..1 + #Divisors(M)] : sol in sols];
    return rs, Eltseq(t)[1..#Divisors(M)];
end function;

procedure update_precision_eta(~nor_eta_ds, Prec, M)
    _<q> := LaurentSeriesRing(Rationals());
    nor_eta := &*[1 - q^n : n in [1..Prec-1]] + O(q^(Prec));
    // eta_ds := [Evaluate(nor_eta, q^d)*q^(d/24) : d in Divisors(M)];
    nor_eta_ds := [Evaluate(nor_eta, q^d) + O(q^Prec) : d in Divisors(M)];
    return;
end procedure;

intrinsic WeaklyHolomorphicBasis(D::RngIntElt,N::RngIntElt : Prec := 100) -> .
{returns a weakly holomorphic basis corresponding to D, N.}
    D0,M,g := get_D0_M_g(D,N);
    update_precision_eta(~nor_eta_ds, Prec, M);
    R<q> := Universe(nor_eta_ds);
    rk := -1;
    dim := 0;
    n := -1;
    while (rk lt dim) do
        print "prec = ", Prec;
        print "n = ", n;
        print "rk = ", rk;
        print "dim = ", dim;
        rs, t := get_integer_prog_solutions(D,N : n := n);
        eta_quotients := [&*[nor_eta_ds[i]^r[i] : i->d in Divisors(M)] * q^(&+[d*r[i] : i->d in Divisors(M)] div 24) : r in rs];
        t_eta_quotient := &*[nor_eta_ds[i]^t[i] : i->d in Divisors(M)] * q^(&+[d*t[i] : i->d in Divisors(M)] div 24);
        k := -Valuation(t_eta_quotient);
        min_v := Minimum([Valuation(eta_quot) : eta_quot in eta_quotients]);
        n := -min_v;
        min_prec := Minimum([RelativePrecision(eta_quot) - min_v + Valuation(eta_quot) : eta_quot in eta_quotients]);
        coeffs := Matrix([AbsEltseq(q^(-min_v)*(R!eta_quo) : FixedLength)[1..min_prec] : eta_quo in eta_quotients]);
        E, T := EchelonForm(coeffs);
        dim := n + &+[d div 4 : d in Divisors(D0)] + 1 - g;
        rk := Rank(E);
        if (dim gt Prec) then
            Prec := dim;
            update_precision_eta(~nor_eta_ds, Prec, M);
        else
            Prec +:= k;
            update_precision_eta(~nor_eta_ds, Prec, M);
            // n *:= 2;
            // update n
            n +:= k;
        end if;
    end while;
    // n div:= 2;
    n -:= k;
    // sanity checks
    assert rk eq dim;
    n_gaps := g - &+[d div 4 : d in Divisors(D0)];
    pole_orders := [PivotColumn(E,i) - n - 1 : i in [1..Rank(E)]];
    assert (n + 1 - #pole_orders) eq n_gaps;
    return E, n, t_eta_quotient;
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
end function;*/

intrinsic FindLambda(Q::AlgMatElt, d ::RngIntElt: bound := 10)->.
    {}
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
    eps_mu := func<k | LegendreSymbol(-1,p)^(l_mu(k) div 2) * &*[Integers() | my_legendre_symbol(eps[i],p) : i in L_mu(k)]>;
    f_1 := function(x)
        a, alpha := Valuation(x,p);
        return IsEven(l_mu(a+1)) select -1/p else my_legendre_symbol(alpha,p) / sqrtp;
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
    km_sqr := -d*scale_sqr*(w*W_kron / h)^2;
    km_sign := -Sign(W_kron);

    return km_sqr, km_sign;
   
end function;

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

    vprintf ShimuraQuotients, 2 : "adding %o log %o\n", -ret, p_prime;
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
intrinsic Kappa0(m::RngIntElt, d::RngIntElt, Q::AlgMatElt, lambda_v::ModTupRngElt) -> FldReElt
{Computing coefficients Kappa0(m) in Schofers formula}
    vprintf ShimuraQuotients, 1:"Kappa0 of %o\n", m;
    Q := ChangeRing(Q, Integers());
    c_Lplus := Content(lambda_v);
    Lplus := RSpaceWithBasis(Matrix(lambda_v div c_Lplus));
    Lminus := Kernel(Transpose(Matrix(lambda_v*Q)));
    L := RSpaceWithBasis(IdentityMatrix(Integers(),3));
    L_quo, L_quo_map := L / (Lplus + Lminus);
    log_coeffs := AssociativeArray();
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

// !! TODO - cache the Kappa0's or do it for a bunch of fs simultaneously
// We use a combination of the two versions of Schofer's formula from [GY] and [Err]
// We write sum log|psi|^2 = -|CM(d)|/4 * sum c_m kappa(-m)
// Note that in [GY] there is no square on the lhs, and 
// in [Err] there is no division by 4 on the rhs,
// but this seems to match with the examples in [Err] !?
intrinsic SchoferFormula(fs::SeqEnum[RngSerLaurElt], d::RngIntElt, D::RngIntElt, N::RngIntElt) -> SeqEnum[Assoc]
{Return the log of the absolute value of f for every f in fs at the CM point with CM d, as a sequence of associative array}
    L, Ldual, disc_grp, to_disc, Qinv := ShimuraCurveLattice(D,N);
    Q := ChangeRing(Qinv^(-1), Integers());
    OK := MaximalOrder(QuadraticField(d));
    is_sqr, cond := IsSquare(d div Discriminant(OK));
    assert is_sqr;
    require cond eq 1 : "Not implemented for non-maximal orders!";
    O := sub<OK | cond>;
    n_d := NumberOfOptimalEmbeddings(O, D, N);
    require n_d gt 0 : "Curve does not have a CM point of discirminant d!";
    W_size := 2^#PrimeDivisors(D*N);
    // Not sure?? Think this what happens to the number of CM points on the full quotient
    sqfree, sq := SquarefreeFactorization(d);
    Ogg_condition := (cond eq 1) or (((sqfree mod 4 eq 1) and (cond eq 2)));
    if ((D*N) mod sqfree eq 0) and Ogg_condition then
        W_size div:= 2;
    end if;
    // n := -Valuation(f);
    bd := 10;
    found_lambda := false;
    while not found_lambda do
        bd *:= 2;
        found_lambda, lambda_v := FindLambda(Q,-d : bound := bd);
    end while;
    assert found_lambda;

    ns := [-Valuation(f) : f in fs];
    n := Maximum(ns);
    log_coeffs := [AssociativeArray() : f in fs];
    for m in [1..n] do
        if &and[Coefficient(f, -m) eq 0 : f in fs] then continue; end if;
        log_coeffs_m := Kappa0(m,d,Q,lambda_v);
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

intrinsic SchoferFormula(f::RngSerLaurElt, d::RngIntElt, D::RngIntElt, N::RngIntElt) -> Assoc
{Return the log of the absolute value of f at the CM point with CM d, as an associative array}
    return SchoferFormula([f], d, D, N)[1];
end intrinsic;

intrinsic AbsoluteValuesAtRationalCMPoint(fs::SeqEnum[RngSerLaurElt], d::RngIntElt, Xstar::ShimuraQuot) -> SeqEnum[ExtReElt]
{Returns the absolute value of f for every f in fs at the rational CM point with CM d.}
    vals := [ExtendedReals() | 1 : f in fs];
    for i->f in fs do
        div_f := DivisorOfBorcherdsForm(f, Xstar);
        in_support := exists(pt){pt : pt in div_f | pt[1] eq d};
        if in_support then
            if pt[2] lt 0 then vals[i] := Infinity(); end if;
            if pt[2] gt 0 then vals[i] := 0; end if;
        end if;
    end for;
    rest_idxs := [i : i in [1..#fs] | vals[i] eq 1];
    if IsEmpty(rest_idxs) then return vals; end if;
    rest_fs := [fs[i] : i in rest_idxs];
    log_coeffs := SchoferFormula(rest_fs, d, Xstar`D, Xstar`N);
    for i->log_coeff in log_coeffs do
        require &and[IsIntegral(log_coeff[p]) :  p in Keys(log_coeff)] : "Increase the precision on the Borcherds forms";
        vals[rest_idxs[i]] := &*[Rationals() | p^(Integers()!log_coeff[p]) : p in Keys(log_coeff)];
    end for;
    return vals;
end intrinsic;

intrinsic BorcherdsForms(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot] : Prec := 100) -> Assoc
{Returns weakly holomorphic modular forms with divisors that are the ramification divisors of each of the double covers in curves,
along with two different hauptmoduls.}
    rams := RamficationPointsOfCovers(Xstar, curves);
    D0,M,g := get_D0_M_g(Xstar`D,Xstar`N);
    n0 := Maximum(2*g-2-&+[d div 4 : d in Divisors(D0)],0);
    E, n, t := WeaklyHolomorphicBasis(Xstar`D, Xstar`N : Prec := Prec);
    k := -Valuation(t);
    E := Submatrix(E, [1..Rank(E)], [1..Ncols(E)]);
    // _<q> := Parent(t);
    R<q> := LaurentSeriesRing(Integers());
    t := R!t;
    fs_E := [q^(-n)*&+[(Integers()!b[i])*q^(i-1) : i in [1..Ncols(E)]] : b in Rows(E)];
    pts := RationalCMPoints(Xstar); // pts <-> infty, 0, rational
    //we do this twice -- we should remember this
    found := false;
    infty_idx := 1;
    while (not found) do
        // infty := pts[1];
        infty := pts[infty_idx];
        non_infty := [pt : pt in pts | pt ne infty];
        for other_pts in CartesianPower(non_infty,2) do
            if other_pts[1] eq other_pts[2] then continue; end if;
            rams[-1] := [other_pts[1]];
            rams[-2] := [other_pts[2]];
            fs := AssociativeArray();
            found := true;
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
                ms := [(d[1] mod 4 eq 0) select d[1] div 4 else d[1] : d in ram];
                min_m := Minimum(ms);
                min_m := Minimum(min_m, -(n0 + k - 1));
                V := RSpace(Integers(),-min_m+1);
                v := &+[div_coeffs[i]*ram[i][2]*V.(m-min_m+1) : i->m in ms];
                r := (-min_m - n0) div k;
                s := (-min_m) - r*k;
                basis_n0 := fs_E[n+2-n0..#fs_E]; // basis for M_{n0-1}^!
                init_basis := fs_E[n+2-n0-k..n+1-n0]; // completing to a basis for M_{n_0+k-1}^!
                // full_basis is a basis for M_{-min_m}^!(4D_0)
                full_basis := [t^r*f + O(q) : f in init_basis[n0+k-s..#init_basis]];
                full_basis cat:= &cat[[t^(r-1-j)*f + O(q) : f in init_basis] : j in [0..r-1]];
                full_basis cat:= [f + O(q) : f in basis_n0];
                coeffs := Matrix([AbsEltseq(q^(-min_m)*f : FixedLength) : f in full_basis]);
                ech_basis := EchelonForm(coeffs);
                ech_fs := [q^min_m*&+[(Integers()!b[i])*q^(i-1) : i in [1..Ncols(ech_basis)]]+O(q) : b in Rows(ech_basis)];
                relevant_ds := [0];
                for d in [1..-min_m] do
                    discs := [4*d div r2: r2 in Divisors(4*d) | IsSquare(r2)];
                    discs := [disc : disc in discs | disc mod 4 in [0,3]];
                    n_d := 0;
                    for disc in discs do
                        S := QuadraticOrder(BinaryQuadraticForms(-disc));
                        n_d +:= NumberOfOptimalEmbeddings(S,Xstar`D,Xstar`N);
                    end for;
                    if (n_d ne 0) then
                        Append(~relevant_ds, d);
                    end if;
                end for;
                cols := Reverse([-d-min_m+1 : d in relevant_ds]);
                target_v := Vector([v[c] : c in cols]); 
                // f := &+[div_coeffs[i]*ram[i][2]*ech_basis[m-min_m+1] : i->m in ms | -m ge n0];
                // coeffs_trunc := Submatrix(coeffs, [1..Nrows(coeffs)],cols);
                coeffs_trunc := Submatrix(ech_basis, [1..Nrows(coeffs)],cols);
                if target_v notin Image(coeffs_trunc) then
                    found := false;
                    break;
                end if;
                sol := Solution(coeffs_trunc, target_v);
                fs[i] := &+[sol[i]*ech_fs[i] : i in [1..#ech_fs]]; // !! Make sure that f has the correct divisor - right a function for Lemma 22
            end for;
            if found then break; end if;
        end for;
        infty_idx +:= 1;
    end while;  
    return fs;
end intrinsic;

intrinsic DivisorOfBorcherdsForm(f::RngSerLaurElt, Xstar::ShimuraQuot) -> SeqEnum
{Return the divisor of the Borcherds form associated to f.}
    N := Valuation(f);
    ells := NumberOfEllipticPointsByCMOrder(Xstar);
    ell_lut := AssociativeArray();
    for q in Keys(ells) do
        for d in Keys(ells[q]) do
            ell_lut[d] := q;
        end for;
    end for;
    divisor := [];
    for m in [N..-1] do
        c_m := Coefficient(f,m);
        if c_m eq 0 then continue; end if;
        sq_divs := [r2 : r2 in Divisors(-4*m) | IsSquare(r2)];
        for r2 in sq_divs do
            d := (4*m) div r2;
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



intrinsic CandidateDiscriminants(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot] : Exclude := {}, bd := 4) -> SeqEnum
    {Returns list of candidate discriminats for Schofer's formula}
    cm_pts := RationalCMPoints(Xstar : Exclude := Exclude);
    cm_pts := Reverse(Sort(cm_pts));

    quad_cm := QuadraticCMPoints(Xstar : Exclude := Exclude, bd := bd);
    quad_cm := Reverse(Sort(quad_cm));

    return [cm_pts, quad_cm];

end intrinsic;

intrinsic AbsoluteValuesAtCMPoints(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot], all_cm_pts ::SeqEnum, fs::Assoc : MaxNum := 7, Prec := 100,  Exclude := {}, Include := {}) -> SchoferTable
{Returns the absolute values of y^2 for all degree 2 covers and two hauptmodules at CM points.}
    cm_pts := all_cm_pts[1];
    quad_cm := all_cm_pts[2];

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
    if #cm_pts gt need then  
    //need to make space for include points, but otherwise fill up with rational points as much as possible
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
    for pt in pt_list_rat do
        d := pt[1];
        vals := AbsoluteValuesAtRationalCMPoint(all_fs, d, Xstar);
        for i->v in vals do
            Append(~table[i], vals[i]);
        end for;
    end for;

    for pt in pt_list_quad do
        d := pt[1];
        norm_val := AbsoluteValuesAtRationalCMPoint(all_fs, d, Xstar);
        for i->v in norm_val do
            Append(~table[i], norm_val[i]);
        end for;
    end for;

    ds := [[pt[1] : pt in pt_list_rat], [ pt[1] : pt in pt_list_quad ]];

    schofer_tab := CreateSchoferTable(table, keys_fs, ds, curves);
    schofer_tab`Xstar := Xstar;
    schofer_tab`BorcherdsForms := fs;

    return schofer_tab;
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
    return s_new, stilde_new;
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

    all_relns := [];
    for j->idx in k_idxs do
        B := kernels[j];
        numcoeffs := #Eltseq(B[1]);
        R<[x]>:=PolynomialRing(Rationals(),#B);
        inf_idx_y2 := Index(table[idx], Infinity());
        inf_idx_s := Index(table[s_idx], Infinity());
        y2vals := Remove(table[idx], inf_idx_y2);
        svals := Remove(table[s_idx],inf_idx_s);
        quad_svals := [s :  i->s in svals   | Type(s) eq RngUPolElt];
        quad_sidxs := [i :  i->s in svals   | Type(s) eq RngUPolElt];
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
    curves := schofer_tab`curves;
    all_fs := [fs[k] : k in keys_fs];
    norm_val := AbsoluteValuesAtRationalCMPoint(all_fs, dnew, Xstar);
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
    UpdateFieldsOfDefn(schofer_tab, dnew);
    return;
end procedure;

function solve_quadratic_constraints(relns)
    R := Universe(relns);
    P := ProjectiveSpace(R);
    S := Scheme(P, relns);
    assert Dimension(S) eq 0;
    return RationalPoints(S);
end function;


intrinsic EquationsOfCovers(schofer_table::SchoferTable, curves::SeqEnum[ShimuraQuot], all_cm_pts::SeqEnum) -> SeqEnum, SeqEnum
{Determine the equations of the covers using the values from Schofers formula}
    R<x> := PolynomialRing(Rationals());
    eqn_list := [ ];
    keys_fs := schofer_table`Keys_fs;
    table := schofer_table`Values;
    ds := schofer_table`Discs;
    ds := ds[1] cat ds[2];
    k_idxs := schofer_table`K_idxs;

    kernels := RationalConstraintsOnEquations(schofer_table, curves);
    relns := QuadraticConstraintsOnEquations(schofer_table, curves, kernels);

    for i->B in kernels do //indexed by k_idxs
        if #relns[i] eq 0 then
            require #B eq 1 : "Try adding quadratic points -- not enough constraints from the rational points";
            v :=  Eltseq(B[1]);
            monic_v := [-v[i]/v[#v] : i in [1..#v-1]];
            f := R!monic_v;
            Append(~eqn_list, f);
        else
            coeffs := solve_quadratic_constraints(relns[i]);
            require #coeffs eq 1 : "We do not have enough constraints coming from quadratic and rational points";
            // while #coeffs gt 2 do
            //     vprint ShimuraQuotients, 1: "Attempting to add more CM points to narrow down quadratic constraints";
            //     candidates := Set([pt[1] : pt in all_cm_pts]) diff Set(ds);
            //     d := Representative(candidates);
            //     add_new_column(schofer_table,d);
            //     new_ds := schofer_table`Discs;
            //     ds := new_ds[1] cat new_ds[2];
            //     table := schofer_table`Values;
            //     flds := schofer_table`FldsOfDefn;
            //     F := flds[keys_fs[k_idxs[k]]][d];
            //     coeffs := solve_quadratic_constraints(relns[i]);
            // end while;
            coeffs := Eltseq(coeffs[1]);
            v := &+[B[j]*coeffs[j] : j in [1..#B]];
            v := Eltseq(v);
            monic_v := [-v[i]/v[#v] : i in [1..#v-1]];
            f := R!monic_v;
            Append(~eqn_list, f);
        end if;
    end for;

    return eqn_list, [k : k in keys_fs | k gt 0];
end intrinsic;

intrinsic EquationsOfCovers(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot] : Prec := 100) -> SeqEnum, SeqEnum
{Determine the equations of the immediate covers of X.}
    fs := BorcherdsForms(Xstar, curves : Prec := Prec);
    d_divs := &cat[[T[1]: T in  DivisorOfBorcherdsForm(f, Xstar)] : f in [fs[-1], fs[-2]]]; //include zero infinity of hauptmoduls
    all_cm_pts := CandidateDiscriminants(Xstar, curves);
    genus_list := [curves[i]`g: i in Xstar`CoveredBy];
    num_vals := Maximum([2*g+5 : g in genus_list]);
    abs_schofer_tab := AbsoluteValuesAtCMPoints(Xstar, curves, all_cm_pts, fs : MaxNum := num_vals, Prec := Prec, Exclude := {}, Include := Set(d_divs));
    ReduceTable(abs_schofer_tab);
    schofer_tab := ValuesAtCMPoints(abs_schofer_tab, curves, all_cm_pts);
    return EquationsOfCovers(schofer_tab, curves, all_cm_pts);
end intrinsic;

intrinsic EquationsAboveP1s(eqn_list::SeqEnum[RngUPolElt[FldRat]], new_keys::SeqEnum[RngIntElt], curves::SeqEnum[ShimuraQuot]) -> SeqEnum, SeqEnum
    {Using Riemann Roch, leverage covered equations to get higher cover equations}
    P1s := [<i, keys> : i->keys in new_keys | Degree(eqn_list[i]) eq 1];
    conics := [<i, keys> : i->keys in new_keys | curves[keys]`g eq 0 and Degree(eqn_list[i]) ne 1];

    curves_above_P1s := AssociativeArray();
    curves_above_conics := AssociativeArray();

    for pair in P1s do
        for c in curves[pair[2]]`CoveredBy do
            curves_above_P1s[c] := pair[1];
        end for;
    end for;

    for pair in conics do
        for c in curves[pair[2]]`CoveredBy do
            if not IsDefined(curves_above_P1s, c) then
                curves_above_conics[c] := pair[1];
            end if;
        end for;
    end for;

    cover_eqns := [];
    cover_keys := [];
    while (not IsEmpty(Keys(curves_above_P1s))) or (not IsEmpty(Keys(curves_above_conics))) do
        for label in Keys(curves_above_P1s) do
            g := curves[label]`g;
            covered_P1 := eqn_list[curves_above_P1s[label]];
            allgplus1covers := { new_keys[i] :  i in [1..#new_keys] | Degree(eqn_list[i]) eq g+1 } meet curves[label]`Covers;
            if #allgplus1covers eq 0 then
                continue;
            end if;
            covered_gplus1_key := Representative(allgplus1covers);
            gplus1idx := Index(new_keys,covered_gplus1_key);
            covered_gplus1 := eqn_list[gplus1idx];
            //if this is empty then it's not hyperelliptic
            c0 := Coefficient(covered_P1, 0);
            c1 := Coefficient(covered_P1, 1);
            _<x> := Parent(covered_gplus1);
            eqn := Evaluate(covered_gplus1, (x^2 - c0)/c1);
            Append(~cover_eqns, eqn);
            Append(~cover_keys, label);
        end for;
        for label in Keys(curves_above_conics) do
            g := curves[label]`g;
            covered_conic := eqn_list[curves_above_conics[label]];
            allgplus1covers := { new_keys[i] :  i in [1..#new_keys] | Degree(eqn_list[i]) eq g+1 } meet curves[label]`Covers;
            if #allgplus1covers eq 0 then
                continue;
            end if;
            covered_gplus1_key := Representative(allgplus1covers);
            gplus1idx := Index(new_keys,covered_gplus1_key);
            covered_gplus1 := eqn_list[gplus1idx];
            //if this is empty then it's not hyperelliptic
            assert Degree(covered_conic) eq 2;
            P2<x,y,z> := ProjectiveSpace(Rationals(),2);
            C := Curve(P2, y^2-z^2*Evaluate(covered_conic, x/z));
            C := Conic(C);
            assert HasRationalPoint(C); // for now not implemented if C does not have a rational point
            C_to_P1 := Inverse(Parametrization(C));
            s_param := C_to_P1(x) / C_to_P1(z);
            _<x> := Parent(covered_gplus1);
            eqn := Evaluate(covered_gplus1, s_param);
            Append(~cover_eqns, eqn);
            Append(~cover_keys, label);
        end for;
        curves_above_P1s := AssociativeArray();
        P1s := [<i, keys> : i->keys in cover_keys | Degree(cover_eqns[i]) eq 1];
        curves_above_P1s := AssociativeArray();
        for pair in P1s do
            for c in curves[pair[2]]`CoveredBy do
                curves_above_P1s[c] := pair[1];
            end for;
        end for;
        curves_above_conics := AssociativeArray();
        conics := [<i, keys> : i->keys in cover_keys | curves[keys]`g eq 0 and Degree(cover_eqns[i]) ne 1];
        curves_above_conics := AssociativeArray();
        for pair in conics do
            for c in curves[pair[2]]`CoveredBy do
                curves_above_conics[c] := pair[1];
            end for;
        end for;
    end while;
    return cover_eqns, cover_keys;

end intrinsic;

intrinsic AllEquationsAboveCovers(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot] : Prec := 100)-> SeqEnum, SeqEnum
    {Get equations of all covers (not just immediate covers)}
    fs := BorcherdsForms(Xstar, curves : Prec := Prec);
    d_divs := &cat[[T[1]: T in  DivisorOfBorcherdsForm(f, Xstar)] : f in [fs[-1], fs[-2]]]; //include zero infinity of hauptmoduls
    all_cm_pts := CandidateDiscriminants(Xstar, curves);
    genus_list := [curves[i]`g: i in Xstar`CoveredBy];
    num_vals := Maximum([2*g+5 : g in genus_list]);
    abs_schofer_tab := AbsoluteValuesAtCMPoints(Xstar, curves, all_cm_pts, fs : MaxNum := num_vals, Prec := Prec, Exclude := {}, Include := Set(d_divs));
    ReduceTable(abs_schofer_tab);
    schofer_tab := ValuesAtCMPoints(abs_schofer_tab, curves, all_cm_pts);
    eqn_list, new_keys := EquationsOfCovers(schofer_tab, curves, all_cm_pts);
    cover_eqns, cover_keys := EquationsAboveP1s(eqn_list, new_keys, curves);
    all_eqns := eqn_list cat cover_eqns;
    all_keys := new_keys cat cover_keys;
    return all_eqns, all_keys;
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
        fraka := mPicR(a@@PicR_to_A);
        B_fraka := QuaternionAlgebra(Rationals(), d, m*Norm(fraka));
        if IsIsomorphic(B_fraka, B) then
            break;
        end if;
    end for;
    assert assigned fraka;
    sigma_a := rec(fraka);
    abs_sig_a := H_R_to_abs^(-1)*sigma_a*H_R_to_abs;
    // _, K_to_abs := IsSubfield(K, abs_H_R);
    _, cc := HasComplexConjugate(abs_H_R);
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
                prod := al_mul(w,prod,D*N);
                al_action[prod] := al_action[prev_prod]*al_action[w];
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

procedure replace_column(schofer_tab, d, dnew)
    //add column associated to dnew remove column associated to d
    table := schofer_tab`Values;
    ds := schofer_tab`Discs;
    Xstar := schofer_tab`Xstar;
    fs := schofer_tab`BorcherdsForms;
    keys_fs := schofer_tab`Keys_fs;
    row_scales := schofer_tab`RowScales;
    all_fs := [fs[k] : k in keys_fs];
    if d in ds[2] then
        deg :=2;
        i := Index(ds[2], d);
        d_idx := #ds[1]+i;
        ds[2][i] := dnew;
    elif d in ds[1] then
        deg := 1;
        i := Index(ds[1], d);
        d_idx := i;
        ds[1][i] := dnew;
    end if;
    norm_val := AbsoluteValuesAtRationalCMPoint(all_fs, dnew, Xstar);
    for i->v in norm_val do
        table[i][d_idx] := norm_val[i]/row_scales[i]^deg;
    end for;
    schofer_tab`Values := table;
    schofer_tab`Discs := [ds[1], ds[2]];
    curves := schofer_tab`curves;
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


    scale_factors :=[];
    for i in k_idxs do
        if exists(j1){j : j->d1 in ratds  | #fldsofdef[keys_fs[i]][d1] eq 1 and Degree(fldsofdef[keys_fs[i]][d1][1]) eq 1 and table[1][j] ne Infinity() and table[1][j] ne 0} then
            //then we have a rational point on X
            d1 := ratds[j1];
            v1 := table[i][j1];
            scale, _ := SquareFree(v1);
            Append(~scale_factors, AbsoluteValue(scale));
        else 
            assert exists(j1){j : j->d1 in ratds  | #fldsofdef[keys_fs[i]][d1] le 2 and {Degree(fldsofdef[keys_fs[i]][d1][k]) : k in [1..#fldsofdef[keys_fs[i]][d1]]} subset {1,2} and table[1][j] ne Infinity() and table[1][j] ne 0};
            assert exists(j2){j : j->d2 in ratds  | #fldsofdef[keys_fs[i]][d2] le 2 and {Degree(fldsofdef[keys_fs[i]][d2][k]) : k in [1..#fldsofdef[keys_fs[i]][d2]]} subset {1,2}  and table[1][j] ne Infinity() and table[1][j] ne 0 and ratds[j1] ne d2};
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
            scale1, _ := SquareFree(v1/d1); //two possibilities
            scale2, _ := SquareFree(v1);
            if IsSquare(scale1*v2/d2) then
                Append(~scale_factors, AbsoluteValue(scale1));
            else
                assert IsSquare(scale2*v2);
                Append(~scale_factors, AbsoluteValue(scale2));
            end if;
        end if;
    end for;
    return scale_factors;

end function;

function find_y2_signs(table, keys_fs, curves, d, j, scale_factors, flds, k_idxs)
    //find signs of y^2 for rational CM point d on each y^2
    //in keys_fs, where j is index of column of d in table
    for k->i in k_idxs do
        if table[i][j] eq Infinity() then continue; end if;
        if table[i][j] eq 0 then continue; end if;

        // fields := FieldsOfDefinitionOfCMPoint(curves[keys_fs[i]], d);
        fields := flds[keys_fs[i]][d];
        Fs_eps := [* <F, eps> : F in fields, eps in [-1,1] *];
        possible_answers := [* *];
        for eps in [-1,1] do
            y2 := scale_factors[k]*eps*table[i][j];
            for F in fields do
                is_sqr, y := IsSquare(F!y2);
                if (is_sqr) then
                    if (Type(F) eq FldRat) or (Degree(F) eq Degree(sub<F|y>)) then
                        Append(~possible_answers, <F,eps,y>);
                    end if;
                end if;
            end for;
        end for;
        // is_sqr := [IsSquare(F!(scale_factors[k]*eps*table[i][j])) : F in fields, eps in [-1,1]];
        // require &or is_sqr : "Value %o at CM point %o does not lie in any subfield.", table[i][j], d;
        // require #[a : a in is_sqr | a] eq 1 : "Too many options for value %o at CM point %o", table[i][j], d;
        assert #possible_answers eq 1;
        // "Wrong number of possible signs or CM fields", possible_answers, table[i][j], d;
        // eps := Fs_eps[Index(is_sqr, true)][2];
        eps := possible_answers[1][2];
        table[i][j] := eps * table[i][j];
    end for;
    return table;
end function;

intrinsic ValuesAtCMPoints(abs_schofer_tab::SchoferTable, curves::SeqEnum[ShimuraQuot], all_cm_pts::SeqEnum) -> SchoferTable
    {}
    ds := abs_schofer_tab`Discs;
    table := abs_schofer_tab`Values;
    keys_fs := abs_schofer_tab`Keys_fs;
    Xstar := abs_schofer_tab`Xstar;
    row_scales := abs_schofer_tab`RowScales;
    flds := abs_schofer_tab`FldsOfDefn;

    ratds := ds[1];
    quadds := ds[2];
    allds := ratds cat quadds;

    s_idx := abs_schofer_tab`sIndex;
    stilde_idx := abs_schofer_tab`sTildeIndex;
    s := table[s_idx];
    stilde := table[stilde_idx];


    table :=[* [* x : i->x in t *] : t in table *];
    s, stilde := find_signs(s, stilde, ds);
    table[s_idx] := s;
    table[stilde_idx] := stilde;
    k_idxs := abs_schofer_tab`K_idxs;

    //Scale the y2 rows of the table
    scale_factors := find_y2_scales(abs_schofer_tab);

    //Next need to go from norms to values on the hauptmoduls for the quad_cm points

    bad_ds := {};
    for i->d in quadds do
        d_idx := #ratds+i;
        good_inds := [];

        while #good_inds ne 1 do
            norm_s := table[s_idx][d_idx];
            norm_stilde := table[stilde_idx][d_idx];
            flds := FieldsOfDefinitionOfCMPoint(Xstar, d);
            assert #flds eq 1;
            assert Degree(flds[1]) eq 2;
            K := flds[1]; //assume fields of definition are exactly quadratic on Xstar
            _<x> := PolynomialRing(Rationals());
            signs := [[1,1], [1,-1],[-1,1],[-1,-1]];
            minpolys := [];
            for eps in signs do
                    trace := 1- eps[1]*norm_stilde +  eps[2]*norm_s;
                    Append(~minpolys, x^2 - trace*x + eps[2]*norm_s);
            end for;
            roots := [Roots(p,K) : p in minpolys];
            good_inds := [i : i->r in roots | #r ne 0 and not(&and[rt[1] in Rationals() : rt in r])];
            //require #good_inds eq 1 :"We need that there is a unique minpoly left after filtering by roots. This is not true at CM point", d;
            Include(~bad_ds, d);
            candidates := Set(all_cm_pts[2]) diff Set(quadds) diff bad_ds;
            assert #candidates ge 1;
            //"No possible choices of CM points left which we can pin down the correct minpoly";
            newd := Representative(candidates);
            replace_column(abs_schofer_tab, d, newd);
            table := abs_schofer_tab`Values;
        end while;
        table[s_idx][d_idx] := minpolys[good_inds[1]];
        norm_s := Coefficient(minpolys[good_inds[1]], 0);
        trace_s := - Coefficient(minpolys[good_inds[1]], 1);
        table[stilde_idx][d_idx] := x^2 - (2- trace_s)*x + (1- trace_s + norm_s);
    end for;
       
    //Find signs on the y2 rows

    //Note that when we are at a quadratic CM point, K the quadratic field the norm is always positive
    //Write v = | y^2(tau)y^2(taubar)|. Then y(tau)y(taubar) = sqrt(eps*norm). Since the fields Q(y(tau)) = Q(y(tau))^sigma
    // where sigma is the unique nontrivial element of Gal(K/Q)
    // y(tau)y(taubar) lies in the fixed field by sigma, i.e. in Q
    // so eps* v is a square in Q, and is positive

    for j->d in ratds do
        table := find_y2_signs(table, keys_fs, curves, d, j, scale_factors, flds, k_idxs);
    end for;

    schofer_table := CreateSchoferTable(table, keys_fs, ds, curves);

    return schofer_table;
end intrinsic;

intrinsic ReduceTable(schofer_tab::SchoferTable)
    {}
    table := schofer_tab`Values;
    sep_ds := schofer_tab`Discs;
    scales := [];
    num_rat_ds := #sep_ds[1];
    for t in table do
        xs := [x : i->x in t | i le num_rat_ds and x notin [0, Infinity()] ];
        //xs := [x : x in t | x notin [0, Infinity()]];
        ps := &join[Set(PrimeDivisors(Numerator(x))) : x in xs];
        ps join:= &join[Set(PrimeDivisors(Denominator(x))) : x in xs];
        ps := [p : p in ps];
        vals := [[Valuation(x,p) : x in xs ] : p in ps];
        mins := [Minimum([<AbsoluteValue(v),v> : v in valp]) : valp in vals];
        scale := &*[Rationals() | p^mins[i][2] : i->p in ps];
        Append(~scales, scale);
    end for;
    degs := [1 : i in sep_ds[1] ] cat [ 2 : i in sep_ds[2]];
    schofer_tab`Values :=  [[x/scales[i]^degs[j] : j->x in t] : i->t in table ];
    schofer_tab`RowScales := scales;
end intrinsic;

intrinsic ValuesAtCMPoints(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot] : MaxNum := 7, Prec := 100) -> SeqEnum, SeqEnum, SeqEnum
{Returns the values of y^2 for all degree 2 covers and two hauptmodules at CM points.}
    fs := BorcherdsForms(Xstar, curves : Prec := Prec);
    d_divs := &cat[[T[1]: T in  DivisorOfBorcherdsForm(f, Xstar)] : f in [fs[-1], fs[-2]]]; //include zero infinity of hauptmoduls
    all_cm_pts := CandidateDiscriminants(Xstar, curves);
    abs_schofer_tab := AbsoluteValuesAtCMPoints(Xstar, curves, all_cm_pts, fs : MaxNum := MaxNum, Prec := Prec, Exclude := {}, Include := Set(d_divs));
    ReduceTable(abs_schofer_tab);
    schofer_tab := ValuesAtCMPoints(abs_schofer_tab, curves, all_cm_pts);
    return schofer_tab;
end intrinsic;

/* Bibliography
[GR] - Gonzales, Rotger - non-elliptic shimura curves of genus one
*/