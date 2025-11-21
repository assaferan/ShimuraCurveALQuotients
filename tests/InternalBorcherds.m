import "BorcherdsProducts.m" : Wpoly, Wpoly2, Wpoly_scaled;

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


procedure test_W()
    // testing the few values we know from Yang
    L, Ldual, disc_grp, to_disc, Qinv := ShimuraCurveLattice(6,1);
    Q := ChangeRing(Qinv^(-1), Integers());
    for d in [-3,-4] do
        _, lambda_v := ElementOfNorm(Q,-d);
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

// This is not (!!!) [Err, Lemma 6.1, p. 845] based on [KRY, Lemmas 2.4 and 2.5]
// [Err only refers to primes not in Sm_mu]
// This is based on [KRY, Lemma 2.6]
function Wpolys_self_dual_KRY_2_4(m,p,mu,Lminus,Q)
    BML := BasisMatrix(Lminus);
    Delta := -Determinant(BML*Q*Transpose(BML));
    assert IsZero(mu); // in the self dual case mu is always zero
    K<sqrt_mp> := QuadraticField(-p);
    _<x> := PolynomialRing(K);
    return sqrt_mp * (1 - KroneckerCharacter(p)(Integers()!m) * x^(Valuation(m,p) + 1));
end function;

// This is based on [KY, Proposition 5.2]
function Wpolys_KY_5_3(m,p,mu,Lminus,Q)
    BML := BasisMatrix(Lminus);
    Qminus := BML*Q*Transpose(BML);
    lat_minus := LatticeWithGram(-Qminus);
    lat_minus_d := Dual(lat_minus : Rescale := false);
    disc_group := lat_minus_d / lat_minus;
    
    Delta := -Determinant(Qminus);
    kappa := SquareFree(Delta);

    // verifications
    assert Valuation(kappa,p) in [0,1];
    assert AbelianInvariants(disc_group) eq [2,-2*kappa];

    disc_group_p := pPrimaryComponent(disc_group, p);
    if (p ne 2) then assert #disc_group_p eq p; end if;
    if (p eq 2) and IsEven(kappa) then assert #disc_group_p eq 8; end if;
    if (p eq 2) and IsOdd(kappa) then assert #disc_group_p eq 4; end if;

    K<sqrt_kappa> := QuadraticField(kappa);
    _<x> := PolynomialRing(K);

    d := Discriminant(Integers(K));
    f := Valuation(d,p);
    a := Valuation(m,p);

    assert a ge -f;

    if (a eq -f) then return 1; end if;

    /*
    norm_form := Matrix([[Norm(x+y) - Norm(x) - Norm(y) : y in [1,sqrt_kappa]] : x in [1,sqrt_kappa]]);
    eps := (p eq 2) select [1,3,5,7] else [1,Integers()!Nonsquare(GF(p))];
    can_form_Q := GramMatrix(MinkowskiReduction(LatticeWithGram(-Qminus) : Canonical));
    a := can_form_Q[1,1] / 2;
    b := can_form_Q[1,2] / 2;
    ZK := Integers(K);
    I := a*ZK + (b-sqrt_kappa)*ZK; // we only do the pricipal ideal case
    can_forms := [GramMatrix(MinkowskiReduction(LatticeWithGram(e*norm_form) : Canonical)) : e in eps];
    assert can_form_Q in can_forms; // Not implemented when this is not the case.
    e := eps[Index(can_forms, can_form_Q)];
    */
    val_e := HasseMinkowskiInvariant(lat_minus, p);
   
    return 1 + val_e*KroneckerCharacter(kappa)(m)*x^(a+f);
end function;
