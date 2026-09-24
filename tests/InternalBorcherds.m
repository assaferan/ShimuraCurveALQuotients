// INTERNAL CHECKS ON THE Wpolys -- and, until 2026-09-15, A TEST THAT RAN NOTHING.
//
// ⚠ HOW THIS FILE FAILED.  It defines test_kronecker_sigma, test_bp_KY and test_W and, for as long
// as git remembers, CALLED NONE OF THEM -- nor did anything else (the apparent outside references
// to "test_W" are test_Whittaker2 and test_WeilRepresentation, prefix collisions).  The suite ran
// it in 0.000 s and reported Success, because a file of definitions asserts nothing.  Compare
// tests/Whittaker2.m and tests/WeilRepresentation.m, which invoke their procedure on the last line.
//
// ⚠ AND WHAT THAT HID.  Because nothing called them, they rotted against the library without a
// murmur.  All three were broken, in three different ways, and none of the three was mathematical:
//
//   * the import named tests/BorcherdsProducts.m, but Wpoly/Wpoly2/Wpoly_scaled live in the
//     LIBRARY, SchoferFormula.m -- so the symbols could never resolve;
//   * ShimuraCurveLattice returns a QuaternionLatticeData record where it used to return 5 values;
//   * ElementOfNorm takes the order and the basis and returns ONE value where it returned two.
//
// This is the repo's own lesson twice over: a passing check is not evidence until you know it
// could have failed, and the exempted classes are where the defects accumulate.
//
// ⚠ These three live in the LIBRARY (SchoferFormula.m), not in tests/BorcherdsProducts.m.
// The old import named the test helper and so could never resolve -- which nothing noticed,
// because the procedures below were never called.  See the header comment.
_ := ClassNumberLU(-4);   // AttachSpec is lazy and `import` compiles its file NOW, so touch
                          // one intrinsic first or references from other packages go unresolved.
import "SchoferFormula.m" : Wpoly, Wpoly2, Wpoly_scaled;

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

// ⚠ RETURNS ITS ASSERTION COUNT.  A caller that only knows "it did not throw" cannot tell a
// thorough run from an empty one; the count is what makes a silently emptied loop go red.
function test_kronecker_sigma(B)
    n := 0;
    kappas := [kappa : kappa in [1..B] | IsSquarefree(kappa)];
    for p in PrimesUpTo(B) do
        for kappa in kappas do
            assert sigmasp_Kudla_Yang_poly(p,0,kappa,true)*EulerFactor(KroneckerCharacter(kappa),p) eq 1;
            n +:= 1;
        end for;
    end for;
    return n;
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


function test_W()   // returns the number of published values checked
    // testing the few values we know from Yang
    // Two API drifts repaired 2026-09-15.  ShimuraCurveLattice returns a QuaternionLatticeData
    // record where it used to return 5 values, and Ldata`Q is that record's integral Gram matrix --
    // checked equal to the old ChangeRing(Qinv^-1, Integers()), not assumed.  ElementOfNorm now
    // takes the order and the basis and returns ONE value, not two; the call below follows
    // tests/EisensteinLocalFactors.m, which uses the current signature.
    // ⚠ Q must be INTEGRAL here, as the original line made it.  Ldata`Q holds the same matrix
    // over the RATIONALS, and `lambda_v*Q` then fails with "incompatible coefficient rings" --
    // equal as values, different as objects.
    n := 0;
    Ldata := ShimuraCurveLattice(6,1);
    Q := ChangeRing(Ldata`Qinv^(-1), Integers());
    for d in [-3,-4] do
        lambda_v := ElementOfNorm(Q, -d, Ldata`O, Ldata`basis_L);
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
            n +:= 4;
        end if;
        if d eq -3 then
            w12<x> := Wpoly_scaled(1,2,mu,Lminus,Q);
            assert w12 eq 1/2*(1-x);
            w13<x> := Wpoly_scaled(1,3,mu,Lminus,Q);
            _<sqrt3> := BaseRing(w13);
            assert w13 eq 1/sqrt3*(1+x);
            n +:= 2;
        end if;
    end for;
    return n;
end function;

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


// ---------------------------------------------------------------------------------------------
// RUN THEM.  ⚠ Count the checks: a silently skipped case must make this red, not green.
nkron := test_kronecker_sigma(20);
nW    := test_W();
assert nkron eq 104;   // 8 primes up to 20 x 13 squarefree kappa in [1..20]
assert nW eq 6;        // the six Wpoly_scaled values Yang states, at d = -4 (four) and -3 (two)
printf "InternalBorcherds: %o sigma identities, %o published Wpoly values...", nkron, nW;

// ⚠ test_bp_KY IS NOT WIRED IN, DELIBERATELY.  It is an exploratory probe, not a test: its
// assertion is commented out in the body and it returns a list of mismatches instead.  Measured
// 2026-09-15, once the import was repaired so it could run at all:
//
//     test_bp_KY(10)  ->  28 mismatches
//     test_bp_KY(20)  -> 118 mismatches
//
// ⚠ and they are not scattered: ALL 118 sit at p = 2 with mu = 0.  Every odd prime agrees, and so
// does p = 2 at mu = 1/2.
//
// ⚠⚠ THE LIBRARY IS NOT THE SIDE THAT IS WRONG.  My first reading of that -- "the discrepancy is
// confined to the Wpoly2 branch" -- was itself wrong, and is retracted.  Arbitrated 2026-09-15 with
// a brute-force representation-density count independent of BOTH sides (this repo's convention is
// Q(x) = 1/2 x G x^T, so count x G x^T = 2m mod 2^(k+1); fixed k = 14, re-confirmed at k = 18, and
// NO repeats-based stopping rule -- the plateau is long enough to defeat one):
//
//     rank 1, Q = [2 kappa], kappa in {1,2,3,5,6,7}, m in {1,2,3,4}
//        -> 24 comparisons, library vs oracle: 0 disagreements
//
// So Wpoly2 is correct at RANK 1 as well as at rank 2, where tests/Whittaker2.m now pins it on all
// five 2-adic shapes production uses.  => WHAT FAILS IS THIS FILE'S RIGHT-HAND SIDE: its
// transcription of [KY, Prop 5.1] at p = 2, the (1-x^2)/EulerFactor(...) * bp_Kudla_Yang_poly
// expression.  The header's guess at "a sqrtp factor" is in the right half of the identity but the
// wrong place -- a missing sqrtp would have moved the odd primes too.
// Left as a probe until the p = 2 form of Prop 5.1 is worked out; do not assert on it.
