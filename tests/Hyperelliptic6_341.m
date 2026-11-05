// X_0^6(341)/<w_3, w_11, w_62>, OF GENUS FOUR: its branch pair on the star is canonical up to the
// star's rational 2-torsion, by three Borcherds forms.
//
// THE SETTING. Y = X_0^6(341)/<w_3, w_11, w_62> has genus four and is a double cover of the star
// C = X_0^6(341)* (genus two), branched at the two rational points P_-24 and P_-88, the images of
// the fixed points of w_6 and w_22 (Ogg's counts: eight fixed points each, no other involution
// outside the subgroup has any). Y is hyperelliptic exactly when P_-24 + P_-88 is a fibre of the
// hyperelliptic map of C, i.e. a canonical divisor.
//
// THE ARGUMENT. A Borcherds form on C has a divisor class of finite order (Borcherds, Duke 97
// (1999), Theorem 4.1 with the Duke 105 (2000) correction), and the class is rational, so it lies
// in Jac(C)(Q)_tors. The three forms banked in data/borcherds/forms_6_341.m have divisors
//     D1 - D2,  D1 - D3,  D2 - D3,   D1 = P_-24 + P_-88,  D2 = P_-43 + P_-123,  D3 = P_-84 + P_-228,
// so the three classes [D_i - D_1] lie in a group of order at most the torsion bound t. With t = 2
// two of D1, D2, D3 are linearly equivalent; two distinct effective divisors of degree two in one
// class on a curve of genus two are both canonical. So two of the three are fibres.
//
// WHAT IS NOT YET CHECKED HERE. If Jac(C)(Q)[2] is nontrivial, D1 could still be the class K + t,
// whose unique effective representative is W1 + W2 for the two rational Weierstrass points; that
// alternative is excluded exactly by enumerating the double covers of C branched at W1 + W2
// against the trace formula, which needs the model of C (in progress). With Jac(C)(Q)[2] = 0 the
// three divisors are all fibres and Y is hyperelliptic outright.
//
// WHAT IS CHECKED. The forms are rebuilt from their eta-quotient coefficient maps, each divisor is
// recomputed from the rebuilt form (the reader refuses a mismatch), the multiple that makes each
// principal part integral is taken explicitly, and the torsion bound is recomputed from the trace
// formula with no model: #Jac(C)(F_p) = (N_1^2 + N_2)/2 - p at odd good p, into which the rational
// torsion injects.

procedure test_hyperelliptic_6_341()
    printf "X_0^6(341)/<w_3,w_11,w_62>: the branch pair is canonical up to 2-torsion...";
    D := 6; N := 341;
    W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
    Xstar := CreateShimuraQuot(D, N, W);
    Xstar`g := GenusShimuraCurveQuotient(D, N, W); Xstar`CurveID := 0;
    assert Xstar`g eq 2;
    H := {1, 3, 11, 33, 62, 186, 682, 2046};
    assert GenusShimuraCurveQuotient(D, N, H) eq 4;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};

    // the branch locus of Y -> C: fixed points of the involutions outside H, by discriminant
    outside := W diff H;
    fixed := [];
    for m in outside do
        e := NumFixedPointsByCMOrder(D, N, m);
        fixed cat:= [ <d, e[d]> : d in Keys(e) | e[d] ne 0 ];
    end for;
    assert Set(fixed) eq { <-24, 8>, <-88, 8> };

    D1 := [<-24, 1>, <-88, 1>]; D2 := [<-43, 1>, <-123, 1>]; D3 := [<-84, 1>, <-228, 1>];
    eff := [D1, D2, D3];

    // (1) the three forms, their divisors, and the integral multiples
    fs := BorcherdsFormsFromFile("data/borcherds/forms_6_341.m", star);
    assert #Keys(fs) eq 3;
    want := {};
    for i in [1..3] do for j in [i+1..3] do
        Include(~want, {<c[1], c[2]> : c in eff[i]} join {<c[1], -c[2]> : c in eff[j]});
    end for; end for;
    assert {Set(DivisorOfBorcherdsForm(fs[k], star)) : k in Keys(fs)} eq want;
    dens := [];
    for k in Sort([k : k in Keys(fs)]) do
        g := qExpansionAtoo(fs[k], 1);
        cf := [Coefficient(g, m) : m in [Valuation(g)..0]];
        n := LCM([Denominator(c) : c in cf]);
        assert forall{c : c in cf | IsIntegral(n*c)};
        assert Coefficient(g, 0) eq 0;
        Append(~dens, n);
    end for;

    // (2) the torsion bound of Jac(C)(Q), from the trace formula
    t := 0;
    for p in [q : q in PrimesInInterval(3, 60) | (D*N) mod q ne 0] do
        N1 := ComputePointsViaTrace(Xstar, p, 1); N2 := ComputePointsViaTrace(Xstar, p, 2);
        t := GCD(t, (N1^2 + N2) div 2 - p);
    end for;
    assert t eq 2;

    // (3) the counting step: three classes in a group of order at most two, so two coincide
    assert #eff gt t;
    printf " [%o divisors, torsion bound %o, principal parts scaled by %o] Done!\n", #eff, t, dens;
end procedure;

test_hyperelliptic_6_341();
