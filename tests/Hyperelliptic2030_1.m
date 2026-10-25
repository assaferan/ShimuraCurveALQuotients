// X_0^2030(1)/<w_10, w_14, w_58>, OF GENUS FOUR, IS HYPERELLIPTIC, by three Borcherds forms on its star.
//
// THE SETTING. Y = X_0^2030(1)/<w_10, w_14, w_58> has genus four and is a double cover of the star
// C = X_0^2030(1)* (genus two), branched at the two rational points P_-8 and P_-280, the images of
// the fixed points of w_2 and w_70, the only involutions outside the subgroup with fixed points
// (Ogg's counts: eight each). Y is hyperelliptic exactly when P_-8 + P_-280 is a fibre of the
// hyperelliptic map of C, that is a canonical divisor: the hyperelliptic involution of Y, being
// central, descends to that of C and swaps the two fixed points of the covering involution, and
// conversely a canonical branch divisor gives the cover y^2 = q(t^2/c) with C : y^2 = q(u).
//
// THE ARGUMENT. A Borcherds form on C has a divisor class of finite order (Borcherds, Duke 97
// (1999), Theorem 4.1 with the Duke 105 (2000) correction), and the class is rational, so it lies
// in Jac(C)(Q)_tors. The three forms banked in data/borcherds/forms_2030_1.m have divisors
//     D1 - D2,  D1 - D3,  D2 - D3,   D1 = P_-8 + P_-280,  D2 = P_-43 + P_-163,  D3 = P_-148 + P_-532,
// and their principal parts are integral as they stand. The torsion of Jac(C)(Q) injects into
// Jac(C)(F_p) for odd good p and the gcd of those orders is 1, so D1 ~ D2 ~ D3 exactly; two
// distinct effective divisors of degree two in one class on a curve of genus two are both
// canonical. So D1 is a fibre and Y is hyperelliptic, with nothing assumed and nothing numerical.
//
// WHAT IS CHECKED. The branch locus from the fixed-point counts; the forms rebuilt from their
// eta-quotient coefficient maps, each divisor recomputed from the rebuilt form (the reader refuses
// a mismatch), the integrality of each principal part; and the torsion bound from the trace
// formula, #Jac(C)(F_p) = (N_1^2 + N_2)/2 - p, with no model of either curve.

procedure test_hyperelliptic_2030_1()
    printf "X_0^2030(1)/<w_10,w_14,w_58> is hyperelliptic...";
    D := 2030; N := 1;
    W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
    Xstar := CreateShimuraQuot(D, N, W);
    Xstar`g := GenusShimuraCurveQuotient(D, N, W); Xstar`CurveID := 0;
    assert Xstar`g eq 2;
    H := {1, 10, 14, 35, 58, 145, 203, 2030};
    assert GenusShimuraCurveQuotient(D, N, H) eq 4;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};

    // the branch locus of Y -> C: fixed points of the involutions outside H, by discriminant
    fixed := [];
    for m in W diff H do
        e := NumFixedPointsByCMOrder(D, N, m);
        fixed cat:= [ <d, e[d]> : d in Keys(e) | e[d] ne 0 ];
    end for;
    assert Set(fixed) eq { <-8, 8>, <-280, 8> };

    D1 := [<-8, 1>, <-280, 1>]; D2 := [<-43, 1>, <-163, 1>]; D3 := [<-148, 1>, <-532, 1>];
    eff := [D1, D2, D3];

    // (1) the three forms, their divisors, and the integrality of the principal parts
    fs := BorcherdsFormsFromFile("data/borcherds/forms_2030_1.m", star);
    assert #Keys(fs) eq 3;
    want := {};
    for i in [1..3] do for j in [i+1..3] do
        Include(~want, {<c[1], c[2]> : c in eff[i]} join {<c[1], -c[2]> : c in eff[j]});
    end for; end for;
    assert {Set(DivisorOfBorcherdsForm(fs[k], star)) : k in Keys(fs)} eq want;
    for k in Keys(fs) do
        g := qExpansionAtoo(fs[k], 1);
        assert forall{m : m in [Valuation(g)..0] | IsIntegral(Coefficient(g, m))};
        assert Coefficient(g, 0) eq 0;
    end for;

    // (2) Jac(C)(Q) is torsion-free, from the trace formula
    t := 0;
    for p in [q : q in PrimesInInterval(3, 60) | (D*N) mod q ne 0] do
        N1 := ComputePointsViaTrace(Xstar, p, 1); N2 := ComputePointsViaTrace(Xstar, p, 2);
        t := GCD(t, (N1^2 + N2) div 2 - p);
    end for;
    assert t eq 1;

    // (3) so the three divisors are linearly equivalent, hence canonical, and Y is hyperelliptic
    printf " [3 divisors pairwise related, torsion bound %o] Done!\n", t;
end procedure;

test_hyperelliptic_2030_1();
