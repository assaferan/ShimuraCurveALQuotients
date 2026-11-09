// X_0^6(403)/<w_2, w_13, w_93>, OF GENUS FOUR, IS HYPERELLIPTIC, by three Borcherds forms on its star.
//
// THE SETTING. Y = X_0^6(403)/<w_2, w_13, w_93> has genus four and is a double cover of the star
// C = X_0^6(403)* (genus two), branched at the two rational points P_-3 and P_-403, the images of
// the fixed points of w_3 and w_403, the only involutions outside the subgroup with fixed points
// (Ogg's counts: eight each). Y is hyperelliptic exactly when P_-3 + P_-403 is a fibre of the
// hyperelliptic map of C, i.e. a canonical divisor (see tests/Hyperelliptic2030_1.m).
//
// The eight points of discriminant -3 are the elliptic points of order three of X_0^6(403); their
// image on C is one point, and the divisors below are divisors on the coarse curve C, the
// stabiliser orders having been divided out when the divisor of a form is read off.
//
// THE ARGUMENT. The three forms banked in data/borcherds/forms_6_403.m have divisors D1 - D2,
// D1 - D3, D2 - D3 for D1 = P_-3 + P_-403, D2 = P_-43 + P_-91, D3 = P_-88 + P_-520, with integral
// principal parts, so the three classes agree up to torsion; the torsion of Jac(C)(Q) injects into
// Jac(C)(F_p) at odd good p and the gcd of those orders is 1, so D1 ~ D2 ~ D3, and two distinct
// linearly equivalent effective divisors of degree two on a curve of genus two are canonical. So
// D1 is a fibre and Y is hyperelliptic.

procedure test_hyperelliptic_6_403()
    printf "X_0^6(403)/<w_2,w_13,w_93> is hyperelliptic...";
    D := 6; N := 403;
    W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
    Xstar := CreateShimuraQuot(D, N, W);
    Xstar`g := GenusShimuraCurveQuotient(D, N, W); Xstar`CurveID := 0;
    assert Xstar`g eq 2;
    H := {1, 2, 13, 26, 93, 186, 1209, 2418};
    assert GenusShimuraCurveQuotient(D, N, H) eq 4;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};

    fixed := [];
    for m in W diff H do
        e := NumFixedPointsByCMOrder(D, N, m);
        fixed cat:= [ <d, e[d]> : d in Keys(e) | e[d] ne 0 ];
    end for;
    assert Set(fixed) eq { <-3, 8>, <-403, 8> };

    D1 := [<-3, 1>, <-403, 1>]; D2 := [<-43, 1>, <-91, 1>]; D3 := [<-88, 1>, <-520, 1>];
    eff := [D1, D2, D3];

    fs := BorcherdsFormsFromFile("data/borcherds/forms_6_403.m", star);
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

    t := 0;
    for p in [q : q in PrimesInInterval(3, 60) | (D*N) mod q ne 0] do
        N1 := ComputePointsViaTrace(Xstar, p, 1); N2 := ComputePointsViaTrace(Xstar, p, 2);
        t := GCD(t, (N1^2 + N2) div 2 - p);
    end for;
    assert t eq 1;
    printf " [3 divisors pairwise related, torsion bound %o] Done!\n", t;
end procedure;

test_hyperelliptic_6_403();
