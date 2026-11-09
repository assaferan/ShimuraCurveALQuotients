// X_0^46(31)* IS HYPERELLIPTIC, with no appeal to Guo-Yang's Assumption 40.
//
// THE ARGUMENT, which is general and is the point of this test. Let X/Q be a curve with a rational
// point whose Jacobian has rational torsion of order at most t. Suppose D_1, ..., D_(t+1) are
// distinct effective divisors of degree two on X, some positive integer multiple of each
// difference D_i - D_j being the divisor of a Borcherds form with integral principal part. Such a
// form has a divisor class of finite order
// (Borcherds, "The Gross-Kohnen-Zagier theorem in higher dimensions", Duke Math. J. 97 (1999),
// proof of Theorem 4.1, with the correction in Duke Math. J. 105 (2000), 183-184, which restores
// the finite-order condition the printed version omits), and the class is Galois-stable, so it
// lies in J(Q)_tors -- and a class some multiple of which is torsion is itself torsion, J(Q) being
// finitely generated, which is what admits the multiple. The t+1 classes [D_i - D_1] therefore lie
// in a group of order at most t, so
// two coincide, and that difference is PRINCIPAL. A function of degree two exists and X, of genus
// at least two, is hyperelliptic. Assumption 40 is not needed: it asks that a particular class be
// trivial, where counting gives that one of several is.
//
// THE CASE. X_0^46(31)* has genus three and no Atkin-Lehner involution left, so none of the cover
// routes reaches it. A search over the 414 degree-two divisors supported on its rational CM points
// and quadratic cycles (2 h 35 m, tools/hyperellipticstar.m, stage 1) returned exactly six, and
// they are the six pairwise differences of
//     A = P_-3 + P_-403,  B = P_-24 + P_-232,  C = P_-115 + P_-372,  D = Z(-292),
// the last a cycle of degree two. Four divisors against a torsion bound of two is the hypothesis
// above with room to spare. The six pole orders at infinity run from 556 to 1760, far beyond the
// window a search over divisors with a single double pole can reach, which is why this needed the
// general divisor shape; and three of the six use the quadratic cycle.
//
// WHAT IS CHECKED HERE, and what is provenance. The search is not repeated: the six forms are
// banked in data/borcherds/forms_46_31.m as their combinations of eta quotients, which is what
// witnesses that forms with these principal parts exist at all, and the test rebuilds them from
// that and recomputes their divisors. The torsion bound is recomputed from the trace formula, with
// no model of the curve. So every hypothesis of the argument is checked here; only the exhaustion
// of the 414 candidates -- which the argument does not use -- is provenance.

procedure test_hyperelliptic_46_31()
    printf "X_0^46(31)* is hyperelliptic, by counting rather than Assumption 40...";
    D := 46; N := 31;
    W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
    Xstar := CreateShimuraQuot(D, N, W);
    Xstar`g := GenusShimuraCurveQuotient(D, N, W); Xstar`CurveID := 0;
    assert Xstar`g eq 3;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};
    // no Atkin-Lehner involution is left, so the cover routes do not apply
    assert IsEmpty(star`CoveredBy);

    A := [<-3, 1>, <-403, 1>];  B := [<-24, 1>, <-232, 1>];
    C := [<-115, 1>, <-372, 1>]; E := [<-292, 1>];
    eff := [A, B, C, E];

    // (1) The six forms, rebuilt from their eta-quotient combinations. Rebuilding is what says a
    // form with that principal part exists; BorcherdsFormsFromFile recomputes each divisor from
    // the rebuilt form and refuses the file if one disagrees.
    fs := BorcherdsFormsFromFile("data/borcherds/forms_46_31.m", star);
    assert #Keys(fs) eq 6;
    // their divisors are exactly the six pairwise differences of the four effective divisors
    want := {};
    for i in [1..4] do
        for j in [i+1..4] do
            Include(~want, {<c[1], c[2]> : c in eff[i]} join {<c[1], -c[2]> : c in eff[j]});
        end for;
    end for;
    assert {Set(DivisorOfBorcherdsForm(fs[k], star)) : k in Keys(fs)} eq want;
    // What makes the class torsion is that an integer multiple of each form is a genuine Borcherds
    // input.  The principal parts here are rational, not integral: the search's integrality option
    // makes the solution integral over an echelonised basis of the pool, which does not make the
    // q-expansion integral.  So the multiple is taken explicitly -- n f has an integral principal
    // part, its lift is the n-th power of the lift of f and so has n times the divisor, and a
    // class whose n-th multiple has finite order has finite order itself.
    dens := [];
    for k in Sort([k : k in Keys(fs)]) do
        g := qExpansionAtoo(fs[k], 1);
        cf := [Coefficient(g, m) : m in [Valuation(g)..0]];
        n := LCM([Denominator(c) : c in cf]);
        assert forall{c : c in cf | IsIntegral(n*c)};
        assert Coefficient(g, 0) eq 0;
        Append(~dens, n);
    end for;

    // (2) The torsion bound, from the trace formula and no model.
    t := JacobianTorsionBound(Xstar : MaxPrime := 60, NumPrimes := 6);
    assert t eq 2;

    // (3) The counting step. More effective divisors than the torsion bound allows classes, so two
    // of them are linearly equivalent and the curve is hyperelliptic. This is the assertion that
    // fails if either hypothesis ever stops holding -- a larger torsion bound, or fewer divisors.
    assert #eff gt t;
    printf " [%o divisors, torsion bound %o, principal parts scaled by %o] Done!\n", #eff, t, dens;
end procedure;

test_hyperelliptic_46_31();
