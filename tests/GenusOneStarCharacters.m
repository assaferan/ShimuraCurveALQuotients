// The hypotheses of the square-coordinate route on a genus-one star quotient, and the census of
// the genus-one Atkin-Lehner quotients that supplies its cover.
//
// A Borcherds form with integral principal part has a character of finite order (Borcherds, "The
// Gross-Kohnen-Zagier theorem in higher dimensions", Duke Math. J. 97 (1999), proof of Theorem
// 4.1, with the correction in Duke Math. J. 105 (2000), 183-184 = arXiv alg-geom/9710002v3),
// equal to the order of the class of its divisor in E(Q); so Guo-Yang's Assumption 40 is
// automatic when E(Q) is torsion-free, and when the torsion has exponent e the e-th power of the
// form is a function.  ⚠ The character is one of the group uniformising the star, NOT of the
// Atkin-Lehner group, so the census below decides nothing about characters on its own.  Its role
// is in the identification argument of GenusOneStar.m: when the torsion is Z/2 on every curve of
// the isogeny class -- the star may be any of them -- and the star has exactly one genus-one
// Atkin-Lehner quotient, that quotient is the cover on which a form of divisor A - P becomes a
// function, and the signs of its square at the CM points follow.  Full rational 2-torsion breaks
// the identification, so this checks the number of nonzero rational 2-torsion points and not only
// the exponent.
//
// The genera come from the Atkin-Lehner signature (GenusShimuraCurveQuotient), the fields of
// definition from the quaternion order (FieldOfDefinitionOfCMPoint) and the torsion from the
// Cremona database; no Borcherds form, Schofer value or stored model is read.

procedure test_genus_one_star_characters()
    printf "Genus-one stars: the hypotheses of the square-coordinate route and its signs...";
    // one base per line: D, N, and whether a genus-one Atkin-Lehner quotient exists
    cases := [<39, 11, false>, <158, 1, false>, <38, 5, false>, <142, 1, false>,
              <14, 17, true>, <94, 3, false>];
    for c in cases do
        D, N, expected := Explode(c);
        W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
        Xstar := CreateShimuraQuot(D, N, W);
        Xstar`g := GenusShimuraCurveQuotient(D, N, W); Xstar`CurveID := 0;
        assert Xstar`g eq 1;
        curves := GetQuotientsAndGenera([Xstar]);
        _ := exists(star){x : x in curves | IsStarCurve(x)};
        etale := EtaleDoubleCovers(star, curves);
        assert (#etale gt 0) eq expected;
        // a double cover of a genus-one curve has genus one when etale and at least two otherwise,
        // so no index-2 quotient of a genus-one star is rational
        assert forall{x : x in [curves[i] : i in star`CoveredBy] | x`g ge 1};
    end for;

    // X_0^14(17)*: one genus-one Atkin-Lehner quotient, by {1, 2, 119, 238}
    D := 14; N := 17;
    W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
    Xstar := CreateShimuraQuot(D, N, W);
    Xstar`g := GenusShimuraCurveQuotient(D, N, W); Xstar`CurveID := 0;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){x : x in curves | IsStarCurve(x)};
    etale := EtaleDoubleCovers(star, curves);
    assert #etale eq 1;
    assert etale[1]`W eq {1, 2, 119, 238};

    // the signs of a square-coordinate at the rational CM points of the star: +1 where the
    // preimages on that quotient are rational, the sign of the discriminant of their field of
    // definition otherwise.  Every rational CM point of the star has preimages of degree at most
    // two there, so the signs are defined at all of them.
    ds := [-4, -8, -35, -43, -51, -67, -84, -168, -280, -340, -408, -427, -595];
    signs := SquareCoordinateSigns(etale[1], ds);
    // regression values (computed here): the fields of definition on X_0^14(17)/{1,2,119,238} are
    // Q at -4, -8, -43, -67, Q(sqrt -7) at -35, -84, -168, -280, -427, Q(sqrt 17) at -51, -340,
    // -408 and Q(sqrt -119) at -595
    assert signs eq [1, 1, -1, 1, 1, 1, -1, -1, -1, 1, 1, -1, -1];
    // the three real quadratic fields give +1 and the imaginary ones -1, so a sign is not merely
    // the degree of the point: both degree-two cases occur
    assert #{signs[i] : i in [1..#ds] | DegreeOfFieldOfDefinitionOfCMPoint(etale[1], ds[i]) eq 2} eq 2;

    // the arithmetic the identification rests on, over EVERY curve of the class (the star may be
    // any of them): torsion Z/2 with a single nonzero rational 2-torsion point at 238a, which is
    // what lets the Atkin-Lehner quotient be identified with the cover carrying psi; and at 429b
    // both the exponent and the number of 2-torsion points fail, so X_0^39(11)* is out of reach
    // of the route for two independent reasons.  Torsion-free classes need no route at all.
    expected := [<"238a1", 2, {1}>, <"282b1", 2, {1}>, <"429b1", 4, {1, 3}>,
                 <"142a1", 1, {0}>, <"158a1", 1, {0}>, <"190a1", 1, {0}>];
    for e in expected do
        lab, exp_e, exp_two := Explode(e);
        cls := IsogenousCurves(EllipticCurve(CremonaDatabase(), lab));
        assert LCM([Exponent(TorsionSubgroup(C)) : C in cls]) eq exp_e;
        two := {};
        for C in cls do
            T := TorsionSubgroup(C);
            Include(~two, #[Pt : Pt in T | 2*Pt eq T!0 and Pt ne T!0]);
        end for;
        assert two eq exp_two;
    end for;
    // so the route applies at 14_17 and nowhere else among the six: 282b has the right torsion but
    // no genus-one Atkin-Lehner quotient, 429b has neither
    printf "Done!\n";
end procedure;

test_genus_one_star_characters();
