// The ingredients of the square-coordinate route on a genus-one star quotient, and the census of
// the genus-one Atkin-Lehner quotients that supplies them.
//
// A Borcherds form with integral principal part has a character of finite order (Borcherds, "The
// Gross-Kohnen-Zagier theorem in higher dimensions", Duke Math. J. 97 (1999), proof of Theorem
// 4.1), equal to the order of the class of its divisor in E(Q); so Guo-Yang's Assumption 40 is
// automatic when E(Q) is torsion-free, and when the torsion has exponent e the e-th power of the
// form is a function.  ⚠ The character is one of the group uniformising the star, NOT of the
// Atkin-Lehner group, so the quotient on which a form with nontrivial character becomes a
// function need not be an Atkin-Lehner quotient: the census below decides nothing about
// characters.  Its role is to say whether an Atkin-Lehner quotient of genus one is available for
// SquareCoordinateSigns to read the fields of definition from.
//
// The genera come from the Atkin-Lehner signature (GenusShimuraCurveQuotient) and the fields of
// definition from the quaternion order (FieldOfDefinitionOfCMPoint); no Borcherds form, Schofer
// value or stored model is read.

procedure test_genus_one_star_characters()
    printf "Genus-one stars: the census of genus-one Atkin-Lehner quotients and the square-coordinate signs...";
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

    printf "Done!\n";
end procedure;

test_genus_one_star_characters();
