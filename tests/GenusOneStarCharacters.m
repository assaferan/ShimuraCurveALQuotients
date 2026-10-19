// The character of a Borcherds form on a genus-one star quotient, decided rather than assumed.
//
// A form psi with prescribed divisor D is a function on the star when D is principal; when D is
// not, psi is a function on X' = X/ker(chi) for an index-2 subgroup of the Atkin-Lehner group,
// and the class of D is a nonzero element of the kernel of J(X*) -> J(X'), which for a double
// cover is nontrivial only when the cover is etale -- and an etale double cover of a genus-one
// curve has genus one (a ramified one has genus at least two, by Riemann-Hurwitz).  So the forms
// on a star with no genus-one index-2 quotient are all functions, and Guo-Yang's Assumption 40
// is not needed there; see the header of GenusOneStar.m.
//
// This test checks the census that decides it on the six genus-one stars the driver has been run
// on, and the ingredients of the construction used when an etale cover does exist.  The genera
// are computed from the Atkin-Lehner signature (GenusShimuraCurveQuotient) and the fields of
// definition from the quaternion order (FieldOfDefinitionOfCMPoint); no Borcherds form, Schofer
// value or stored model is read, so this is arithmetic of the objects, not of a construction.

procedure test_genus_one_star_characters()
    printf "Characters of Borcherds forms on genus-one stars: the etale-double-cover census...";
    // one base per line: D, N, and whether an etale double cover exists
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
        // the genus-one quotients are exactly the etale double covers, and a double cover of the
        // star that is not etale has genus at least two
        assert forall{x : x in [curves[i] : i in star`CoveredBy] | x`g eq 1 or x`g ge 2};
    end for;

    // X_0^14(17)*: one etale double cover, the quotient by {1, 2, 119, 238}
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
    // preimages on the cover are rational, the sign of the discriminant of their field of
    // definition otherwise.  Every rational CM point of the star has preimages of degree at most
    // two on an etale double cover, so the signs are defined at all of them.
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
