// Regression test for the CM-point supply on even-level bases.
//
// CandidateDiscriminants applies a coprime-to-level filter to the CM points it offers Schofer's
// formula. The filter is a blunt instrument: on an even-level base it drops EVERY even discriminant,
// including the zeros and poles of the two hauptmoduls, which the pipeline requires as anchors
// (AbsoluteValuesAtCMPoints takes them as its Include / must-use set). Because that Include set is
// selected FROM the candidate list, dropping them upstream silently emptied it, and the build then
// failed downstream in ValuesAtCMPoints ("Sequence index 0 should be in the range 1 to 4") because no
// surviving discriminant was a zero of the hauptmodul.
//
// X0^15(2) is the witness: its four divisor discriminants are -40, -12, -120, -12 -- all even, so all
// dropped -- leaving only 2 rational candidates against a demand of 7.
//
// The fix is the Keep parameter: discriminants listed there are exempt from the filter, so the anchors
// are admitted without relaxing the filter globally. This test pins that behaviour directly, without
// building Borcherds forms (which SchoferIsometry.m / the X0_* pipeline tests cover end-to-end).

procedure test_CMSupply()
    printf "Testing CM-point supply (hauptmodul anchors) on X0^15(2)...";
    D := 15; N := 2;
    Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
    Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W);
    Xstar`CurveID := 0;
    curves := GetQuotientsAndGenera([Xstar]);
    assert exists(star){c : c in curves | IsStarCurve(c)};

    // The divisor discriminants of the two hauptmoduls of X0^15(2) (DivisorOfBorcherdsForm on
    // fs[-1], fs[-2]); hard-coded so the test needs no Borcherds-form computation.
    anchors := {-40, -12, -120};

    base := CandidateDiscriminants(star, curves);
    base_rat := {p[1] : p in base[1]};

    kept := CandidateDiscriminants(star, curves : Keep := anchors);
    kept_rat := {p[1] : p in kept[1]};

    // [1] Every anchor is offered as a candidate once it is pinned via Keep.  This is the property the
    //     pipeline depends on; before the fix the intersection was empty.
    assert anchors subset kept_rat;

    // [2] Keep only ADDS: it must never drop a discriminant the unfiltered call already offered.
    assert base_rat subset kept_rat;

    // [3] Keep is targeted, not a global relaxation: nothing beyond the anchors is admitted.
    assert kept_rat subset (base_rat join anchors);

    // [4] The filter really is what was hiding them -- i.e. this base genuinely needs the fix.
    //     (All four X0^15(2) divisor discriminants are even, hence dropped by a coprime-to-N=2 filter.)
    assert IsEmpty(anchors meet base_rat);
    assert &and[IsEven(d) : d in anchors];

    printf "Done!\n";
end procedure;

test_CMSupply();
