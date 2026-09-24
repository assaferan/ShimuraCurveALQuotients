// Regression test for the CM-point supply on even-level bases.
//
// ⚠ CONTEXT CHANGED 2026-09-07: the coprime-to-level filter is now **OFF BY DEFAULT**
// (`CMCOPRIME=1` re-enables it). Everything below describes why `Keep` exists and still pins its
// behaviour, but at this base the anchors now arrive without it -- see check [4].
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

    // [4] ⚠ REWRITTEN 2026-09-07, WHEN THE FILTER BECAME OFF BY DEFAULT.
    //     This used to assert the OPPOSITE -- `IsEmpty(anchors meet base_rat)` -- to show that the
    //     coprime-to-level filter was what hid the anchors, so the base genuinely needed `Keep`.
    //     That premise was the OLD DEFAULT. With the filter off, the anchors are offered by the
    //     plain call, which is exactly what the flip was for, and the old assertion necessarily
    //     fails. It failed in CI on 206a0cb3; this is a changed premise, NOT a regression.
    //     What the pipeline actually depends on is that the anchors ARE AVAILABLE, so assert that
    //     directly -- a stronger statement than the old one, and now true without `Keep`.
    //     ⚠ AND IT MUST HOLD IN BOTH REGIMES. Asserting only the new default made the test fail
    //     under `CMCOPRIME=1` -- i.e. it then described one mode and broke in the other, which is
    //     how this assertion got stale in the first place. So branch on the flag and state the
    //     invariant for each.
    if GetEnv("CMCOPRIME") ne "" then
        // filter ON (the escape hatch): the anchors are hidden, which is exactly why Keep exists.
        assert IsEmpty(anchors meet base_rat);
    else
        // filter OFF (the default since 2026-09-07): the anchors arrive from the plain call.
        assert anchors subset base_rat;
    end if;
    assert &and[IsEven(d) : d in anchors];
    //     ⚠ `Keep` is therefore a NO-OP at this base under the default, but it is NOT dead: it is
    //     what makes the anchors survive when the filter is re-enabled with `CMCOPRIME=1`. Checks
    //     [1]-[3] still pin that behaviour, and they are what to look at if that escape hatch is
    //     ever used.

    printf "Done!\n";
end procedure;

test_CMSupply();
