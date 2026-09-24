// The exact m = 0 multipliers, full panel: M0MultiplierExact against the measured ground truth
// on X0^15(2) -- all NINE forms, not just the two the old closed form could score.
//
// The expected values are the independently measured ground truth (gtsweep + Hauptmodul
// consistency; memory m0-ground-truth-15-2), re-derived exactly by the cusp-class assembly
// (vvdata/weyl-campaign: the 15_2 full panel 9/9, zeros to 1e-71).  The old M0Multiplier
// closed form (tests/M0Multiplier.m) pins form -1's value 4; this test pins the rest.

procedure test_m0_multiplier_exact_15_2()
    printf "  M0MultiplierExact on X0^15(2), full panel...";
    D := 15; N := 2;
    Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
    Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};
    fs := BorcherdsForms(star, curves : Prec := 100);
    keys := Sort([k : k in Keys(fs)]);
    assert keys eq [-2, -1] cat [9..15];

    Ld := ShimuraCurveLattice(D, N);
    mults := M0MultiplierExact([fs[k] : k in keys], Ld, D, N);

    expected := AssociativeArray();
    expected[-2] := 2;  expected[-1] := 4;
    expected[9]  := 0;  expected[10] := 0;  expected[11] := 4;
    expected[12] := 2;  expected[13] := 4;  expected[14] := -2;
    expected[15] := 2;
    for i->k in keys do
        assert mults[i] eq expected[k];
    end for;
    printf " ok\n";
end procedure;

test_m0_multiplier_exact_15_2();
