// The m = 0 multipliers by exact algebraic arithmetic (M0MultipliersAlgebraic) against the
// measured ground truth on X0^15(2), all nine forms -- the same expected values as
// tests/M0MultiplierExact.m (memory m0-ground-truth-15-2: gtsweep + Hauptmodul consistency) -- and
// against the numerical routine M0MultipliersBySupport on the same forms.  The algebraic routine
// shares no step with the numerical one: no Fourier transform, no eta evaluation, no rational
// snap; agreement of the two on every form is the check that both are right.

procedure test_m0_multipliers_algebraic_15_2()
    printf "  M0MultipliersAlgebraic on X0^15(2), full panel...";
    D := 15; N := 2;
    Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
    Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};
    fs := BorcherdsForms(star, curves : Prec := 100);
    keys := Sort([k : k in Keys(fs)]);
    assert keys eq [-2, -1] cat [9..15];

    Ld := ShimuraCurveLattice(D, N);
    t := Realtime();
    arrs := M0MultipliersAlgebraic([fs[k] : k in keys], Ld, D, N);
    t := Realtime(t);

    expected := AssociativeArray();
    expected[-2] := 2;  expected[-1] := 4;
    expected[9]  := 0;  expected[10] := 0;  expected[11] := 4;
    expected[12] := 2;  expected[13] := 4;  expected[14] := -2;
    expected[15] := 2;
    for i->k in keys do
        assert Keys(arrs[i]) eq {N} and arrs[i][N] eq expected[k];
    end for;
    num := M0MultipliersBySupport([fs[k] : k in keys], Ld, D, N);
    for i->k in keys do
        assert num[i][N] eq arrs[i][N];
    end for;
    printf " ok (%o s)\n", t;
end procedure;

test_m0_multipliers_algebraic_15_2();
