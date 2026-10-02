// The second part of the m = 0 term of Schofer's formula: at a prime level N and a CM point of
// discriminant d off the divisor, the correction is log N * ((1/2) c_eta(0) - S(f, d)), where
// S(f, d) is the sum of the coefficients of f at oo at the exponents -k^2 Q0, Q0 = |d| or |d|/4
// (M0InfinityPoleSum; paper/level-prime-kappa.tex, Proposition prop:mult and Remark rem:xsum).
//
// Two checks.  (1) The sum itself, on hand-made series: it must pick up exactly the exponents
// k^2 Q0 and nothing else, for d = 1 mod 4 and d = 0 mod 4 alike.  (2) On X_0^15(2), the base of
// the only outside measurement of the multiplier (Guo-Yang arXiv:1510.06193v1, Table 45), the sum
// is zero for every form at every discriminant of Table 45 at which the term fires and the point
// is off the divisor; the pipeline's value there is therefore exactly (1/2) c_eta(0) log 2, as
// Table 45 requires.  Points on the divisor are excluded because the sum is allowed to be nonzero
// there (the oo-pole is what puts them on the divisor).

procedure test_pole_sum_series()
    printf "  the oo-pole sum on hand-made series...";
    R<q> := LaurentSeriesRing(Rationals());
    f := 3*q^-28 + 5*q^-7 - 2*q^-6 + q^-1 + O(q);
    assert M0InfinityPoleSum(f, -7) eq 8;     // k = 1 (q^-7) and k = 2 (q^-28)
    assert M0InfinityPoleSum(f, -28) eq 8;    // Q0 = 28/4 = 7 again: the same two exponents
    assert M0InfinityPoleSum(f, -3) eq 0;     // exponents 3, 12, 27: none present
    assert M0InfinityPoleSum(f, -24) eq -2;   // Q0 = 6: q^-6
    assert M0InfinityPoleSum(f, -4) eq 1;     // Q0 = 1: q^-1 (q^-4, q^-9, ... absent)
    assert M0InfinityPoleSum(2*q^3 + O(q^5), -7) eq 0;   // holomorphic at oo
    printf " ok\n";
end procedure;

procedure test_pole_sum_15_2()
    printf "  the oo-pole sum on X0^15(2) at the Table 45 discriminants...";
    D := 15; N := 2;
    Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
    Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};
    fs := BorcherdsForms(star, curves : Prec := 100);
    table45 := [-7, -15, -60, -52, -88, -120, -132, -148, -168, -228, -232, -280, -312, -340,
                -372, -408, -520, -708, -760];
    firing := [d : d in table45 | GCD(N, FundamentalDiscriminant(d)) ne N];
    assert firing eq [-7, -15, -60];
    checked := 0;
    for k in Keys(fs) do
        f := fs[k];
        foo := qExpansionAtoo(f, 1);
        on_divisor := {pair[1] : pair in DivisorOfBorcherdsForm(f, star)};
        for d in firing do
            if d in on_divisor then continue; end if;
            assert M0InfinityPoleSum(foo, d) eq 0;
            checked +:= 1;
        end for;
    end for;
    // 9 forms x 3 discriminants, less the four forms with c_oo(-15) != 0, whose divisor contains
    // both tau_{-15} and tau_{-60} (d = 4m/r^2 with m = -15, r = 2 and 1)
    assert checked eq 19;
    printf " ok (%o pairs)\n", checked;
end procedure;

printf "Testing the oo-pole part of the m = 0 term...\n";
test_pole_sum_series();
test_pole_sum_15_2();
printf "Done!\n";
