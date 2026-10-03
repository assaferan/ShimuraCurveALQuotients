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

// (3) Values at NON-FUNDAMENTAL discriminants against Guo-Yang's Table 45 (arXiv:1510.06193v1),
// which gives the Hauptmodul s of X_0^15(2)/W at each CM point.  The Borcherds product the code
// returns as fs[-2] has |value| = (1280/9) |s(s-2)| at every single-point cycle: this is first
// CHECKED at the fundamental discriminants -7, -15, -52 (s = 1/4, 5/4, 1), so the identification is
// not assumed, and then REQUIRED at the conductor-2 points -28, -60 (s = 9/4, -1/12) and at the
// conductor-4 point -240 (s = -25/12), where 2 splits in Q(sqrt d).  Before 2026-10-03 the value at
// -240 came out exactly halved (and irrational): the m > 0 sum used the class number of the order
// where the formula takes the field's.  -48 (conductor 4, 2 INERT in Q(sqrt -3), s = -1/4) is the
// one point where the m = 0 term must NOT fire: the zero coset's local factor at an inert
// conductor plane has no pole (campaign level-p2/inertplane.m), and with the term the value came
// out 4 log 2 too large.
procedure test_values_at_conductor_discriminants_15_2()
    printf "  values at conductor-2 and conductor-4 discriminants on X0^15(2) against Table 45...";
    D := 15; N := 2;
    Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
    Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};
    fs := BorcherdsForms(star, curves : Prec := 100);
    svals := AssociativeArray();                      // Guo-Yang v1 Table 45, column s
    svals[-7] := 1/4; svals[-15] := 5/4; svals[-52] := 1;
    svals[-28] := 9/4; svals[-60] := -1/12; svals[-240] := -25/12; svals[-48] := -1/4;
    want := Set(Keys(svals));
    cm := CandidateDiscriminants(star, curves : Keep := want);
    rat := cm[1]; quad := cm[2];
    for t in [<-28, 2, 1>, <-240, 4, 1>, <-48, 4, 1>] do   // <d, conductor, degree>; conductor 4 is never offered
        if not exists{u : u in rat | u[1] eq t[1]} then Append(~rat, t); end if;
    end for;
    tab, _ := AbsoluteValuesAtCMPoints(star, curves, [rat, quad], fs : MaxNum := 60, Prec := 100, Exclude := {}, Include := want);
    ks := Sort([k : k in Keys(fs)]);
    assert ks[1] eq -2;
    // the value prints as a formal sum "aLog2+bLog3..."; compare against the expected one in that form
    function logstring(q)   // q a positive rational -> "aLog2+bLog3..." in the code's format
        f := Factorization(Numerator(q)); g := Factorization(Denominator(q));
        terms := Sort([<t[1], t[2]> : t in f] cat [<t[1], -t[2]> : t in g]);
        s := "";
        for t in terms do
            s cat:= (t[2] gt 0 and #s gt 0 select "+" else "") cat (t[2] eq 1 select "" else (t[2] eq -1 select "-" else IntegerToString(t[2]))) cat "Log" cat IntegerToString(t[1]);
        end for;
        return s;
    end function;
    C := 1280/9;
    nchecked := 0;
    for d in [-7, -15, -52, -28, -60, -240, -48] do
        i := Index(tab`Discs, d);
        error if i eq 0, Sprintf("d = %o was not evaluated", d);
        got := Sprint(tab`Values[1][i]);
        exp := logstring(C * AbsoluteValue(svals[d] * (svals[d] - 2)));
        error if got ne exp,
            Sprintf("X0^15(2), d = %o: fs[-2] value is %o, Table 45 gives s = %o hence %o", d, got, svals[d], exp);
        nchecked +:= 1;
    end for;
    assert nchecked eq 7;
    printf " ok (%o discriminants, conductors 1, 2 and 4)\n", nchecked;
end procedure;

printf "Testing the oo-pole part of the m = 0 term...\n";
test_pole_sum_series();
test_pole_sum_15_2();
test_values_at_conductor_discriminants_15_2();
printf "Done!\n";
