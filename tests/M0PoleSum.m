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
// which gives the Hauptmodul s of X_0^15(2)/W at each CM point.  Every Borcherds form of the model
// set is a function on the genus-0 star curve with known divisor (DivisorOfBorcherdsForm, Guo-Yang
// Lemma 25), so its value at a CM point is C * prod |s(d) - s(d_i)|^(m_i) over the zeros and poles
// d_i of the divisor other than tau_-12 (where s = oo).  The constant C of each form is read at the
// fundamental discriminant -7 and the formula is then REQUIRED at the fundamental points -15, -52,
// the conductor-2 points -28, -60 (2 splits) and the conductor-4 points -240 (2 splits) and -48
// (2 inert in Q(sqrt -3)) -- 46 values, 18 of them at conductor 4.  At conductor 4 the m = 0 term
// is the fibre sum M0FibreCorrection (paper/kappa0-proof-standalone.tex, prop:fibre, lem:Wcond):
// (m + c_oo(-15)) log 2 at -240 and (2/3 (m + c_oo(-3)) + 1/3 b) log 2 at -48, with m the
// multiplier and b the cusp-0 coefficient of q^(-3/4).  Before 2026-10-03 the term was applied as
// a rule ("fire iff 2 splits"), right at -240 and -48 only for the forms with c_oo(-15) = 0,
// resp. m + c_oo(-3) = 0 and b = 0 -- which included the one form this test then checked, the
// W = {1,3,5,15} cover with value (1280/9)|s(s-2)|, kept below as check (3a); and the m > 0 sum
// used the class number of the order where the formula takes the field's (every conductor-4 value
// halved).
// (Until 2026-10-03 this test called its form fs[-2]: the table's rows follow Keys(fs), not the
// sorted keys, and row 1 is key 11.)
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
    svals[-40] := 0; svals[-120] := 2;                // zeros of the forms; -12 is their common pole
    points := [-7, -15, -52, -28, -60, -240, -48];
    cm := CandidateDiscriminants(star, curves : Keep := Set(points));
    rat := cm[1]; quad := cm[2];
    for t in [<-28, 2, 1>, <-240, 4, 1>, <-48, 4, 1>] do   // <d, conductor, degree>; conductor 4 is never offered
        if not exists{u : u in rat | u[1] eq t[1]} then Append(~rat, t); end if;
    end for;
    tab, _ := AbsoluteValuesAtCMPoints(star, curves, [rat, quad], fs : MaxNum := 60, Prec := 100, Exclude := {}, Include := Set(points));
    for d in points do error if Index(tab`Discs, d) eq 0, Sprintf("d = %o was not evaluated", d); end for;
    // a table cell is a formal sum of logs of primes; the forms' values here are rational numbers
    function cell_value(x)
        error if x eq LogSum(0) or x eq LogSum(Infinity()), "a divisor point was not skipped";
        v := Rationals()!1;
        for q in Keys(x`log_coeffs) do
            c := x`log_coeffs[q];
            error if not IsIntegral(c), Sprintf("irrational value %o", x);
            v *:= (Rationals()!q)^(Integers()!c);
        end for;
        return v;
    end function;

    // (3a) the W = {1,3,5,15} cover: divisor (-40) + (-120) - 2(-12), value (1280/9) |s(s-2)|
    row := Index(tab`Keys_fs, 11);
    assert row gt 0;
    C := 1280/9;
    for d in points do
        got := cell_value(tab`Values[row][Index(tab`Discs, d)]);
        exp := C * AbsoluteValue(svals[d] * (svals[d] - 2));
        error if got ne exp,
            Sprintf("X0^15(2), d = %o: the W = {1,3,5,15} form's value is %o, Table 45 gives s = %o hence %o", d, got, svals[d], exp);
    end for;

    // (3b) every form, from its divisor: C read at -7, required everywhere else off the divisor
    nchecked := 0;
    for r->k in tab`Keys_fs do
        div_f := DivisorOfBorcherdsForm(fs[k], star);
        for pr in div_f do
            error if pr[1] ne -12 and not IsDefined(svals, pr[1]),
                Sprintf("form %o has a divisor point of discriminant %o not in the table", k, pr[1]);
        end for;
        model := func< d | &*[Rationals() | AbsoluteValue(svals[d] - svals[pr[1]])^(Integers()!pr[2]) : pr in div_f | pr[1] ne -12] >;
        on_div := {pr[1] : pr in div_f};
        assert -7 notin on_div;
        Ck := cell_value(tab`Values[r][Index(tab`Discs, -7)]) / model(-7);
        for d in points do
            if d eq -7 or d in on_div then continue; end if;
            got := cell_value(tab`Values[r][Index(tab`Discs, d)]);
            exp := Ck * model(d);
            error if got ne exp,
                Sprintf("X0^15(2), form %o (divisor %o), d = %o: value %o, Table 45 and the divisor give %o",
                        k, div_f, d, got, exp);
            nchecked +:= 1;
        end for;
    end for;
    assert nchecked eq 46;
    printf " ok (%o values of 9 forms at conductors 1, 2 and 4)\n", nchecked;

    // (3c) a conductor prime OUTSIDE the level: d = -588 = 14^2 * (-3), where 7 does not divide 30.
    // Table 45 is silent, but the nine values (now sums over the three star points of discriminant
    // -588) give N(s), N(s - 2) and N((s + 1/12)(s - 5/4)) through the divisors, and the monic cubic
    // H in Q[X] with those absolute values at 0, 2 and (as a product) at -1/12, 5/4 -- the class
    // polynomial of s at d = -588 -- must exist: a rationality condition on the discriminant of a
    // quadratic.  With the log 7 term of M0FibreCorrection (standalone lem:unimod) it does, and the
    // cubic field it cuts out must be the field of definition of the CM points, which
    // FieldsOfDefinitionOfCMPoint computes by class field theory (discriminant -588); without the
    // term no rational cubic exists (campaign level-p2/classpoly.log).
    printf "  the class polynomial of s at d = -588 (conductor prime 7 outside the level)...";
    d := -588;
    OK := MaximalOrder(QuadraticField(d)); O := sub<OK | 14>;
    npts := NumberOfOptimalEmbeddings(O, D, N) div 8;
    assert npts eq 3;
    Ld := ShimuraCurveLattice(D, N);
    vals := SchoferFormula([fs[k] : k in tab`Keys_fs], d, D, N, Ld : PointDegree := npts);
    norms := AssociativeArray();
    for r->k in tab`Keys_fs do
        div_f := DivisorOfBorcherdsForm(fs[k], star);
        model7 := &*[Rationals() | AbsoluteValue(svals[-7] - svals[pr[1]])^(Integers()!pr[2]) : pr in div_f | pr[1] ne -12];
        Ck := cell_value(tab`Values[r][Index(tab`Discs, -7)]) / model7;
        norms[k] := cell_value(vals[r]) / Ck^npts;          // prod over the divisor of N|s - s_i|^{m_i}
    end for;
    A := norms[-1]; B := norms[-2]; CE := norms[14];          // N|s|, N|s-2|, N|(s+1/12)(s-5/4)|
    assert norms[11] eq A*B and norms[9] eq B*CE and norms[10] eq A*CE and norms[12] eq A*B*CE;
    P<X> := PolynomialRing(Rationals());
    found := [];
    for e0, e1, e2 in [1, -1] do
        // H = X^3 + a X^2 + b X + c with H(0) = e0 A, H(2) = e1 B: c and b = (e1 B - 8 - 4a - c)/2 are
        // determined by a, and H(-1/12) H(5/4) = e2 CE is a quadratic in a
        c := e0*A;
        Pa<a> := PolynomialRing(Rationals());
        b := (e1*B - 8 - 4*a - c)/2;
        quad := (-1/1728 + a/144 - b/12 + c)*(125/64 + 25*a/16 + 5*b/4 + c) - e2*CE;
        for root in Roots(quad) do
            aa := root[1]; bb := Evaluate(b, aa);
            Append(~found, X^3 + aa*X^2 + bb*X + c);
        end for;
    end for;
    error if IsEmpty(found), "no rational monic cubic H has the norms the values prescribe at d = -588";
    flds := FieldsOfDefinitionOfCMPoint(star, d);
    want_disc := {Discriminant(MaximalOrder(F)) : F in flds | Type(F) ne FldRat};
    assert want_disc eq {-588};
    good := [H : H in found | IsIrreducible(H) and Discriminant(MaximalOrder(NumberField(H))) eq -588];
    error if #good ne 1,
        Sprintf("expected exactly one cubic cutting out the field of definition (discriminant -588), found %o among %o", #good, found);
    // regression pin (computed by this test on 2026-10-03, not an independent value): the roots have
    // 7-adic valuation -1/3, i.e. all three points reduce to the pole tau_-12 at the prime above 7
    assert good[1] eq X^3 - 191/54*X^2 + 343/432*X - 6889/48384;
    printf " ok\n";
end procedure;

printf "Testing the oo-pole part of the m = 0 term...\n";
test_pole_sum_series();
test_pole_sum_15_2();
test_values_at_conductor_discriminants_15_2();
printf "Done!\n";
