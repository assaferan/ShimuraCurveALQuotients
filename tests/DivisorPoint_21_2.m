// A point ON the divisor, second base: X_0^21(2) at tau_{-4}, the pole of the Hauptmodul s (Guo-Yang,
// arXiv:1510.06193v1, appendix table "CM-values of X_0^21(2)": s(-4) = oo).  Every form has a pole
// there, so single values are infinite; a quotient f_a^deg(b) / f_b^deg(a) is finite and the divisors
// force it: with |f_k| = C_k prod |s - s_i|^{m_i} and deg(k) = sum m_i, the quotient tends to
// C_a^deg(b) / C_b^deg(a) at s = oo.  C_k is read at one table point off the form's divisor and
// checked at the other table points coprime to the level.
//
// What this adds to tests/M0PoleSum.m (3d): there both quotients see the same gap of the -log m rule
// (standalone rem:divisor), log 3 - log(3/4) = 2 log 2.  Here the singular pairs at -4 are
// x = +-lambda_0 (m = 1), x = +-3 lambda_0 (m = 9, the pole of form 13 at q^-9) and x = +-lambda_0/2
// (m = 1/4), so the quotients through form 13 test the gap log 9 - log 1 = 2 log 3.
// Forms 12..15 have half-integral principal parts, outside Theorem B's hypothesis c_eta(-m) in Z;
// their doubles are inputs and the formula is linear in F, so everything is done with 2 f.
// Provenance of the expected values: the s-values are Guo-Yang's; the C_k are computed by this test
// from them and from the forms' divisors (not independent values); the quotients at -4 are then
// forced by the divisors alone.
// Points where the level prime 2 is RAMIFIED in the CM field (d_0 even: -84, -168, -232, ...) are left
// out: the local factor at a prime dividing gcd(d_0, N) is the known gap of SchoferFormula.m
// (CMNONCOPRIME); the conductor-2 points (-16, -28, -60, -100, -112, ...) stay, their terms are the
// fibre sum of prop:fibre.

procedure test_divisor_point_21_2()
    printf "  X0^21(2): quotients of forms at the pole tau_-4 of s against the divisors...";
    D := 21; N := 2;
    sv := AssociativeArray();
    for t in [<-7,-7>,<-15,-5/3>,<-16,1>,<-28,1/9>,<-60,9>,<-84,-3>,<-100,1/5>,<-112,25>,<-120,-1/3>,
              <-148,37/9>,<-168,0>,<-228,-25/3>,<-232,-32>,<-280,-35/9>,<-312,-8/3>,<-372,-3/4>,
              <-408,-75>,<-420,21>,<-532,-19/4>,<-708,25/48>,<-840,-16/3>] do
        sv[t[1]] := Rationals()!t[2];
    end for;
    Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
    Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};
    fs := BorcherdsForms(star, curves : Prec := 100);
    ks := Sort([k : k in Keys(fs)]);
    assert ks eq [-2, -1] cat [9..15];
    Ld := ShimuraCurveLattice(D, N);
    function cell_value(x)
        v := Rationals()!1;
        for q in Keys(x`log_coeffs) do
            c := x`log_coeffs[q];
            if not IsIntegral(c) then return false; end if;
            v *:= (Rationals()!q)^(Integers()!c);
        end for;
        return v;
    end function;
    pts := Sort([d : d in Keys(sv) | IsOdd(FundamentalDiscriminant(d))]);
    V := AssociativeArray();
    for d in pts do V[d] := SchoferFormula([fs[k] : k in ks], d, D, N, Ld); end for;
    C2 := AssociativeArray(); degs := AssociativeArray(); dvs := AssociativeArray();
    nchecked := 0;
    for i->k in ks do
        dv := DivisorOfBorcherdsForm(fs[k], star); dvs[k] := dv;
        for pr in dv do
            error if pr[1] ne -4 and not IsDefined(sv, pr[1]),
                Sprintf("form %o has a divisor point of discriminant %o not in the table", k, pr[1]);
        end for;
        model2 := func< d | &*[Rationals() | AbsoluteValue(sv[d] - sv[pr[1]])^(Integers()!(2*pr[2])) : pr in dv | pr[1] ne -4] >;
        off := [d : d in pts | not exists{pr : pr in dv | pr[1] eq d}];
        assert #off ge 2;
        cv := cell_value(2*V[off[1]][i]);
        error if cv cmpeq false, Sprintf("form %o: the double has an irrational value at d = %o", k, off[1]);
        C2[k] := cv / model2(off[1]);
        degs[k] := Integers()!&+[Rationals() | pr[2] : pr in dv | pr[1] ne -4];
        for d in off[2..#off] do
            cv := cell_value(2*V[d][i]);
            error if cv cmpeq false, Sprintf("form %o: the double has an irrational value at d = %o", k, d);
            error if cv ne C2[k]*model2(d),
                Sprintf("X0^21(2), form %o (divisor %o), d = %o: value %o of 2f, the table and the divisor give %o",
                        k, dv, d, cv, C2[k]*model2(d));
            nchecked +:= 1;
        end for;
    end for;
    // regression pins (computed by this test on 2026-10-04 from the table and the divisors): C_k^2
    assert C2[-2] eq 144 and C2[-1] eq 441 and C2[9] eq 63504 and C2[13] eq 21609/16777216;
    v4 := SchoferFormula([fs[k] : k in ks], -4, D, N, Ld);
    nq := 0;
    for i->a in ks, j->b in ks do
        if a ge b then continue; end if;
        na := 2*degs[b]; nb := 2*degs[a];
        got := cell_value(na*v4[i] - nb*v4[j]);
        want := C2[a]^degs[b] / C2[b]^degs[a];
        error if got cmpeq false or got ne want,
            Sprintf("X0^21(2) at d = -4: f_%o^%o / f_%o^%o is %o, the divisors give %o", a, na, b, nb, got, want);
        nq +:= 1;
    end for;
    assert nq eq 36;
    printf " ok (%o off-divisor values, %o quotients at the pole)\n", nchecked, nq;
end procedure;

printf "Testing a point on the divisor on X0^21(2)...\n";
test_divisor_point_21_2();
printf "Done!\n";
