// ODD-N m=0 ORACLE, one base. Measures the true m=0 multiplier of every Borcherds form of X_0(D,N)*
// by the kernel-consistency method confirmed on X0^15(2) (266 independent checks) and X0^10(11)
// (reproduces Guo-Yang Table A.2).
//
// Requirements for the measurement to be IDENTIFIABLE: the CM point set must contain both discs where
// the m=0 term fires (p | N/(N,d_fund), i.e. p does not divide d_fund) and discs where it does not.
// If every disc fires, the missing term is a uniform row scaling that y^2 = f(s) absorbs and nothing
// can be measured -- which is exactly why this term went unnoticed. This worktree is configured with
// coprime_to_level := false so non-coprime (non-firing) discs enter the pool, and with the m=0 term
// OFF for odd N so the fitted delta is the FULL multiplier.
//
// The build is EXPECTED to fail at QuadraticConstraintsOnEquations -- that failure IS the signal, and
// the [RCE] dump printed before it carries the data. Feed the dump to fit_signs.py with p = N.
//
// usage: magma -b DD:=<D> NN:=<N> oracle_one.m

AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 0);

D := StringToInteger(DD); N := StringToInteger(NN);
printf "\n[ORACLE] ===== X0^%o(%o) =====\n", D, N;
t0 := Cputime();
try
    Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
    Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W);
    Xstar`CurveID := 0;
    curves := GetQuotientsAndGenera([Xstar]);
    assert exists(star){c : c in curves | IsStarCurve(c)};

    // Dump each form's PRINCIPAL PARTS before building. The measured delta is only half the data: a
    // closed formula for the multiplier must be a function of these coefficients, so without them the
    // sweep yields numbers with nothing to fit against. BorcherdsForms is cached, so computing it here
    // costs nothing extra -- AllEquationsAboveCovers below reuses it.
    fs := BorcherdsForms(star, curves : Prec := 100);
    for k in Sort([kk : kk in Keys(fs)]) do
        foo := qExpansionAtoo(fs[k], 1);
        f0 := qExpansionAt0(fs[k], 1);
        poo := [<m, Coefficient(foo, -m)> : m in [1..-Valuation(foo)] | Coefficient(foo, -m) ne 0];
        p0 := [<j, Coefficient(f0, -j)> : j in [1..-Valuation(f0)] | Coefficient(f0, -j) ne 0];
        printf "[PP] %o %o form=%o oo=%o zero=%o\n", D, N, k, poo, p0;
    end for;

    eqns, ws := AllEquationsAboveCovers(star, curves : Prec := 100);
    printf "[ORACLE] %o %o BUILD_OK (term already consistent; multipliers may all be 0)\n", D, N;
catch e
    msg := Join(Split(Sprintf("%o", e`Object), "\n"), " | ");
    if #msg gt 200 then msg := msg[1..200]; end if;
    printf "[ORACLE] %o %o BUILD_FAIL %o\n", D, N, msg;
end try;
printf "[ORACLE] %o %o elapsed %o s\n", D, N, Cputime(t0);
exit;
