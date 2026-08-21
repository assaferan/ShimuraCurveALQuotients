// Size the vector-valued oracle run on a base, WITHOUT guessing the numerical budget.
//   magma DD:=14 NN:=5 vvdata/poles145.m
// Prints, per Borcherds form: the pole order at oo and the eta-quotient exponent spread, then the
// implied (Prec, NumSamples) from the two calibrated rules
//   * aliasing error ~ exp(4 pi sqrt(p K) - 2 pi K y),  y = 1   ->  choose K with that <= 1e-100
//   * Prec >= 2 pi p / ln 10 (dynamic range) + monomial cancellation + margin.
// The X0^15(2) calibration point: p = 30, exponents to +-21, K = 192 gave 8.7e-103, P = 200 converged.
AttachSpec("ShimuraQuotients.spec");

D := 14; N := 5;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;

t0 := Cputime();
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
printf "curves built (%o s)\n", Cputime(t0);

fs := BorcherdsForms(star, curves : Prec := 100);
ks := Sort([k : k in Keys(fs)]);
printf "BorcherdsForms done (%o s): %o forms, keys = %o\n", Cputime(t0), #ks, ks;

Ld := ShimuraCurveLattice(D, N);
M := IsOdd(D*N) select 4*D*N else 2*D*N;
printf "X0^%o(%o): M = %o, |disc grp| = %o, index Gamma_0(M) = %o\n",
       D, N, M, #Ld`disc_grp, Index(Gamma0(M));

printf "\n%-6o %-10o %-12o %-8o\n", "form", "pole ord", "max|exp|", "#terms";
pmax := 0;
emax := 0;
for k in ks do
    f := fs[k];
    q := qExpansionAtoo(f, 1 : RelPrec := true);
    p := -Valuation(q);
    exps := Exponents(f);
    e := IsEmpty(exps) select 0 else Maximum([Maximum([Abs(x) : x in r]) : r in exps]);
    printf "%-6o %-10o %-12o %-8o\n", k, p, e, #exps;
    pmax := Maximum(pmax, p);
    emax := Maximum(emax, e);
end for;

// --- the two rules, evaluated at the worst form -------------------------------------------------
RR := RealField(20);
p := RR!pmax;
y := RR!1;
Kneed := 0;
for K in [16*j : j in [1..64]] do
    err := 4*Pi(RR)*Sqrt(p*K) - 2*Pi(RR)*K*y;          // natural log of the aliasing error
    if err/Log(RR!10) le -100 then Kneed := K; break; end if;
end for;
// cancellation loss: on 15_2 (max exponent 21) it measured ~66 digits; scale linearly in the spread
cancel := Ceiling(66 * emax / 21);
Pneed := Ceiling(2*Pi(RR)*p/Log(RR!10)) + cancel + 50;

printf "\nworst form: pole order %o, exponent spread %o\n", pmax, emax;
printf "  aliasing rule  -> NumSamples K = %o", Kneed;
if Kneed gt 0 then
    printf "   (log10 err = %o)\n",
           ChangePrecision((4*Pi(RR)*Sqrt(p*Kneed) - 2*Pi(RR)*Kneed)/Log(RR!10), 6);
else
    printf "   *** NOT REACHED below K = 1024 ***\n";
end if;
printf "  precision rule -> Prec = %o  (dynamic range %o + cancellation %o + margin 50)\n",
       Pneed, Ceiling(2*Pi(RR)*p/Log(RR!10)), cancel;
printf "  15_2 reference: p = 30, spread 21, K = 192, P = 200, |G| = 1800, 19 min for 9 forms\n";
printf "  cost scales with |G| * K * (#cosets): |G| ratio = %o\n", RR!#Ld`disc_grp/1800;
printf "TOTAL %o s\n", Cputime(t0);
quit;
