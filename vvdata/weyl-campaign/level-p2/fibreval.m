// Values of the forms of X_0^15(2) at the seven Table-45 points with the m = 0-type terms for
// NONZERO cosets removed: FIBRE_NOOUTER (no outer term), FIBRE_SKIP (no Yang conductor term, no
// mu != 0 subtraction).  Prints the pairs (gamma, m, mu, x) with Q(x) = m met by the m > 0 loop.
AttachSpec("ShimuraQuotients.spec");
D := 15; N := 2;
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
svals := AssociativeArray();
svals[-7] := 1/4; svals[-15] := 5/4; svals[-52] := 1;
svals[-28] := 9/4; svals[-60] := -1/12; svals[-240] := -25/12; svals[-48] := -1/4;
want := Set(Keys(svals));
cm := CandidateDiscriminants(star, curves : Keep := want);
rat := cm[1]; quad := cm[2];
for t in [<-28, 2, 1>, <-240, 4, 1>, <-48, 4, 1>] do
    if not exists{u : u in rat | u[1] eq t[1]} then Append(~rat, t); end if;
end for;
tab, _ := AbsoluteValuesAtCMPoints(star, curves, [rat, quad], fs : MaxNum := 60, Prec := 100, Exclude := {}, Include := want);
ks := Sort([k : k in Keys(fs)]);
printf "form keys %o\n", ks;
C := 1280/9;
for d in [-7, -15, -52, -28, -60, -240, -48] do
    i := Index(tab`Discs, d);
    printf "VALUE d %o  fs[-2]: %o   truth: %o\n", d, tab`Values[1][i], C * AbsoluteValue(svals[d] * (svals[d] - 2));
    for j in [2..#ks] do printf "VALUE d %o  fs[%o]: %o\n", d, ks[j], tab`Values[j][i]; end for;
end for;
