// F = fs[-2]/fs[-1] at d = -12 on X_0^15(2): Table 45 forces 1/20 (F = (s-2)/(20 s), s(-12) = oo)
AttachSpec("ShimuraQuotients.spec");
SetColumns(0); SetVerbose("ShimuraQuotients", 1);
D := 15; N := 2;
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
Ld := ShimuraCurveLattice(D, N);
etas := [fs[-2], fs[-1], fs[10]];
vals := SchoferFormula(etas, -12, D, N, Ld);
printf "\nfs[-2](-12) = %o\nfs[-1](-12) = %o\nfs[10](-12) = %o\nF = fs[-2]/fs[-1] = %o   (truth 1/20 = -2Log2-Log5)\nF2 = fs[10]/fs[-1]^3 = %o   (truth 9/(2^16 5) = -16Log2+2Log3-Log5)\n", vals[1], vals[2], vals[3], vals[1] - vals[2], vals[3] - 3*vals[2];
