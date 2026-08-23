// The degree-weighting experiment on X0^10(23): raw Schofer LogSums for the
// hauptmodul rows at the two anchor discs and the failing quadratic disc -68.
// Offline analysis then tests which power of 23 in the -68 cell makes a
// sign-candidate minpoly admissible in the star field.
AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 1);
D := 10; N := 23;
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
Ld := ShimuraCurveLattice(D, N);
for d in [-40, -43, -68] do
    vals := AbsoluteValuesAtRationalCMPoint([fs[-1], fs[-2]], d, star, Ld);
    printf "RAW d=%o  s: %o\n", d, vals[1];
    printf "RAW d=%o  st: %o\n", d, vals[2];
end for;
quit;
