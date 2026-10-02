AttachSpec("ShimuraQuotients.spec");
D := 21; N := 1;
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
for k in Sort([k : k in Keys(fs)]) do
    printf "form %o: divisor %o\n", k, DivisorOfBorcherdsForm(fs[k], Xstar);
end for;
