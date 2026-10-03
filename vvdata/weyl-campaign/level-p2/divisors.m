AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
D := 15; N := 2;
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
for k in Sort([k : k in Keys(fs)]) do
    printf "DIV form %o: %o\n", k, DivisorOfBorcherdsForm(fs[k], star);
    if k gt 0 then printf "    curve %o: W = %o genus %o\n", k, curves[k]`W, curves[k]`g; end if;
    foo := qExpansionAtoo(fs[k], 1); f0 := qExpansionAt0(fs[k], 1);
    printf "    principal part at oo: %o\n    principal part at 0 (q^(1/60)): %o\n",
        [<m, Coefficient(foo, m)> : m in [Valuation(foo)..-1] | Coefficient(foo, m) ne 0],
        [<m, Coefficient(f0, m)> : m in [Valuation(f0)..-1] | Coefficient(f0, m) ne 0];
end for;
