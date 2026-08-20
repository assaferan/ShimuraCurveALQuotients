// Test Yifan Yang, "Special values of hypergeometric functions...", Theorem B:
// the Borcherds form psi_F has weight c_0(0), the constant term of the eta = 0 component.
// If the m=0 multiplier the pipeline computes is really c_0(0), it should just BE this number.
AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 0);
D := StringToInteger(DD); N := StringToInteger(NN);
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
for k in Sort([k : k in Keys(fs)]) do
    foo := qExpansionAtoo(fs[k], 1);
    f0  := qExpansionAt0(fs[k], 1);
    printf "[C00] %o %o form=%-4o c_oo(0)=%-8o c_0(0)=%-8o\n",
           D, N, k, Coefficient(foo, 0), Coefficient(f0, 0);
end for;
exit;
