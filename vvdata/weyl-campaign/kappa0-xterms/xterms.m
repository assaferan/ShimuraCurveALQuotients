// For the forms of X_0^D(N): the pole order at the cusp 0 as a fraction of the cusp width M.
// If it is below the smallest Q(x) over x != 0 on the CM line at a firing discriminant (at least
// 3/4 in every case, and |d|/4 or |d| at a given d), F_f has no principal part at the cosets
// [x + nu] and the x != 0 terms of the m = 0 correction vanish.
AttachSpec("ShimuraQuotients.spec");
D := StringToInteger(D_s); N := StringToInteger(N_s);
M := IsOdd(D*N) select 4*D*N else 2*D*N;
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
ks := Sort([k : k in Keys(fs)]);
printf "X_0^%o(%o), M = %o\n", D, N, M;
for k in ks do
    f := fs[k];
    foo := qExpansionAtoo(f, 80); f0 := qExpansionAt0(f, 80);
    v0 := Valuation(f0);
    nz := [j : j in [1..-v0] | Coefficient(f0, -j) ne 0 and j ge 3*M/4];
    printf "form %o: ord_oo = %o, ord_0 = %o = %o * M; cusp-0 exponents >= 3/4 with nonzero coefficient: %o\n",
           k, Valuation(foo), v0, -v0/M, [j/M : j in nz];
end for;
printf "XTERMS_DONE\n";
