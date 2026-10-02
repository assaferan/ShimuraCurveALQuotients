// Do the x != 0 terms of the m = 0 correction vanish on the 15_2 forms?
// They are nonzero only if some form has c_oo(-k^2|d|) != 0 at a firing discriminant d
// (d odd fundamental here, N = 2), with the cusp-0 side cancelling it in c_0(-k^2|d|).
AttachSpec("ShimuraQuotients.spec");
D := 15; N := 2; M := 60;
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
ks := Sort([k : k in Keys(fs)]);
exps := [3, 7, 15, 23, 31, 35, 39, 47, 55, 60, 28, 63, 12, 27, 75];
for k in ks do
    f := fs[k];
    foo := qExpansionAtoo(f, 80); f0 := qExpansionAt0(f, 80);
    printf "form %o: ord_oo = %o, ord_0 (in q^(1/%o)) = %o, f0 const = %o\n",
           k, Valuation(foo), M, Valuation(f0), Coefficient(f0, 0);
    for m in exps do
        coo := Coefficient(foo, -m);
        c0 := Coefficient(f0, -M*m);
        if coo ne 0 or c0 ne 0 then
            printf "   m = %o: c_oo(-m) = %o, cusp-0 part c_eta(-m) = %o, c_0(-m) = %o\n",
                   m, coo, c0, coo + c0;
        end if;
    end for;
end for;
