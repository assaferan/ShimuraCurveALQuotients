// Dump the infinity-side principal parts of the two Hauptmodul forms, WITHOUT
// running ValuesAtCMPoints -- which is where 39_2 dies inside M0MultiplierExact.
AttachSpec("ShimuraQuotients.spec");
D := 39; N := 2;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
printf "BASE %o %o\n", D, N;
for k in [-1, -2] do
    e := qExpansionAtoo(fs[k], 1);
    v := Valuation(e);
    for n in [v..-1] do
        c := Coefficient(e, n);
        if c ne 0 then printf "PP %o %o %o\n", k, -n, c; end if;
    end for;
end for;
printf "DONE\n"; quit;
