// Generic panel dump: BASE + MONO lines (enum32k input); DD/NN required.
AttachSpec("ShimuraQuotients.spec");
D := StringToInteger(DD); N := StringToInteger(NN);
keystr := "-2,-1,9,10,11,12,13,14,15";
if assigned KEYS then keystr := KEYS; end if;
keys0 := [StringToInteger(s) : s in Split(keystr, ",")];
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
avail := Sort([k : k in Keys(fs)]);
printf "AVAIL = %o\n", avail;
keys := [k : k in keys0 | k in avail];
printf "USEDKEYS %o\n", &cat[Sprintf("%o,", k) : k in keys];
error if #keys eq 0, "no default keys available";
R := Parent(fs[keys[1]]); ds := R`ds;
M := IsOdd(D*N) select 4*D*N else 2*D*N;
printf "BASE %o %o  M = %o  ds = %o\n", D, N, M, ds;
monos := {@ @};
for k in keys do for r in Exponents(fs[k]) do Include(~monos, r); end for; end for;
printf "#monomials = %o\n", #monos;
for ri->r in monos do
    printf "MONO %o %o\n", ri, [r[i] : i in [1..#ds]];
end for;
quit;
