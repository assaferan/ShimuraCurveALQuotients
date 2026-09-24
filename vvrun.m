// Research driver for the vector-valued constant terms.
//   magma DD:=15 NN:=2 PR:=200 KK:=192 vvrun.m
// Prints, per Borcherds form of the base: the principal-part check error (the correctness gate),
// c_eta(0) at every isotropic coset, and the resulting m = 0 multiplier (1/2) c_eta(0).
AttachSpec("ShimuraQuotients.spec");

D := 15; N := 2;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
PREC := 200; if assigned PR then PREC := StringToInteger(PR); end if;
KS := 192;   if assigned KK then KS := StringToInteger(KK); end if;

t0 := Cputime();
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
ks := Sort([k : k in Keys(fs)]);
forms := [fs[k] : k in ks];
Ld := ShimuraCurveLattice(D, N);
M := IsOdd(D*N) select 4*D*N else 2*D*N;
printf "X0^%o(%o): %o forms, M = %o, |G| = %o, keys = %o\n", D, N, #ks, M, #Ld`disc_grp, ks;

consts, isoelts, errs := VVConstantTerms(forms, Ld, M : Prec := PREC, NumSamples := KS);
printf "#isotropic = %o (2N-1 = %o)\n\n", #isoelts, 2*N-1;
printf "%-6o %-14o %-22o %o\n", "form", "max|PPerr|", "c_0(0)", "c_eta(0) at nonzero isotropic";
for i->k in ks do
    printf "%-6o %-14o %-22o %o\n", k, ChangePrecision(errs[i], 6),
           ChangePrecision(consts[i][1], 10), [ChangePrecision(c, 10) : c in consts[i][2..#consts[i]]];
end for;
printf "\nmultipliers (1/2)c_eta(0): %o\n",
       [ChangePrecision(Real(c[2])/2, 10) : c in consts];
printf "TOTAL %o s\n", Cputime(t0);
quit;
