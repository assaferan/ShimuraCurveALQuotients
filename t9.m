// Structure of the constant term across the ISOTROPIC cosets: is c_eta(0) the same for every
// nonzero isotropic eta?  If so the multiplier is just (1/2) c_eta(0) and the 4(N-1) disappears.
AttachSpec("ShimuraQuotients.spec");
Attach("vvlib.m");

D := 6; N := 5;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
PREC := 120; if assigned PR then PREC := StringToInteger(PR); end if;
KS := 96;    if assigned KK then KS := StringToInteger(KK); end if;
FORM := 10;  if assigned FM then FORM := StringToInteger(FM); end if;

Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
f := fs[FORM]; Rq := Parent(f); M := Rq`M; ds := Rq`ds;
foo := qExpansionAtoo(f,100); f0 := qExpansionAt0(f,100);
printf "X0^%o(%o) form %o: v_oo = %o, v_0 = %o\n", D, N, FORM, Valuation(foo), Valuation(f0);

Ld := ShimuraCurveLattice(D, N);
CC := ComplexField(PREC); ii := CC.1; twopii := 2*Pi(CC)*ii;
S, Td, elts, i0 := VVRho(Ld, CC : Dual := true);
n := #elts;
Q := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
vs := [ChangeRing(g@@Ld`to_disc, Rationals()) : g in elts];
nmv := [ (vs[i]*Q, vs[i])/(2*dn^2) : i in [1..n] ];
iso := [i : i in [1..n] | IsIntegral(nmv[i])];
reps := VVCosetReps(M); words := [VVSTWord(g) : g in reps];
U := ZeroMatrix(CC, #reps, n);
for k in [1..#reps] do
    v := VVRhoInvE0(S, Td, words[k], i0);
    for i in [1..n] do U[k][i] := v[i][1]; end for;
end for;

taus := [CC!(t/KS) + ii : t in [0..KS-1]];
A := ZeroMatrix(CC, KS, #reps);
for t in [1..KS] do
    for k in [1..#reps] do
        z := taus[t]; fac := CC!1; w := words[k];
        for r := #w to 1 by -1 do
            if w[r][1] eq "S" then fac /:= Sqrt(z); z := -1/z; else z := z + w[r][2]; end if;
        end for;
        A[t][k] := fac * VVEtaEval(f, z);
    end for;
end for;
Fv := A*U;
coef := function(i, nn)
    s := CC!0;
    for t in [1..KS] do s +:= Fv[t][i] * Exp(-twopii*CC!nn*taus[t]); end for;
    return s/KS;
end function;

printf "\n#isotropic = %o (2N-1 = %o)\n", #iso, 2*N-1;
printf "  idx   eta (scaled dual coords)      c_eta(0)\n";
tot := CC!0;
for i in iso do
    c := coef(i, 0);
    tot +:= c;
    printf "  %-5o %-30o %o\n", i, Eltseq(vs[i]), ChangePrecision(c, 12);
end for;
printf "  SUM = %o\n", ChangePrecision(tot, 12);
quit;
