// Diagnostic: isolate what drives the error on the deep-pole forms of 15_2.
// One rho build, then sweep the number of sample points K, reporting WHERE the principal-part
// check fails (component, exponent, numeric vs predicted) instead of just the max.
AttachSpec("ShimuraQuotients.spec");
Attach("vvlib.m");

D := 15; N := 2;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
PREC := 120; if assigned PR then PREC := StringToInteger(PR); end if;
FORM := -2;  if assigned FM then FORM := StringToInteger(FM); end if;
KLIST := [48, 64, 96];

t0 := Cputime();
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
Rq := Parent(fs[FORM]); M := Rq`M; ds := Rq`ds;
f := fs[FORM];
foo := qExpansionAtoo(f, 100); f0 := qExpansionAt0(f, 100);
printf "form %o: v_oo = %o, v_0 = %o, M = %o\n", FORM, Valuation(foo), Valuation(f0), M;
printf "exponent vectors: %o\n", Exponents(f);
printf "coeffs: %o\n", [f`coeffs[r] : r in Exponents(f)];

Ld := ShimuraCurveLattice(D, N);
CC := ComplexField(PREC); ii := CC.1; twopii := 2*Pi(CC)*ii;
S, Td, elts, i0 := VVRho(Ld, CC : Dual := true);
n := #elts;
Q := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
vs := [ChangeRing(g@@Ld`to_disc, Rationals()) : g in elts];
nmv := [ (vs[i]*Q, vs[i])/(2*dn^2) : i in [1..n] ];
res := [ IsIntegral(M*nmv[i]) select (Integers()!(M*nmv[i]) mod M) else -1 : i in [1..n] ];
iso := [i : i in [1..n] | IsIntegral(nmv[i])];
reps := VVCosetReps(M); words := [VVSTWord(g) : g in reps];
t2 := Cputime();
U := ZeroMatrix(CC, #reps, n);
for k in [1..#reps] do
    v := VVRhoInvE0(S, Td, words[k], i0);
    for i in [1..n] do U[k][i] := v[i][1]; end for;
end for;
printf "rho vectors (%o s)\n", Cputime(t2);

for KS in KLIST do
    taus := [CC!(t/KS) + ii : t in [0..KS-1]];
    A := ZeroMatrix(CC, KS, #reps);
    mx := 0.0;
    for t in [1..KS] do
        for k in [1..#reps] do
            z := taus[t]; fac := CC!1; w := words[k];
            for r := #w to 1 by -1 do
                if w[r][1] eq "S" then fac /:= Sqrt(z); z := -1/z; else z := z + w[r][2]; end if;
            end for;
            A[t][k] := fac * VVEtaEval(f, z);
            mx := Maximum(mx, Abs(A[t][k]));
        end for;
    end for;
    Fv := A*U;
    maxF := Maximum([Abs(Fv[t][i]) : t in [1..KS], i in [1..n]]);
    coef := function(i, nn)
        s := CC!0;
        for t in [1..KS] do s +:= Fv[t][i] * Exp(-twopii*CC!nn*taus[t]); end for;
        return s/KS;
    end function;
    best := <0.0, 0, 0>;
    v0 := Minimum(0, Valuation(f0));
    for i in [1..n] do
        for j in [1..-v0] do
            if (j mod M) ne res[i] then continue; end if;
            nn := -(Rationals()!j)/M;
            pred := CC!Coefficient(f0, -j);
            if i eq i0 and IsIntegral(nn) then pred +:= Coefficient(foo, Integers()!nn); end if;
            e := Abs(coef(i,nn) - pred);
            if e gt best[1] then best := <e, i, nn>; end if;
        end for;
    end for;
    for nn in [Valuation(foo)..-1] do
        if (-M*nn) le -v0 then continue; end if;
        e := Abs(coef(i0,nn) - CC!Coefficient(foo, nn));
        if e gt best[1] then best := <e, i0, nn>; end if;
    end for;
    printf "K=%-4o max|f|gamma| = 10^%-6o max|F| = 10^%-6o | worst PP: i=%o nn=%o err=%o (num=%o pred=%o) | c_0(0)=%o sum_iso=%o\n",
        KS, Round(Log(10,mx)), Round(Log(10,maxF)), best[2], best[3], ChangePrecision(best[1],6),
        ChangePrecision(coef(best[2], best[3]), 8),
        (best[2] eq i0 and IsIntegral(best[3])) select Coefficient(foo, Integers()!best[3]) else 0,
        ChangePrecision(coef(i0,0),10), ChangePrecision(&+[CC|coef(i,0) : i in iso],10);
end for;
printf "TOTAL %o s\n", Cputime(t0);
quit;
