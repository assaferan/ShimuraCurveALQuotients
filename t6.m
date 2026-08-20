// Step 6: for EVERY Borcherds form of a base, build F_f numerically and report
//   (A) max deviation of the numeric principal part from the code's prediction
//         c_eta(n) = [eta=0][n in Z] Coefficient(foo,n)  +  [M*nm(eta) = -M*n mod M] Coefficient(f0, M*n)
//       -- i.e. GY Lemma 24 in the pipeline's own normalisation.  This is the correctness check.
//   (B) the DELIVERABLE: c_0(0), and sum over isotropic eta of c_eta(0).
AttachSpec("ShimuraQuotients.spec");
Attach("vvlib.m");

D := 6; N := 1;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
PREC := 60;  if assigned PR then PREC := StringToInteger(PR); end if;
KS := 64;    if assigned KK then KS := StringToInteger(KK); end if;
YY := 1;     if assigned YV then YY := StringToInteger(YV); end if;

t0 := Cputime();
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
ks := Sort([k : k in Keys(fs)]);
M := Parent(fs[ks[1]])`M;
printf "X0^%o(%o): %o forms, M = %o, keys = %o  (%o s)\n", D, N, #ks, M, ks, Cputime(t0);

Ld := ShimuraCurveLattice(D, N);
CC := ComplexField(PREC);  ii := CC.1;
t1 := Cputime();
S, Td, elts, i0 := VVRho(Ld, CC : Dual := true);
n := #elts;
// spot-check unitarity on a few rows instead of forming S*S^-1 (O(n^3) is hopeless at n = 1800)
for i in [1, 2, n] do
    r := &+[S[i][k]*ComplexConjugate(S[i][k]) : k in [1..n]];
    assert Abs(r - 1) lt 1e-20;
end for;
printf "rho: |G| = %o  (%o s)\n", n, Cputime(t1);

Q := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
vs := [ChangeRing(g@@Ld`to_disc, Rationals()) : g in elts];
nmv := [ (vs[i]*Q, vs[i])/(2*dn^2) : i in [1..n] ];
res := [ Integers()!(M*nmv[i]) mod M : i in [1..n] ];       // the bucket index j of eta
iso := [i : i in [1..n] | IsIntegral(nmv[i])];

reps := VVCosetReps(M);
words := [VVSTWord(g) : g in reps];
t2 := Cputime();
U := ZeroMatrix(CC, #reps, n);
for k in [1..#reps] do
    v := VVRhoInvE0(S, Td, words[k], i0);
    for i in [1..n] do U[k][i] := v[i][1]; end for;
end for;
printf "%o coset rho-vectors (%o s)\n", #reps, Cputime(t2);

y := CC!YY;
taus := [CC!(t/KS) + y*ii : t in [0..KS-1]];

printf "\n%-6o %-8o %-8o %-14o %-26o %o\n", "form", "val_oo", "val_0", "max|PP err|", "c_0(0)", "sum_iso c_eta(0)";
for kk in ks do
    f := fs[kk];
    foo := qExpansionAtoo(f, 80);  f0 := qExpansionAt0(f, 80);
    tt := Cputime();
    Fv := ZeroMatrix(CC, KS, n);
    for t in [1..KS] do
        A := ZeroMatrix(CC, 1, #reps);
        for k in [1..#reps] do A[1][k] := VVSlashEval(f, words[k], taus[t]); end for;
        row := A*U;
        for i in [1..n] do Fv[t][i] := row[1][i]; end for;
    end for;
    coef := function(i, nn)
        s := CC!0;
        for t in [1..KS] do s +:= Fv[t][i] * Exp(-2*Pi(CC)*ii*CC!nn*taus[t]); end for;
        return s/KS;
    end function;
    // (A) principal part check
    err := 0.0;
    lo := Minimum(Valuation(foo), Valuation(f0) - 1);
    for i in [1..n] do
        for j in [1..-Minimum(0,Valuation(f0))] do
            if (j mod M) ne res[i] then continue; end if;
            nn := -(Rationals()!j)/M;
            pred := CC!Coefficient(f0, -j);
            if i eq i0 and IsIntegral(nn) then pred +:= Coefficient(foo, Integers()!nn); end if;
            err := Maximum(err, Abs(coef(i,nn) - pred));
        end for;
        if i eq i0 then    // integer n not covered above (|n| beyond the 0-block range)
            for nn in [Valuation(foo)..-1] do
                if (-M*nn) le -Minimum(0,Valuation(f0)) then continue; end if;
                pred := CC!Coefficient(foo, nn);
                err := Maximum(err, Abs(coef(i,nn) - pred));
            end for;
        end if;
    end for;
    c00 := coef(i0, 0);
    ciso := &+[CC | coef(i, 0) : i in iso];
    printf "%-6o %-8o %-8o %-14o %-26o %o   (%o s)\n", kk, Valuation(foo), Valuation(f0),
           err, c00, ciso, Cputime(tt);
end for;
printf "\nTOTAL %o s\n", Cputime(t0);
quit;
