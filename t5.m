// Step 5: assemble F_f = sum_gamma (f|gamma) rho(gamma^{-1}) e_0 numerically, Fourier-extract, and
// VALIDATE the principal part against the pipeline's two scalar q-expansions (GY Lemma 24).
// If that matches, the n = 0 coefficient c_0(0) is the deliverable.
AttachSpec("ShimuraQuotients.spec");
Attach("vvlib.m");

D := 6; N := 1;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
PREC := 60;  if assigned PR then PREC := StringToInteger(PR); end if;
YY := 1;  if assigned YV then YY := StringToRational(YV); end if;
KS := 64;    if assigned KK then KS := StringToInteger(KK); end if;
WHICH := 1;  if assigned WH then WHICH := StringToInteger(WH); end if;

t0 := Cputime();
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
ks := Sort([k : k in Keys(fs)]);
printf "X0^%o(%o): %o forms in %o s; keys = %o\n", D, N, #ks, Cputime(t0), ks;
f := fs[ks[WHICH]];
R := Parent(f);  M := R`M;
foo := qExpansionAtoo(f, 80);  f0 := qExpansionAt0(f, 80);
printf "M = %o  disc = %o  val(oo) = %o  val(0) = %o\n", M, R`disc, Valuation(foo), Valuation(f0);

Ld := ShimuraCurveLattice(D, N);
CC := ComplexField(PREC);  ii := CC.1;

t1 := Cputime();
S, Td, elts, i0 := VVRho(Ld, CC : Dual := true);
n := #elts;
Sinv := S^(-1);
printf "rho built (|G| = %o) in %o s\n", n, Cputime(t1);

reps := VVCosetReps(M);
words := [VVSTWord(g) : g in reps];
assert &and[VVWordMatrix(words[k]) eq reps[k] : k in [1..#reps]];
printf "%o cosets\n", #reps;

t2 := Cputime();
U := ZeroMatrix(CC, #reps, n);          // row k = rho(gamma_k^{-1}) e_0
for k in [1..#reps] do
    v := VVRhoInvE0(S, Sinv, Td, words[k], i0);
    for i in [1..n] do U[k][i] := v[i][1]; end for;
end for;
printf "rho vectors in %o s\n", Cputime(t2);

// sample F on the horocycle Im = y
t3 := Cputime();
y := CC!YY;
Fvals := ZeroMatrix(CC, KS, n);
for t in [0..KS-1] do
    tau := CC!(t/KS) + y*ii;
    A := ZeroMatrix(CC, 1, #reps);
    for k in [1..#reps] do A[1][k] := VVSlashEval(f, words[k], tau); end for;
    row := A*U;
    for i in [1..n] do Fvals[t+1][i] := row[1][i]; end for;
end for;
printf "sampled %o points in %o s\n", KS, Cputime(t3);

// Fourier: component eta has exponents in n0_eta + Z with n0_eta = -Q(eta) mod 1
Q := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
vs := [ChangeRing(g@@Ld`to_disc, Rationals()) : g in elts];
nmv := [ (vs[i]*Q, vs[i])/(2*dn^2) : i in [1..n] ];
coef := function(i, nn)     // nn rational, nn = -nm_i mod 1
    s := CC!0;
    for t in [0..KS-1] do
        tau := CC!(t/KS) + y*ii;
        s +:= Fvals[t+1][i] * Exp(-2*Pi(CC)*ii*CC!nn*tau);
    end for;
    return s/KS;
end function;

// ---- report: eta = 0 component's principal part + constant term ----
printf "\n--- eta = 0 component (nm = %o) ---\n", nmv[i0];
printf "  n     numeric c_0(n)                      Coefficient(foo,n)\n";
for nn in [Valuation(foo)..4] do
    c := coef(i0, nn);
    printf "  %-5o %-38o %o\n", nn, c, Coefficient(foo, nn);
end for;

// ---- a few nonzero eta with the smallest positive nm, to see the 0-block ----
ord := Sort([<nmv[i] - Floor(nmv[i]), i> : i in [1..n] | i ne i0]);
seen := {};
printf "\n--- sample eta != 0 ---\n";
cnt := 0;
for z in ord do
    r := z[1]; i := z[2];
    if r in seen then continue; end if;
    Include(~seen, r);
    cnt +:= 1; if cnt gt 4 then break; end if;
    jj := Integers()!(M*r);            // the 0-block index: nm(eta) = j/M mod 1
    printf "  eta idx %-5o nm = %-10o  (j = %o)\n", i, r, jj;
    for nn in [-r, 1-r, 2-r] do
        printf "      c_eta(%-8o) = %-38o   Coefficient(f0, %o) = %o\n",
               nn, coef(i, nn), Integers()!(M*nn), Coefficient(f0, Integers()!(M*nn));
    end for;
end for;
printf "\nTOTAL %o s\n", Cputime(t0);
quit;
