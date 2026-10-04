// The m = 0 multipliers at a COMPOSITE squarefree level, by support class (prop:composite of the
// standalone): M0MultipliersBySupport on X_0^6(35), M = 420, |L^v/L| = 88200, on arbitrary
// weight-0 eta quotients of level 420 (no Borcherds form exists at any composite-level base yet).
// What is checked: the routine runs at this size; the constant terms are constant on each support
// class {5}, {7}, {5,7} (lem:iso / the rho^S restriction law of ctdgap.m); the three classes carry
// DIFFERENT values, so a single multiplier -- the prime-level rule -- would be wrong here.
// Run from a code tree carrying M0MultipliersBySupport (branch kappa0-proof).
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
SetVerbose("ShimuraQuotients", 2);
D := StringToInteger(D); N := StringToInteger(N);
M := IsOdd(D*N) select 4*D*N else 2*D*N;
Ld := ShimuraCurveLattice(D, N);
ds := Divisors(M);
ex := func< pairs | [Integers() | IsDefined(pairs, d) select pairs[d] else 0 : d in ds] >;
A1 := AssociativeArray(); A1[1] := 1; A1[2] := -1;
A2 := AssociativeArray(); A2[5] := 2; A2[7] := -1; A2[35] := -1;
A3 := AssociativeArray(); A3[3] := 1; A3[4] := 1; A3[12] := -2;
R := EtaQuotientsRing(M, D*N);
f1 := EtaQuotient(R, ex(A1)); f2 := EtaQuotient(R, ex(A2)); f3 := EtaQuotient(R, ex(A3));
fs := [f1, f2, f3, f1 + 2*f2 - f3];
printf "X_0^%o(%o), M = %o, |L^v/L| = %o\n", D, N, M, #Ld`disc_grp;
t := Cputime();
arrs := M0MultipliersBySupport(fs, Ld, D, N);
printf "\n%o s\n", Cputime(t);
for i->a in arrs do
    printf "form %o: %o\n", i, [<p, a[p]> : p in Sort(SetToSequence(Keys(a)))];
end for;
