// The m = 0 multipliers at a COMPOSITE squarefree level, by support class (prop:composite of the
// standalone): M0MultipliersBySupport on X_0^6(35), M = 420, |L^v/L| = 88200, on
// weight-0 eta quotients of level 420 with INTEGRAL pole orders at the cusps (24th powers of eta ratios; the first probe,
// compmult.m, used fractional orders and every constant term vanished), (no Borcherds form exists at any composite-level base yet).
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
A1 := AssociativeArray(); A1[1] := 24; A1[2] := -24;
A2 := AssociativeArray(); A2[5] := 24; A2[10] := -24;
A3 := AssociativeArray(); A3[7] := 24; A3[14] := -24;
R := EtaQuotientsRing(M, D*N);
f1 := EtaQuotient(R, ex(A1)); f2 := EtaQuotient(R, ex(A2)); f3 := EtaQuotient(R, ex(A3));
A4 := AssociativeArray(); A4[1] := 12; A4[35] := 12; A4[5] := -12; A4[7] := -12; f4 := EtaQuotient(R, ex(A4));
fs := [f1, f2, f3, f4];
printf "X_0^%o(%o), M = %o, |L^v/L| = %o\n", D, N, M, #Ld`disc_grp;
t := Cputime();
arrs := M0MultipliersBySupport(fs, Ld, D, N : Prec := 200);   // Prec 80 failed the two-point check at the 24th powers (word 422, W = 210, depth 5041)
printf "\n%o s\n", Cputime(t);
for i->a in arrs do
    printf "form %o: %o\n", i, [<p, a[p]> : p in Sort(SetToSequence(Keys(a)))];
end for;
