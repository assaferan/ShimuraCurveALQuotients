// The m = 0 multipliers at a COMPOSITE squarefree level, by support class (prop:composite of the
// standalone): M0MultipliersBySupport on X_0^6(35), M = 420, |L^v/L| = 88200, on WEIGHT-1/2 eta
// quotients of level 420 (sum of exponents 1, integral exponents at oo), the weight of a Borcherds input.  ⚠ compmult.m, compmult2*.m, compmult3.m used WEIGHT-0 quotients: the
// coset sum slashes in weight 1/2, so their "constant" carried a (c tau + d)^(-1/2) and the two-point
// check rightly failed wherever the q^0 coefficient was nonzero (and passed trivially where it was 0).
// What is checked: the routine runs at this size; the constant terms are constant on each support
// class {5}, {7}, {5,7}; whether the three classes carry DIFFERENT values, as prop:composite needs
// them to be able to.  No Borcherds form exists at any composite-level base yet.
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
SetVerbose("ShimuraQuotients", 2);
D := StringToInteger(D); N := StringToInteger(N);
M := IsOdd(D*N) select 4*D*N else 2*D*N;
Ld := ShimuraCurveLattice(D, N);
ds := Divisors(M);
// weight 1/2 (sum of exponents 1) with integral exponents at oo (sum d r_d = 0 mod 24), as a Borcherds
// input; f1, f2 have a simple pole at oo, f3 has order 0 there (found by a brute-force search).
A1 := AssociativeArray(); A1[1] := 2; A1[2] := 1; A1[14] := -2;                      // eta(t)^2 eta(2t) / eta(14t)^2
A2 := AssociativeArray(); A2[1] := 3; A2[2] := 3; A2[6] := -2; A2[7] := -3;          // eta(t)^3 eta(2t)^3 / (eta(6t)^2 eta(7t)^3)
A3 := AssociativeArray(); A3[5] := 2; A3[10] := -1;                                  // eta(5t)^2 / eta(10t)
ex := func< pairs | [Integers() | IsDefined(pairs, d) select pairs[d] else 0 : d in ds] >;
R := EtaQuotientsRing(M, D*N);
fs := [EtaQuotient(R, ex(A)) : A in [A1, A2, A3]];
Append(~fs, fs[1] + fs[2]);
printf "X_0^%o(%o), M = %o, |L^v/L| = %o\n", D, N, M, #Ld`disc_grp;
t := Cputime();
arrs := M0MultipliersBySupport(fs, Ld, D, N);
printf "\n%o s\n", Cputime(t);
for i->a in arrs do
    printf "form %o: %o\n", i, [<p, a[p]> : p in Sort(SetToSequence(Keys(a)))];
end for;
