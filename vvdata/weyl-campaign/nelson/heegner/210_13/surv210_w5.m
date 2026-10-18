AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
D := 210; N := 13; W := {Integers() | 1,2,21,42,65,130,1365,2730};
X := CreateShimuraQuot(D, N, W); X`g := 3;
E := EllipticCurve("2730u1");
star := func<x, y | x*y div GCD(x, y)^2>;
Wp := W join {star(5, w) : w in W};
ok := false;
for m in Wp diff W do
    fp := NumFixedPointsByCMOrder(D, N, m);
    for d in Keys(fp) do
        if fp[d] eq 0 then continue; end if;
        F<a> := AbsoluteField(FieldOfDefinitionOfCMPoint(X, d));
        if IsCoercible(E(F), [9/2, 1/832*(-105*a - 2183), 1]) then P1 := E(F)![9/2, 1/832*(-105*a - 2183), 1]; P2 := E(F)![-63, 1/13*(-30*a + 433), 1]; ok := true; printf "field from disc %o: %o\n", d, DefiningPolynomial(F); break m; end if;
    end for;
end for;
assert ok;
B := Divisor(Place(P1)) + Divisor(Place(P2));
covers := DoubleCoversWithBranchDivisor(E, B);
printf "covers: noether %o\n", [c`noether : c in covers];
primes := [p : p in PrimesInInterval(3, 800) | (D*N) mod p ne 0];
target := [ComputePointsViaTrace(X, p, 1) : p in primes];
hits := MatchCoverToCounts(E, covers, target, primes : ExtraPrimes := PrimeDivisors(D*N));
printf "consistent (cover, kappa): %o\n", hits;
for h in hits do
    c := TwistedCover(covers[h[1]], h[2]);
    if c`noether ne 5 then printf "cover %o: not hyperelliptic\n", h[1]; continue; end if;
    H := HyperellipticModelOfCover(c);
    PrintFile("/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/heeg/octic210_w5.m", Sprintf("Px<T> := PolynomialRing(Rationals()); F8 := %o;", HyperellipticPolynomials(H)) : Overwrite := true);
    printf "raw model saved, %o digits\n", Max([#IntegerToString(Numerator(x)) : x in Coefficients(HyperellipticPolynomials(H))]);
    A, mA := AutomorphismGroup(H); printf "automorphism group order %o\n", #A;
end for;
printf "DONE\n";
