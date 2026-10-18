// The b = -5 survivor of the 1326_1 sweep: pin its twist with many primes and write its model over the conic.
AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
load "/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/heeg/overconic.m";
D := 1326; N := 1; W := {Integers() | 1,17,78,1326};
X := CreateShimuraQuot(D, N, W); X`g := 3;
F<a> := AbsoluteField(FieldOfDefinitionOfCMPoint(X, -408));
printf "F = %o\n", DefiningPolynomial(F);
E := EllipticCurve("1326b2");
P1 := E(F)![1/336973*(-3819*a + 4276685), 1/54252653*(-7018932*a + 697001641), 1];
P2 := E(F)![1/217156*(3819*a + 3550145), 1/50597348*(10429599*a + 1232055397), 1];
printf "heights %o %o\n", Height(P1), Height(P2);
B := Divisor(Place(P1)) + Divisor(Place(P2));
covers := DoubleCoversWithBranchDivisor(E, B);
printf "covers: noether %o\n", [c`noether : c in covers];
primes := [p : p in PrimesInInterval(3, 700) | (D*N) mod p ne 0];
time target := [ComputePointsViaTrace(X, p, 1) : p in primes];
SetVerbose("ExactDoubleCover", 1);
time hits := MatchCoverToCounts(E, covers, target, primes : ExtraPrimes := PrimeDivisors(D*N));
printf "consistent (cover, kappa): %o\n", hits;
for h in hits do
    c := TwistedCover(covers[h[1]], h[2]);
    printf "cover %o twisted by %o: genus %o, noether %o\n", h[1], h[2], c`genus, c`noether;
    if c`noether eq 5 then
        time H := HyperellipticModelOfCover(c);
        printf "model degree %o, coefficient size %o digits\n", Degree(HyperellipticPolynomials(H)), Max([#IntegerToString(Numerator(a)) : a in Coefficients(HyperellipticPolynomials(H))]);
        try
            time Hr := ReducedModel(H);
            printf "REDUCED MODEL: y^2 = %o\n", HyperellipticPolynomials(Hr);
        catch err
            printf "ReducedModel failed: %o\n", err`Object; Hr := H;
        end try;
        PrintFile("/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/heeg/octic1326.m", Sprintf("Px<T> := PolynomialRing(Rationals()); F8 := %o;", HyperellipticPolynomials(Hr)) : Overwrite := true);
        // normalise with the Klein four group: find the involutions as Mobius maps
        A, mA := AutomorphismGroup(Hr);
        invs := [mA(a) : a in A | Order(a) eq 2];
        printf "automorphism group order %o, %o involutions\n", #A, #invs;
        Pt<t> := PolynomialRing(Rationals());
        mats := [];
        for phi in invs do
            eqs := DefiningPolynomials(phi);   // images of (x, y, z) on the weighted model
            // x-coordinate map: eqs[1]/eqs[3] as a Mobius transformation in x = X/Z
            num := eqs[1]; den := eqs[3];
            R := Parent(num);
            a := MonomialCoefficient(num, R.1); b := MonomialCoefficient(num, R.3);
            cc := MonomialCoefficient(den, R.1); dd := MonomialCoefficient(den, R.3);
            if [a,b,cc,dd] eq [1,0,0,1] then continue; end if;   // the hyperelliptic involution
            Append(~mats, Matrix(Rationals(), 2, 2, [a,b,cc,dd]));
        end for;
        printf "reduced involutions (Mobius matrices): %o\n", mats;
        PrintFile("/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/heeg/octic1326.m", Sprintf("mats := %o;", mats));
    end if;
end for;
printf "DONE\n";
