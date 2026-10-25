AttachSpec("/home/assaferan/GitHub/CHIMP/CHIMP.spec");
R<x> := PolynomialRing(Rationals());
C := HyperellipticCurve(6*x^5 - 12*x^4 + 29*x^3 + 9*x^2 - 33*x + 9, x^2 + x);
printf "odd conductor: %o = %o\n", OddConductor(C), Factorisation(OddConductor(C));
time e2 := EvenConductorExponent_Genus2(C, 2 : Proof := true);
printf "PROVED conductor exponent at 2: %o   => conductor = 2^%o * 3^2 * 11^2 * 31^2\n", e2, e2;
quit;
