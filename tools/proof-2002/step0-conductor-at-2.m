AttachSpec("/home/assaferan/GitHub/CHIMP/CHIMP.spec");
R<x> := PolynomialRing(Rationals());
C := HyperellipticCurve(16*x^6+4*x^5+180*x^4-456*x^3+196*x^2+64*x+4, x^2+x);
printf "odd conductor: %o = %o\n", OddConductor(C), Factorisation(OddConductor(C));
time e2 := EvenConductorExponent_Genus2(C, 2 : Proof := true);
printf "PROVED conductor exponent at 2: %o   => conductor = 2^%o * 7^2 * 11^2 * 13^2 = %o\n", e2, e2, 2^e2*1002001;
quit;
