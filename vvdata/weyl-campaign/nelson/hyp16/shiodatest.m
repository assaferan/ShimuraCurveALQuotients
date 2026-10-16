// Can a genus-3 hyperelliptic curve be rebuilt from coordinate-free invariants, the way finish.m
// rebuilds a genus-2 one from Igusa-Clebsch?  Tested on the stored model of X_0^35(1), whose
// octic is known, so a failure here is Magma's reconstruction and not the numerics.
SetColumns(0);
P<x> := PolynomialRing(Rationals());
f0 := P![ -7/256, 0, -69/320, 0, -277/80000, 0, -1/40000, 0, -7/100000000 ];
C0 := HyperellipticCurve(f0);
printf "stored X_0^35(1): %o\n  genus %o\n", C0, Genus(C0);

SI, wt := ShiodaInvariants(f0);
printf "Shioda invariants (weights %o):\n", wt;
for i in [1..#SI] do printf "   J%o = %o\n", wt[i], SI[i]; end for;

// the absolute (coordinate-free AND scale-free) invariants are what a numerical computation can
// hope to recognise, since the fitted octic is only determined up to GL_2 and a scalar
DI := DiscriminantFromShiodaInvariants(SI);
printf "discriminant from invariants: %o\n", DI;

printf "\nnormalised (absolute) invariants:\n";
AI := ShiodaInvariants(f0 : normalize := true);
for i in [1..#AI] do printf "   j%o = %o\n", i, AI[i]; end for;

printf "\nreconstruction:\n";
try
    Cr := HyperellipticCurveFromShiodaInvariants(SI);
    printf "   from unnormalised: %o\n", Cr;
    printf "   isomorphic to stored over Q: %o\n", IsIsomorphic(C0, Cr);
catch e
    printf "   unnormalised reconstruction failed: %o\n", e`Object;
end try;
try
    Cr2 := HyperellipticCurveFromShiodaInvariants(AI);
    printf "   from normalised: %o\n", Cr2;
    b, m := IsIsomorphic(C0, Cr2);
    printf "   isomorphic to stored over Q: %o\n", b;
catch e
    printf "   normalised reconstruction failed: %o\n", e`Object;
end try;
printf "SHIODATEST DONE\n";
quit;
