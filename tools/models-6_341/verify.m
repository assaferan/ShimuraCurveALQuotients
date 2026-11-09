// THE GENUS-FOUR MODEL OF X_0^6(341)/<w_3,w_11,w_62>, CHECKED (given the candidate star model).
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
P<x> := PolynomialRing(Rationals());
C  := HyperellipticCurve(P![9,-33,9,29,-12,6], P![0,1,1]);              // the star candidate
Cp := HyperellipticCurve(P![9,39,0,-2,-1,-2], P![0,1,1]);               // C' from the enumeration
g  := P![0,24,-47,118,37,-132,36];                                      // x -> 1/x of 4f + h^2
assert IsIsomorphic(HyperellipticCurve(g), C);
e := 3; x0 := 1/3; c := -3;
q1 := P ! &+[ Coefficient(g, i) * (e*x - x0)^i * (x - 1)^(6-i) : i in [0..6] ];   // u = (x - x0)/(x - e)
printf "q1(u) = %o   (degree %o, q1(0) = %o)\n", q1, Degree(q1), Coefficient(q1, 0);
assert IsIsomorphic(HyperellipticCurve(q1), C);
Y0 := HyperellipticCurve(Evaluate(q1, x^2/c));
Y  := ReducedMinimalWeierstrassModel(Y0);
printf "Y : %o\n   genus %o, primes of the discriminant %o\n", Y, Genus(Y), PrimeDivisors(Integers() ! Discriminant(Y));
F  := HyperellipticPolynomials(Y0);
h1 := P ! [ Coefficient(F, 2*i) : i in [0..5] ];
printf "Y/(x -> -x) = {y^2 = h1(u)} isomorphic to the star over Q:  %o\n", IsIsomorphic(HyperellipticCurve(h1), C);
printf "Y/(x,y -> -x,-y) = {v^2 = u h1(u)} isomorphic to C' over Q: %o\n", IsIsomorphic(HyperellipticCurve(x*h1), Cp);
H := {1, 3, 11, 33, 62, 186, 682, 2046};
Yq := CreateShimuraQuot(6, 341, H);  Yq`g := GenusShimuraCurveQuotient(6, 341, H);
primes := [ p : p in PrimesUpTo(150) | 2046 mod p ne 0 ];
sing := []; bad := [];
for p in primes do
    tr := ComputePointsViaTrace(Yq, p, 1);
    try n := #Points(ChangeRing(Y, GF(p))); catch err n := -1; end try;
    if n eq -1 then Append(~sing, p); elif n ne tr then Append(~bad, <p, n, tr>); end if;
end for;
printf "%o good primes below 150: Y singular at %o, count mismatches %o\n", #primes, sing, bad;
Cps := IntegralModel(SimplifiedModel(Cp));
printf "C': conductor %o\n", Factorisation(Conductor(Cps));
lo, hi := RankBounds(Jacobian(Cps));
printf "    rank bounds %o %o, torsion of order %o\n", lo, hi, #TorsionSubgroup(Jacobian(Cps));
quit;
