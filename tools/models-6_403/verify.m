// THE GENUS-FOUR MODEL OF X_0^6(403)/<w_2,w_13,w_93>, CHECKED (given the candidate star model).
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
P<x> := PolynomialRing(Rationals());
C  := HyperellipticCurve(P![9,-21,33,-47,15,9,1], P![0,1,1]);              // the star (candidate model)
Cp := HyperellipticCurve(P![-351,-513,-187,75,25,-11,1], P![0,1,1]);        // C' from the enumeration
g  := 4*HyperellipticPolynomials(C) + P![0,1,1]^2;
e := 1; x0 := -7; c := 1;
q1 := P ! &+[ Coefficient(g, i) * (e*x - x0)^i * (x - 1)^(6-i) : i in [0..6] ];   // u = (x - x0)/(x - e)
printf "q1(u) = %o   (degree %o, q1(0) = %o)\n", q1, Degree(q1), Coefficient(q1, 0);
assert IsIsomorphic(HyperellipticCurve(q1), C);
Y0 := HyperellipticCurve(Evaluate(q1, x^2/c));
Ynice := HyperellipticCurve(P![1296,0,2745,0,9388,0,2630,0,316,0,9]);
printf "y^2 = 9x^10 + 316x^8 + 2630x^6 + 9388x^4 + 2745x^2 + 1296 is y^2 = q1(t^2) over Q: %o\n", IsIsomorphic(Ynice, Y0);
Y := ReducedMinimalWeierstrassModel(Y0);
printf "Y : %o\n   genus %o, primes of the discriminant %o\n", Y, Genus(Y), PrimeDivisors(Integers() ! Discriminant(Y));
F  := HyperellipticPolynomials(Y0);
h1 := P ! [ Coefficient(F, 2*i) : i in [0..5] ];
printf "Y/(x -> -x) = {y^2 = h1(u)} isomorphic to the star over Q:  %o\n", IsIsomorphic(HyperellipticCurve(h1), C);
printf "Y/(x,y -> -x,-y) = {v^2 = u h1(u)} isomorphic to C' over Q: %o\n", IsIsomorphic(HyperellipticCurve(x*h1), Cp);
H := {1, 2, 13, 26, 93, 186, 1209, 2418};
Yq := CreateShimuraQuot(6, 403, H);  Yq`g := GenusShimuraCurveQuotient(6, 403, H);
primes := [ p : p in PrimesUpTo(150) | 2418 mod p ne 0 ];
sing := []; bad := [];
for p in primes do
    tr := ComputePointsViaTrace(Yq, p, 1);
    try n := #Points(ChangeRing(Ynice, GF(p))); catch err n := -1; end try;
    if n eq -1 then Append(~sing, p); elif n ne tr then Append(~bad, <p, n, tr>); end if;
end for;
printf "%o good primes below 150: Y singular at %o, count mismatches %o\n", #primes, sing, bad;
Cps := IntegralModel(SimplifiedModel(Cp));
printf "C': conductor %o\n", Factorisation(Conductor(Cps));
lo, hi := RankBounds(Jacobian(Cps));
printf "    rank bounds %o %o, torsion of order %o\n", lo, hi, #TorsionSubgroup(Jacobian(Cps));
quit;
