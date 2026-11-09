// THE GENUS-FOUR MODEL OF X_0^2418(1)/<w_2,w_13,w_93>, CHECKED (the 2002 stepH-verify, transposed).
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
P<x> := PolynomialRing(Rationals());
C  := HyperellipticCurve(P![4,30,-6,-5,-9,-3,9], P![0,1,1]);              // the star (candidate model)
Cp := HyperellipticCurve(P![9,-27,31,-1,-5,-9,-3], P![0,1,1]);           // C' from the enumeration
Y  := HyperellipticCurve(P![169,0,172,0,5190,0,5068,0,3481,0,2304]);
Ym := ReducedMinimalWeierstrassModel(Y);
printf "Y : %o\n   genus %o, primes of the discriminant %o\n   reduced minimal model: %o\n", Y, Genus(Y),
       PrimeDivisors(Integers() ! Discriminant(Y)), Ym;
g  := 4*HyperellipticPolynomials(C) + P![0,1,1]^2;
e := -4/3; x0 := 4/3;
// u = (x - x0)/(x - e):  x = (e u - x0)/(u - 1);  q1(u) = (u-1)^6 g((e u - x0)/(u - 1))
q1 := P ! &+[ Coefficient(g, i) * (e*x - x0)^i * (x - 1)^(6-i) : i in [0..6] ];
printf "q1(u) = %o   (degree %o, q1(0) = %o)\n", q1, Degree(q1), Coefficient(q1, 0);
assert IsIsomorphic(HyperellipticCurve(q1), C);
printf "y^2 = q1(t^2) is Y over Q: %o\n", IsIsomorphic(HyperellipticCurve(Evaluate(q1, x^2)), Y);
F  := HyperellipticPolynomials(Y);
h1 := P ! [ Coefficient(F, 2*i) : i in [0..5] ];
printf "Y/(x -> -x) = {y^2 = h1(u)} isomorphic to the star over Q:  %o\n", IsIsomorphic(HyperellipticCurve(h1), C);
printf "Y/(x,y -> -x,-y) = {v^2 = u h1(u)} isomorphic to C' over Q: %o\n", IsIsomorphic(HyperellipticCurve(x*h1), Cp);
printf "fixed points of x -> -x on Y: %o\n", [ pt : pt in Points(Y : Bound := 1) | pt[1] eq 0 ];
H := {1, 2, 13, 26, 93, 186, 1209, 2418};
Yq := CreateShimuraQuot(2418, 1, H);  Yq`g := GenusShimuraCurveQuotient(2418, 1, H);
primes := [ p : p in PrimesUpTo(150) | 2418 mod p ne 0 ];
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
