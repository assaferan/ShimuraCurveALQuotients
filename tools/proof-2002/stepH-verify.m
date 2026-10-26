// THE GENUS-FOUR MODEL Y_2, CHECKED  (paper, Section 3, parts (3) and (4)).
//
//   magma -b tools/proof-2002/stepH-verify.m < /dev/null
//
// Prints: the quintic model q1 of the star with its rational Weierstrass point at infinity and the
// fibre at infinity at u = 0; that y^2 = q1(t^2/2) is Y_2; that the two involutions of Y_2 have
// quotients isomorphic over Q to the star C and to C'; the fixed points of wbar on Y_2; the counts
// of Y_2 against the trace formula at every good prime below 150; and the conductor of C'.

AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
P<x> := PolynomialRing(Rationals());
C  := HyperellipticCurve(P![ 4, 64, 196, -456, 180, 4, 16 ], P![ 0, 1, 1 ]);
Cp := HyperellipticCurve(P![ 1, 9, 29, -44, 10, -2 ], P![ 0, 1, 1 ]);
Y2 := HyperellipticCurve(P![ 65536, 0, -16384, 0, 185856, 0, -279360, 0, 98289, 0, -1001 ]);
printf "Y_2 : %o\n   genus %o, primes of the discriminant %o\n", Y2, Genus(Y2),
       PrimeDivisors(Integers() ! Discriminant(Y2));

// the quintic model of the star:  u = 8/(8x+1),  q1(u) = 64 u^6 g(-1/8 + 1/u)
g  := 4*HyperellipticPolynomials(C) + P![0,1,1]^2;
q1 := P ! (64 * x^6 * Evaluate(g, -1/8 + 1/x));
printf "q1(u) = %o   (degree %o, q1(0) = %o)\n", q1, Degree(q1), Coefficient(q1, 0);
assert IsIsomorphic(HyperellipticCurve(q1), C);
b, m := IsIsomorphic(HyperellipticCurve(Evaluate(q1, x^2/2)), Y2);
printf "y^2 = q1(t^2/2) is Y_2 over Q: %o\n", b;

F  := HyperellipticPolynomials(Y2);
h1 := P ! [ Coefficient(F, 2*i) : i in [0..5] ];                 // F(x) = h1(x^2)
printf "Y_2/(x -> -x) = {y^2 = h1(u)} isomorphic to the star over Q:   %o\n", IsIsomorphic(HyperellipticCurve(h1), C);
printf "Y_2/(x,y -> -x,-y) = {v^2 = u h1(u)} isomorphic to C' over Q:  %o\n", IsIsomorphic(HyperellipticCurve(x*h1), Cp);
printf "fixed points of x -> -x on Y_2: %o\n", [ pt : pt in Points(Y2 : Bound := 1) | pt[1] eq 0 ];

H := {1, 2, 13, 26, 77, 154, 1001, 2002};
Yq := CreateShimuraQuot(2002, 1, H);  Yq`g := GenusShimuraCurveQuotient(2002, 1, H);
primes := [ p : p in PrimesUpTo(150) | 2002 mod p ne 0 ];
sing := []; bad := [];
for p in primes do
    tr := ComputePointsViaTrace(Yq, p, 1);
    try n := #Points(ChangeRing(Y2, GF(p))); catch err n := -1; end try;
    if n eq -1 then Append(~sing, p); elif n ne tr then Append(~bad, <p, n, tr>); end if;
end for;
printf "%o good primes below 150: Y_2 singular at %o, count mismatches %o\n", #primes, sing, bad;

Cps := IntegralModel(SimplifiedModel(Cp));
printf "C': conductor %o (Magma's Conductor; the exponent at 2 is Ogg's formula there)\n", Factorisation(Conductor(Cps));
lo, hi := RankBounds(Jacobian(Cps));
printf "    rank bounds %o %o, torsion of order %o\n", lo, hi, #TorsionSubgroup(Jacobian(Cps));
printf "    a_p for p <= 50: %o\n", [ <p, p + 1 - #Points(ChangeRing(Cp, GF(p)))> : p in primes | p le 50 ];
quit;
