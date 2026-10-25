// Can the second candidate C_A be the star?  The certificate puts the branch pair B of Y -> star in
// {K, K + t_A}.  B ~ K needs a rational theta characteristic (none: the Weierstrass points form Galois
// orbits of sizes 2 and 4, so no odd-size subset is stable); B ~ K + t_A means B = W1 + W2, the
// conjugate pair, and a rational cover branched there needs a rational L with L^2 = O(W1 + W2),
// i.e. [W1 + W2 - 2 P0] in 2 J_A(Q) for a rational point P0.
SetColumns(0);
P<x> := PolynomialRing(Rationals());
CA := HyperellipticCurve(P![144,-120,265,-82,133,-12,36]);
JA := Jacobian(CA);
f := HyperellipticPolynomials(CA);
printf "Weierstrass factorisation: %o\n", Factorisation(f);
printf "Galois group of the quartic factor: %o\n", GaloisGroup([g[1] : g in Factorisation(f) | Degree(g[1]) eq 4][1]);
pts := Points(CA : Bound := 5000); printf "rational points: %o\n", pts;
T, mT := TorsionSubgroup(JA); tA := mT(T.1); printf "torsion point t_A = %o\n", tA;
q := [g[1] : g in Factorisation(f) | Degree(g[1]) eq 2][1];
W12 := JA![q, 0];                       // [W1 + W2 - (oo+ + oo-)]
printf "W1 + W2 - (oo+ + oo-) = %o, equals t_A: %o, 2x = 0: %o\n", W12, W12 eq tA, 2*W12 eq JA!0;
P0 := pts[1]; printf "base point P0 = %o\n", P0;
// [P - P0] for the rational points; [oo- - oo+]
inf := [pt : pt in pts | pt[3] eq 0];
dinf := JA![inf[1], inf[2]] ;            // oo_a - oo_b? Magma: J![P, Q] = [P - Q]
printf "[oo1 - oo2] = %o, order %o\n", dinf, Order(dinf);
// the class [W1 + W2 - 2 P0] = W12 + [(oo+ + oo-) - 2 P0] = W12 + [oo+ - P0] + [oo- - P0]
c := W12 + (JA![inf[1], P0]) + (JA![inf[2], P0]);
printf "class [W1 + W2 - 2 P0] = %o, order %o\n", c, Order(c);
printf "divisible by 2 in J_A(Q): %o\n", IsDivisibleBy(c, 2);
printf "rank bounds of J_A(Q): %o %o\n", RankBounds(JA);
// for comparison, the same on C: [W1 + W2 - 2 P0] with W1, W2 the rational Weierstrass points
C := HyperellipticCurve(P![9,-33,9,29,-12,6], P![0,1,1]); Cs := SimplifiedModel(C); J := Jacobian(Cs);
cp := Points(Cs : Bound := 100);
Wc := [pt : pt in cp | pt[2] eq 0]; Q0 := [pt : pt in cp | pt[2] ne 0][1];
cc := (J![Wc[1], Q0]) + (J![Wc[2], Q0]);
printf "on C: [W1 + W2 - 2 P0] = %o, order %o, divisible by 2: %o\n", cc, Order(cc), IsDivisibleBy(cc, 2);
quit;
