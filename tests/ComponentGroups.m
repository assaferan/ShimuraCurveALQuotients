// Component groups at p | D of the Jacobian of a star, the reduction of CM points to them, and the refinement of the
// torsion bound (ComponentGroups.m), against values this repository does not produce.  About six minutes.
//
// [1] OUTSIDE CHECK. Two stars of genus one whose curve is known: X_0^15(38)* = 570e4 and X_0^6(323)* = 1938h1 (the
//     traces of Frobenius to the Sturm bound and the component groups pick the curve within its class; the labels
//     are this repository's result, but the Tamagawa numbers below are Cremona's for those curves).  The component
//     group of the Jacobian at a prime of multiplicative reduction is cyclic of order the Tamagawa number, so the
//     critical group of the quotient graph must be Z/8 and Z/4 for 570e4 at 3 and 5, and Z/14 and Z/3 for 1938h1
//     at 2 and 3.  The rational 2-torsion point of 1938h1 lies off the identity component at 2 and on it at 3
//     (Cremona's tables, or LocalInformation), which is what the injection criterion must reflect: the 2-part can
//     inject into Phi_3 only through the torus, and the criterion refuses it.
// [2] A DIVISOR OF KNOWN CLASS. On X_0^15(38)* the pairing P_-228 + P_-760 - P_-12 - P_-40 is principal
//     (paper/hyperelliptic-15-38.tex), so its class is zero at 3 and at 5, while the other pairing
//     P_-228 + P_-12 - P_-760 - P_-40 is not even torsion and has a nonzero class at 3.
// [3] THE REFINED BOUND. On X_0^6(649)* the rational torsion has order dividing 6, the three fibre differences
//     have class zero in Phi_3 = Z/102, and the 3-part injects at 3, so the bound on the order of a torsion fibre
//     difference is 2 (paper/hyperelliptic-6-649.tex; code-derived, a regression value).
// [4] THE INJECTION CRITERION on its own: at p = 2 every odd l injects and 2 does not; at p = 3, l = 3 injects
//     (p >= 3) and l = 2 does not (2 | p - 1); at p = 7, l = 2 and 3 do not and l = 5, 7 do.

printf "Testing the component groups at p | D, the reduction of CM points and the torsion refinement...";

W1538 := {m : m in Divisors(570) | GCD(m, 570 div m) eq 1};
G3 := ComponentGroupData(15, 38, W1538, 3);
G5 := ComponentGroupData(15, 38, W1538, 5);
E := EllipticCurve("570e4");
assert ComponentGroup(G3) eq [TamagawaNumber(E, 3)] and ComponentGroup(G5) eq [TamagawaNumber(E, 5)];
assert IsTotallyDegenerate(G3, 1) and IsTotallyDegenerate(G5, 1);
pairing := [<-228, 1>, <-760, 1>, <-12, -1>, <-40, -1>];
other := [<-228, 1>, <-12, 1>, <-760, -1>, <-40, -1>];
assert ClassInComponentGroup(G3, pairing) eq [0] and ClassInComponentGroup(G5, pairing) eq [0];
assert ClassInComponentGroup(G3, other) ne [0];
cnt, npts := CMSpecialisation(G3, -228);
assert npts eq 1 and &+[cnt[v] : v in Keys(cnt)] eq 1;
B, principal, rec := TorsionOrderBoundOfCMDivisor(15, 38, W1538, other : Bound := 4, Primes := [3]);
assert not principal;

W6323 := {m : m in Divisors(1938) | GCD(m, 1938 div m) eq 1};
H2 := ComponentGroupData(6, 323, W6323, 2);
H3 := ComponentGroupData(6, 323, W6323, 3);
E := EllipticCurve("1938h1");
assert ComponentGroup(H2) eq [TamagawaNumber(E, 2)] and ComponentGroup(H3) eq [TamagawaNumber(E, 3)];
assert IsTorsionPartInjective(H2, 3) and not IsTorsionPartInjective(H2, 2);
assert IsTorsionPartInjective(H3, 3) and not IsTorsionPartInjective(H3, 2);
cnt, npts := CMSpecialisation(H2, -7752);
assert npts eq 4 and Sort([cnt[v] : v in Keys(cnt)]) eq [2, 2];
nv, edges, thick := QuotientGraph(H2);
assert nv eq 11 and &+thick eq 14 + 3;                       // the cycle has total thickness 14; three tree edges of thickness 1

W6649 := {m : m in Divisors(3894) | GCD(m, 3894 div m) eq 1};
B, principal, rec := TorsionOrderBoundOfCMDivisor(6, 649, W6649, [<-24, 1>, <-132, 1>, <-40, -1>, <-187, -1>] : Bound := 6, Primes := [3]);
assert B eq 2 and principal and rec[1][2] eq [102] and rec[1][3] eq [0];

W7 := {m : m in Divisors(14) | GCD(m, 14 div m) eq 1};
G := ComponentGroupData(14, 1, W7, 7);
assert [IsTorsionPartInjective(G, l) : l in [2, 3, 5, 7]] eq [false, false, true, true];

printf " Done!\n";
