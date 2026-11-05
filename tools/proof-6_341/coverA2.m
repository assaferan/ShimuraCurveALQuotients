SetColumns(0);
P<x> := PolynomialRing(Rationals());
CA := HyperellipticCurve(P![144,-120,265,-82,133,-12,36]); JA := Jacobian(CA);
f := HyperellipticPolynomials(CA); F := Factorisation(f);
q4 := [g[1] : g in F | Degree(g[1]) eq 4][1];
G := GaloisGroup(q4); printf "Galois group of the quartic factor: order %o, transitive %o\n", #G, IsTransitive(G);
pts := Points(CA : Bound := 5000); P0 := pts[1]; inf := [pt : pt in pts | pt[3] eq 0];
q := [g[1] : g in F | Degree(g[1]) eq 2][1];
c := (JA![q, 0]) + (JA![inf[1], P0]) + (JA![inf[2], P0]);
printf "class c = [W1 + W2 - 2 P0] = %o\n", c;
// independent check through the Mordell-Weil group
MW, mw, flag1, flag2 := MordellWeilGroupGenus2(JA);
printf "Mordell-Weil group: %o   (proved: %o %o)\n", MW, flag1, flag2;
v := c @@ mw;
printf "coordinates of c in MW: %o;  c in 2 J(Q): %o\n", Eltseq(v), forall{i : i in [1..#Generators(MW)] | IsEven(Eltseq(v)[i]) or (Order(MW.i) ne 0 and Order(MW.i) eq 1)};
// explicit: is c - 2*s = 0 for some s? equivalently v in 2 MW
printf "v in 2*MW: %o\n", v in sub<MW | [2*MW.i : i in [1..#Generators(MW)]]>;
// the same for C with its torsion point t
C := HyperellipticCurve(P![9,-33,9,29,-12,6], P![0,1,1]); Cs := SimplifiedModel(C); J := Jacobian(Cs);
MWc, mwc, f1, f2 := MordellWeilGroupGenus2(J);
T, mT := TorsionSubgroup(J); t := mT(T.1);
printf "on C: Mordell-Weil %o (proved %o %o); t in 2*MW: %o\n", MWc, f1, f2, (t @@ mwc) in sub<MWc | [2*MWc.i : i in [1..#Generators(MWc)]]>;
quit;
