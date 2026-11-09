// Certificates of hyperellipticity for quotients of X_0^6(385).  Read by ReadHyperellipticCertificates.
//
// Y = X_0^6(385)/H, genus 3, shape C: a double cover of the star E = 2310q1 branched at P_-40, P_-1540
// (rational points) and Z_-9240 (a place of degree two).  P_-40 + P_-1540 - Z_-9240 is torsion by the
// sharp dual test at level 4620, precision 4000 (W of dimension 1), and zero because its class in
// Phi_2 = Z/20 is zero and E(Q)_tors = Z/2 leaves the identity component at 2; the Prym of Y -> E is
// simple (L-polynomial irreducible at 19), so Y is the untwisted cover.  Provenance: PR #152,
// paper/hyperelliptic-6-385.tex, tests/Hyperelliptic6_385.m and tests/_offline/Hyperelliptic6_385.m.
c := AssociativeArray();
c["D"] := 6; c["N"] := 385;
c["W"] := [1, 5, 6, 11, 14, 21, 30, 55, 66, 70, 105, 154, 231, 330, 770, 1155];
c["shape"] := "C";
c["E"] := "2310q1";
c["divisors"] := [ [<-40, 1>, <-1540, 1>], [<-9240, 1>] ];
c["dualtests"] := [ <4620, 4000, "sharp", "02387f1c3ec26965c0eb9a19f73f44bf"> ];
c["relations"] := [ <1, 2, 1> ];
c["phi"] := [2];
c["tamagawa"] := [2, 3];
c["prym"] := 19;
c["argument"] := "paper/hyperelliptic-15-38.tex, Sections 1-4; paper/hyperelliptic-6-385.tex";
certificates := [c];
