// Certificates of hyperellipticity for quotients of X_0^714(5).  Read by ReadHyperellipticCertificates.
//
// Y = X_0^714(5)/H, genus 3, shape C: a double cover of the star E = 3570x2 branched at P_-340, P_-840
// (rational points) and Z_-14280 (a place of degree two).  P_-340 + P_-840 - Z_-14280 is torsion by the
// sharp dual test at level 7140, precision 4000 (W of dimension 1), and zero because its class in
// Phi_3, Phi_7 and Phi_17 is zero and every nonzero point of E(Q)_tors = (Z/2)^2 leaves the identity
// component at one of them; the Prym of Y -> E is simple (L-polynomial irreducible at 11), so Y is the
// untwisted cover.  Provenance: PR #153, paper/hyperelliptic-714-5.tex, tests/Hyperelliptic714_5.m and
// tests/_offline/Hyperelliptic714_5.m; the basis took 99 h on lovelace.
c := AssociativeArray();
c["D"] := 714; c["N"] := 5;
c["W"] := [1, 3, 10, 14, 17, 30, 35, 42, 51, 105, 170, 238, 510, 595, 714, 1785];
c["shape"] := "C";
c["E"] := "3570x2";
c["divisors"] := [ [<-340, 1>, <-840, 1>], [<-14280, 1>] ];
c["dualtests"] := [ <7140, 4000, "sharp", "e66db0fecb0b683c3fb0c2fb1ef9eaeb"> ];
c["relations"] := [ <1, 2, 1> ];
c["phi"] := [3, 7, 17];
c["tamagawa"] := [3, 7, 17];
c["prym"] := 11;
c["argument"] := "paper/hyperelliptic-15-38.tex, Sections 1-4; paper/hyperelliptic-714-5.tex";
certificates := [c];
