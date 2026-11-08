// Certificates of hyperellipticity for quotients of X_0^6(253).  Read by ReadHyperellipticCertificates.
//
// X_0^6(253)/<w_2, w_23, w_33>, genus 3, shape C: a double cover of the star E = 1518p1 branched at the
// rational CM points P_-1012, P_-88 and the place Z_-552 of degree two.  P_-1012 + P_-88 - Z_-552 is torsion
// by the sharp dual test at level 3036, precision 4000, zero by its class in Phi_2 and Phi_3 and the
// reduction of E(Q)_tors, and E(Q)[2] = 0.  Found 2026-11-08 (decision 8); #128's model is not used.
c := AssociativeArray();
c["D"] := 6; c["N"] := 253;
c["W"] := [1, 2, 23, 33, 46, 66, 759, 1518];
c["shape"] := "C";
c["E"] := "1518p1";
c["divisors"] := [ [<-1012, 1>, <-88, 1>], [<-552, 1>] ];
c["dualtests"] := [ <3036, 4000, "sharp", "7809579f985de8da83986129afc53899"> ];
c["relations"] := [ <1, 2, 1> ];
c["phi"] := [2, 3];
c["tamagawa"] := [2, 3];
c["argument"] := "paper/hyperelliptic-15-38.tex, Sections 1-4 (shape C)";
certificates := [c];
