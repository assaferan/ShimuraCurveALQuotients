// Certificates of hyperellipticity for quotients of X_0^6(247).  Read by ReadHyperellipticCertificates.
//
// X_0^6(247)/<w_2, w_19, w_39>, genus 3, shape C: a double cover of the star E = 1482i1 branched at the
// rational CM points P_-52, P_-3 and the place Z_-456 of degree two.  P_-52 + P_-3 - Z_-456 is torsion by
// the sharp dual test at level 2964, precision 4000, zero by its class in Phi_2 and Phi_3 and the
// reduction of E(Q)_tors, and E(Q)[2] = 0.  Found 2026-11-08 (decision 8); #128's model is not used.
c := AssociativeArray();
c["D"] := 6; c["N"] := 247;
c["W"] := [1, 2, 19, 38, 39, 78, 741, 1482];
c["shape"] := "C";
c["E"] := "1482i1";
c["divisors"] := [ [<-52, 1>, <-3, 1>], [<-456, 1>] ];
c["dualtests"] := [ <2964, 4000, "sharp", "c643f46ea259ad22b0acf6ab64ef361e"> ];
c["relations"] := [ <1, 2, 1> ];
c["phi"] := [2, 3];
c["tamagawa"] := [2, 3];
c["argument"] := "paper/hyperelliptic-15-38.tex, Sections 1-4 (shape C)";
certificates := [c];
