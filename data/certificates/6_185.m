// Certificates of hyperellipticity for quotients of X_0^6(185).  Read by ReadHyperellipticCertificates.
//
// X_0^6(185)/<w_2, w_15, w_111>, genus 3, shape C: a double cover of the star E = 1110m1 branched at the
// rational CM points P_-555, P_-40 and the place Z_-4440 of degree two.  P_-555 + P_-40 - Z_-4440 is
// torsion by the sharp dual test at level 2220, precision 4000 (W of dimension 1), and zero by its class
// in Phi_2 and Phi_3 and the reduction of E(Q)_tors; E(Q)[2] = 0, so there is no twist to exclude.
// Found 2026-11-08 (decision 8); #128's model of this curve is not used.
c := AssociativeArray();
c["D"] := 6; c["N"] := 185;
c["W"] := [1, 2, 15, 30, 111, 185, 222, 370];
c["shape"] := "C";
c["E"] := "1110m1";
c["divisors"] := [ [<-555, 1>, <-40, 1>], [<-4440, 1>] ];
c["dualtests"] := [ <2220, 4000, "sharp", "c66c2be2af43b7684bf605f3d34fbc3d"> ];
c["relations"] := [ <1, 2, 1> ];
c["phi"] := [2, 3];
c["tamagawa"] := [2, 3];
c["argument"] := "paper/hyperelliptic-15-38.tex, Sections 1-4 (shape C)";
certificates := [c];
