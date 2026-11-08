// Certificates of hyperellipticity for quotients of X_0^6(301).  Read by ReadHyperellipticCertificates.
//
// X_0^6(301)/<w_3, w_14, w_43>, genus 3, shape C: a double cover of the star E = 1806l1 branched at the
// places Z_-1204 (over Q(sqrt -7)) and Z_-1032 (over Q(sqrt 129)), each of degree two.  Z_-1204 - Z_-1032
// is torsion by the sharp dual test at level 3612, precision 4000, zero by its class in Phi_2 and Phi_3 and
// the reduction of E(Q)_tors, and E(Q)[2] = 0.  Found 2026-11-09 (decision 8); #128's model is not used.
c := AssociativeArray();
c["D"] := 6; c["N"] := 301;
c["W"] := [1, 3, 14, 42, 43, 129, 602, 1806];
c["shape"] := "C";
c["E"] := "1806l1";
c["divisors"] := [ [<-1204, 1>], [<-1032, 1>] ];
c["dualtests"] := [ <3612, 4000, "sharp", "db9edec77c06f391f26fdd28b0c8bbbd"> ];
c["relations"] := [ <1, 2, 1> ];
c["phi"] := [2, 3];
c["tamagawa"] := [2, 3];
c["argument"] := "paper/hyperelliptic-15-38.tex, Sections 1-4 (shape C)";
certificates := [c];
