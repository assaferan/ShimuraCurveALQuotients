// Certificates of hyperellipticity for quotients of X_0^38(15).  Read by ReadHyperellipticCertificates.
//
// X_0^38(15)/<w_3, w_10, w_38>, genus 3, shape C: a double cover of the star E = 190a1 (the star's form is
// old at 3) branched at the rational CM points P_-24, P_-20 and the place Z_-2280 of degree two.
// P_-24 + P_-20 - Z_-2280 is torsion by the sharp dual test at level 1140, precision 4000, zero by its class
// in Phi_2 and Phi_19 and the reduction of E(Q)_tors, and E(Q)[2] = 0.  Not one of the 34: the filters
// leave it open, and the special fibre at 3 cannot rule it out, since Z = X_0^38(5)/<w_10, w_38> is
// hyperelliptic mod 3 although not over Q.  Found 2026-11-08.
c := AssociativeArray();
c["D"] := 38; c["N"] := 15;
c["W"] := [1, 3, 10, 30, 38, 95, 114, 285];
c["shape"] := "C";
c["E"] := "190a1";
c["divisors"] := [ [<-24, 1>, <-20, 1>], [<-2280, 1>] ];
c["dualtests"] := [ <1140, 4000, "sharp", "1fbd9ca124791e88a17abc214aebf629"> ];
c["relations"] := [ <1, 2, 1> ];
c["phi"] := [2, 19];
c["tamagawa"] := [2, 19];
c["argument"] := "paper/hyperelliptic-15-38.tex, Sections 1-4 (shape C)";
certificates := [c];
