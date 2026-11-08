// Certificates of hyperellipticity for quotients of X_0^57(10).  Read by ReadHyperellipticCertificates.
//
// X_0^57(10)/<w_5, w_6, w_19>, genus 3, shape C: a double cover of the star E = 285a1 (the form is old at 2)
// branched at the rational CM points P_-760, P_-4 and the place Z_-1140 of degree two.
// P_-760 + P_-4 - Z_-1140 is torsion by the sharp dual test at level 1140, precision 4000 (P_-4 has weight
// 1/4: it is fixed by w_2), and zero by its class in Phi_3 and Phi_19 and the reduction of E(Q)_tors; the
// Prym is simple (irreducible at 7), so Y is the untwisted cover.  Found 2026-11-08 (decision 8); #128's
// model of this curve is not used.
c := AssociativeArray();
c["D"] := 57; c["N"] := 10;
c["W"] := [1, 5, 6, 19, 30, 95, 114, 570];
c["shape"] := "C";
c["E"] := "285a1";
c["divisors"] := [ [<-760, 1>, <-4, 1>], [<-1140, 1>] ];
c["dualtests"] := [ <1140, 4000, "sharp", "1fbd9ca124791e88a17abc214aebf629"> ];
c["relations"] := [ <1, 2, 1> ];
c["phi"] := [3, 19];
c["tamagawa"] := [3, 19];
c["prym"] := 7;
c["argument"] := "paper/hyperelliptic-15-38.tex, Sections 1-4 (shape C)";
certificates := [c];
