// Certificates of hyperellipticity for quotients of X_0^55(7).  Read by ReadHyperellipticCertificates.
//
// X_0^55(7)/<w_35, w_55>, genus 3, shape C: a double cover of the star E = 385a1 branched at the places
// Z_-1540 (over Q(sqrt -35)) and Z_-20 (over Q(sqrt 5)), each of degree two.  Z_-1540 - Z_-20 is torsion
// by the sharp dual test at level 1540, precision 8000, on the a_3(E)-eigenspace of T(9) (the whole space
// W has dimension 7), and zero by its class in Phi_5 and Phi_11 and the reduction of E(Q)_tors; the Prym is
// simple (irreducible at 13), so Y is the untwisted cover.  Found 2026-11-08 (decision 8); #128's model
// of this curve is not used.
c := AssociativeArray();
c["D"] := 55; c["N"] := 7;
c["W"] := [1, 35, 55, 77];
c["shape"] := "C";
c["E"] := "385a1";
c["divisors"] := [ [<-1540, 1>], [<-20, 1>] ];
c["dualtests"] := [ <1540, 8000, "sharp", "c19850cc5f43b50729b50a6f779ebcbe", 3> ];
c["relations"] := [ <1, 2, 1> ];
c["phi"] := [5, 11];
c["tamagawa"] := [5, 11];
c["prym"] := 13;
c["argument"] := "paper/hyperelliptic-15-38.tex, Sections 1-4 (shape C)";
certificates := [c];
