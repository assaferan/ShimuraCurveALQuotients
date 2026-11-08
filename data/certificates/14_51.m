// Certificates of hyperellipticity for quotients of X_0^14(51).  Read by ReadHyperellipticCertificates.
//
// X_0^14(51)/<w_6, w_14, w_34>, genus 3, shape C: a double cover of the star E = 238a1 (the star's form is
// old at 3) branched at the rational CM points P_-1428, P_-408, P_-168, P_-8.  P_-1428 + P_-408 - P_-168 - P_-8
// is torsion by the sharp dual test at level 1428, precision 4000, zero by its class in Phi_2 and Phi_7 and
// the reduction of E(Q)_tors; the Prym of Y -> E is simple (irreducible at 5), so Y is the untwisted cover.
// Not one of the 34: the filters leave it open, and the special fibre at 3 cannot rule it out, since
// Z = X_0^14(17)/<w_14, w_34> is hyperelliptic mod 3 although not over Q.  Found 2026-11-08.
c := AssociativeArray();
c["D"] := 14; c["N"] := 51;
c["W"] := [1, 6, 14, 21, 34, 51, 119, 714];
c["shape"] := "C";
c["E"] := "238a1";
c["divisors"] := [ [<-1428, 1>, <-408, 1>], [<-168, 1>, <-8, 1>] ];
c["dualtests"] := [ <1428, 4000, "sharp", "71f057cdb3a63253525461f61f8d435e"> ];
c["relations"] := [ <1, 2, 1> ];
c["phi"] := [2, 7];
c["tamagawa"] := [2, 7];
c["prym"] := 5;
c["argument"] := "paper/hyperelliptic-15-38.tex, Sections 1-4 (shape C)";
certificates := [c];
