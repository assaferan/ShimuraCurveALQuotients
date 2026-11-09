// Certificates of hyperellipticity for quotients of X_0^3390(1).  Read by ReadHyperellipticCertificates.
//
// X_0^3390(1)/<w_5, w_6, w_113>, genus 4, shape B: a double cover of the star (genus 2) branched at P_-3 +
// P_-40, related to the fibres P_-43 + P_-67 and P_-123 + P_-168 by the sharp dual test at level 6780,
// precision 1100 (W of dimension 51); the star's Jacobian is torsion-free, so the branch pair is
// canonical.  Provenance: PR #155, paper/hyperelliptic-3390-1.tex, tests/Hyperelliptic3390_1.m and
// tests/_offline/Hyperelliptic3390_1.m.
c := AssociativeArray();
c["D"] := 3390; c["N"] := 1;
c["W"] := [1, 5, 6, 30, 113, 565, 678, 3390];
c["shape"] := "B";
c["divisors"] := [ [<-3, 1>, <-40, 1>], [<-43, 1>, <-67, 1>], [<-123, 1>, <-168, 1>] ];
c["dualtests"] := [ <6780, 1100, "sharp", "65988ce2d6dccb5f13d9c6dd61e2f125"> ];
c["relations"] := [ <1, 2, 1>, <1, 3, 1>, <2, 3, 1> ];
c["torsion"] := [60, 20];
c["argument"] := "paper/hyperelliptic-candidates.tex, Proposition prop:normalform (genus four over a genus-two star)";
certificates := [c];
