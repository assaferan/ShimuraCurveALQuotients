// Certificates of hyperellipticity for quotients of X_0^2418(1).  Read by ReadHyperellipticCertificates.
//
// X_0^2418(1)/<w_2, w_13, w_93>, genus 4, shape B: a double cover of the star (genus 2) branched at P_-312
// + P_-403, which the three Borcherds forms of data/borcherds/forms_2418_1.m relate to P_-19 + P_-67 and
// P_-195 + P_-267; the star's Jacobian is torsion-free, so the branch pair is canonical.  Provenance: PR
// #142/#143, tests/Hyperelliptic2418_1.m.
c := AssociativeArray();
c["D"] := 2418; c["N"] := 1;
c["W"] := [1, 2, 13, 26, 93, 186, 1209, 2418];
c["shape"] := "B";
c["divisors"] := [ [<-312, 1>, <-403, 1>], [<-19, 1>, <-67, 1>], [<-195, 1>, <-267, 1>] ];
c["forms"] := "data/borcherds/forms_2418_1.m";
c["torsion"] := [60, 20];
c["argument"] := "paper/hyperelliptic-candidates.tex, Proposition prop:normalform (genus four over a genus-two star)";
certificates := [c];
