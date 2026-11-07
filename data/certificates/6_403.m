// Certificates of hyperellipticity for quotients of X_0^6(403).  Read by ReadHyperellipticCertificates.
//
// X_0^6(403)/<w_2, w_13, w_93>, genus 4, shape B: a double cover of the star (genus 2) branched at P_-3 +
// P_-403 (P_-3 an elliptic point of order three), which the three Borcherds forms of
// data/borcherds/forms_6_403.m relate to P_-43 + P_-91 and P_-88 + P_-520; the star's Jacobian is torsion-
// free, so the branch pair is canonical.  Provenance: PR #142/#143, tests/Hyperelliptic6_403.m.
c := AssociativeArray();
c["D"] := 6; c["N"] := 403;
c["W"] := [1, 2, 13, 26, 93, 186, 1209, 2418];
c["shape"] := "B";
c["divisors"] := [ [<-3, 1>, <-403, 1>], [<-43, 1>, <-91, 1>], [<-88, 1>, <-520, 1>] ];
c["forms"] := "data/borcherds/forms_6_403.m";
c["torsion"] := [60, 20];
c["argument"] := "paper/hyperelliptic-candidates.tex, Proposition prop:normalform (genus four over a genus-two star)";
certificates := [c];
