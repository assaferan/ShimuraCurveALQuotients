// Certificates of hyperellipticity for quotients of X_0^2030(1).  Read by ReadHyperellipticCertificates.
//
// X_0^2030(1)/<w_10, w_14, w_58>, genus 4, shape B: a double cover of the star (genus 2) branched at P_-8
// + P_-280, which the three Borcherds forms of data/borcherds/forms_2030_1.m relate to P_-43 + P_-163 and
// P_-148 + P_-532; the star's Jacobian is torsion-free, so the branch pair is canonical.  Provenance: PR
// #142/#143, tests/Hyperelliptic2030_1.m.
c := AssociativeArray();
c["D"] := 2030; c["N"] := 1;
c["W"] := [1, 10, 14, 35, 58, 145, 203, 2030];
c["shape"] := "B";
c["divisors"] := [ [<-8, 1>, <-280, 1>], [<-43, 1>, <-163, 1>], [<-148, 1>, <-532, 1>] ];
c["forms"] := "data/borcherds/forms_2030_1.m";
c["torsion"] := [60, 20];
c["argument"] := "paper/hyperelliptic-candidates.tex, Proposition prop:normalform (genus four over a genus-two star)";
certificates := [c];
