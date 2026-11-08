// Certificates of hyperellipticity for quotients of X_0^2002(1).  Read by ReadHyperellipticCertificates.
//
// Y = X_0^2002(1)/<w_2, w_13, w_77>, genus 4, shape B: a double cover of the star (genus 2) branched at
// P_-11 + P_-91, related to P_-148 + P_-232 by the sharp dual test at level 4004, precision 2200; the
// star's Jacobian is torsion-free, so the branch pair is canonical.  The banked dual certificate
// data/dual/cert_2002_1.m (an explicit principal part pairing to zero with all of S_{3/2}(4004),
// tests/_offline/DualCertificate2002_1.m) witnesses the same relation.  Provenance: PR #139/#140,
// paper/star-2002-proof.tex, tests/DualCertificate2002_1_divisor.m.
c := AssociativeArray();
c["D"] := 2002; c["N"] := 1;
c["W"] := [1, 2, 13, 26, 77, 154, 1001, 2002];
c["shape"] := "B";
c["divisors"] := [ [<-11, 1>, <-91, 1>], [<-148, 1>, <-232, 1>] ];
c["dualtests"] := [ <4004, 2200, "sharp", "df5e45ff3c1a9361ce129ab1b2dbdee6"> ];
c["relations"] := [ <1, 2, 1> ];
c["torsion"] := [60, 20];
c["argument"] := "paper/hyperelliptic-candidates.tex, Proposition prop:normalform (genus four over a genus-two star)";
certificates := [c];
