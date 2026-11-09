// Certificates of hyperellipticity for quotients of X_0^6(379).  Read by ReadHyperellipticCertificates.
//
// X_0^6(379)/W_full, genus 3, shape A with the relation from the sharp dual test: the difference of the
// fibres P_-123 + P_-163 and P_-52 + P_-168 (rational CM points) is torsion, certified against the
// weight-3/2 cusp forms of level 4548 to precision 800, and the Jacobian is torsion-free over Q, so the
// two divisors are linearly equivalent.  The banked dual certificate data/dual/cert_6_379.m (an explicit
// principal part pairing to zero with all of S_{3/2}(4548), tests/_offline/DualCertificate6_379.m) is a
// witness for the same relation.  Provenance: PR #146, paper/hyperelliptic-6-379.tex,
// tests/DualCertificate6_379_divisor.m.
c := AssociativeArray();
c["D"] := 6; c["N"] := 379;
c["W"] := [1, 2, 3, 6, 379, 758, 1137, 2274];
c["shape"] := "A";
c["divisors"] := [ [<-123, 1>, <-163, 1>], [<-52, 1>, <-168, 1>] ];
c["dualtests"] := [ <4548, 800, "sharp", "67c92a4e1615e667faf5bb430b05c0cb"> ];
c["relations"] := [ <1, 2, 1> ];
c["torsion"] := [120, 30];
c["argument"] := "paper/torsion-counting-hyperelliptic.tex, Theorem 1.2; paper/hyperelliptic-6-379.tex";
certificates := [c];
