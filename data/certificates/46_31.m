// Certificates of hyperellipticity for quotients of X_0^46(31).  Read by ReadHyperellipticCertificates
// (HyperellipticCertificates.m) and checked by VerifyHyperellipticCertificate.
//
// X_0^46(31)/W_full, genus 3, shape A (torsion counting): the four effective degree-two divisors
//     P_-3 + P_-403,  P_-24 + P_-232,  P_-115 + P_-372,  Z(-292)  (a cycle of degree two)
// are pairwise related by the six Borcherds forms banked in data/borcherds/forms_46_31.m, and the
// rational torsion of the Jacobian has order at most 2.  Provenance: the divisors were found by the
// search of tools/hyperellipticstar.m over the 414 degree-two divisors supported on the rational
// CM points and quadratic cycles (PR #133); the argument is Theorem 1.2 of
// paper/torsion-counting-hyperelliptic.tex.  tests/Hyperelliptic46_31.m checks the same inputs.
c := AssociativeArray();
c["D"] := 46; c["N"] := 31;
c["W"] := [1, 2, 23, 31, 46, 62, 713, 1426];
c["shape"] := "A";
c["divisors"] := [ [<-3, 1>, <-403, 1>], [<-24, 1>, <-232, 1>], [<-115, 1>, <-372, 1>], [<-292, 1>] ];
c["forms"] := "data/borcherds/forms_46_31.m";
c["torsion"] := [60, 6];          // MaxPrime and NumPrimes of JacobianTorsionBound
c["argument"] := "paper/torsion-counting-hyperelliptic.tex, Theorem 1.2";
certificates := [c];
