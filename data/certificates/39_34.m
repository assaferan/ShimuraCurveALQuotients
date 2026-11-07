// Certificates of hyperellipticity for quotients of X_0^39(34).  Read by ReadHyperellipticCertificates.
//
// X_0^39(34)/W_full, genus 3, shape A with relations from the 2 || N dual tests (D odd): the fibres
// F = P_-84 + P_-132, F1 = P_-15 + P_-52, F7 = P_-1092 + P_-372 (rational CM points) and the
// quadratic cycles Q1 = Z(-3432), Q2 = Z(-1560): F7 - F, Q1 - F, Q2 - F are certified against
// V = (U(4) + 2) W+ at level 2652, precision 4000, and F1 - F (odd discriminants) by pushing it to
// X_0^39(17)* and testing against W+.  The Jacobian is A_f with f of level 663, old at 2, and its
// torsion has order at most 4, against five divisors.  Provenance: PR #149,
// paper/hyperelliptic-39-34.tex, tests/Hyperelliptic39_34.m and tests/_offline/Hyperelliptic39_34.m.
c := AssociativeArray();
c["D"] := 39; c["N"] := 34;
c["W"] := [1, 2, 3, 6, 13, 17, 26, 34, 39, 51, 78, 102, 221, 442, 663, 1326];
c["shape"] := "A";
c["divisors"] := [ [<-84, 1>, <-132, 1>], [<-15, 1>, <-52, 1>], [<-1092, 1>, <-372, 1>], [<-3432, 1>], [<-1560, 1>] ];
c["dualtests"] := [ <2652, 4000, "plus2", "3664b6901e637a006f7bfd25f161ae9e">, <2652, 4000, "odd", "3664b6901e637a006f7bfd25f161ae9e"> ];
c["relations"] := [ <3, 1, 1>, <4, 1, 1>, <5, 1, 1>, <2, 1, 2> ];
c["torsion"] := [60, 20];
c["argument"] := "paper/torsion-counting-hyperelliptic.tex, Theorem 1.2; paper/hyperelliptic-39-34.tex";
certificates := [c];
