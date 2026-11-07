// Certificates of hyperellipticity for quotients of X_0^77(10).  Read by ReadHyperellipticCertificates.
//
// X_0^77(10)/W_full, genus 3, shape A with relations from the 2 || N dual test (D odd): the fibres
// F2 = P_-16 + P_-100, F3 = P_-280 + P_-660, F7 = P_-760 + P_-1540, F8 = P_-840 + P_-1320 (rational
// CM points) and the quadratic cycle Q = Z(-4180) differ from F2 by torsion, certified against
// V = (U(4) + 2) W+ at level 1540, precision 8000; the Jacobian is A_f with f of level 385, old at 2,
// and its torsion has order at most 2.  Provenance: PR #148, paper/hyperelliptic-77-10.tex,
// tests/Hyperelliptic77_10.m and tests/_offline/Hyperelliptic77_10.m; the five divisors were found
// by the search of tools/proof-77-10/scan.m.
c := AssociativeArray();
c["D"] := 77; c["N"] := 10;
c["W"] := [1, 2, 5, 7, 10, 11, 14, 22, 35, 55, 70, 77, 110, 154, 385, 770];
c["shape"] := "A";
c["divisors"] := [ [<-16, 1>, <-100, 1>], [<-280, 1>, <-660, 1>], [<-760, 1>, <-1540, 1>], [<-840, 1>, <-1320, 1>], [<-4180, 1>] ];
c["dualtests"] := [ <1540, 8000, "plus2", "c19850cc5f43b50729b50a6f779ebcbe"> ];
c["relations"] := [ <1, 2, 1>, <1, 3, 1>, <1, 4, 1>, <1, 5, 1> ];
c["torsion"] := [60, 20];
c["argument"] := "paper/torsion-counting-hyperelliptic.tex, Theorem 1.2; paper/hyperelliptic-77-10.tex";
certificates := [c];
