// Certificates of hyperellipticity for quotients of X_0^21(58).  Read by ReadHyperellipticCertificates.
//
// X_0^21(58)/W_full, genus 3, shape A with relations from the 2 || N dual test (D odd): the four
// fibres F2 = P_-16 + P_-120, F4 = P_-88 + P_-100, F5 = P_-112 + P_-372, F6 = P_-168 + P_-228
// (rational CM points) have torsion differences, certified against V = (U(4) + 2) W+ at level 2436,
// precision 4000; the Jacobian is A_f with f of level 609, old at 2, and its torsion has order at
// most 2.  Provenance: PR #147, paper/hyperelliptic-21-58.tex, tests/Hyperelliptic21_58.m and
// tests/_offline/Hyperelliptic21_58.m.
c := AssociativeArray();
c["D"] := 21; c["N"] := 58;
c["W"] := [1, 2, 3, 6, 7, 14, 21, 29, 42, 58, 87, 174, 203, 406, 609, 1218];
c["shape"] := "A";
c["divisors"] := [ [<-16, 1>, <-120, 1>], [<-88, 1>, <-100, 1>], [<-112, 1>, <-372, 1>], [<-168, 1>, <-228, 1>] ];
c["dualtests"] := [ <2436, 4000, "plus2", "12739adb5a711b79c0cf2e125210db4b"> ];
c["relations"] := [ <1, 2, 1>, <1, 3, 1>, <1, 4, 1>, <2, 3, 1> ];
c["torsion"] := [60, 20];
c["argument"] := "paper/torsion-counting-hyperelliptic.tex, Theorem 1.2; paper/hyperelliptic-21-58.tex";
certificates := [c];
