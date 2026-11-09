// Certificates of hyperellipticity for quotients of X_0^15(46).  Read by ReadHyperellipticCertificates.
//
// X_0^15(46)/<w_5, w_6, w_23>, genus 3, shape C: a double cover of the star E = 345f1 (the form is old at 2)
// branched at the rational CM points P_-1380, P_-60, P_-40, P_-15.  Of the three pairings, only
// {P_-1380, P_-60} | {P_-40, P_-15} is certified torsion: by the 2 || N dual test for odd discriminants at
// level 1380, precision 4000 (the divisor pushed to X_0^15(23)* and tested against the plus space; the
// sharp test does not certify it even at precision 10000); it is zero by its class in Phi_3 and Phi_5 and
// the reduction of E(Q)_tors, and E(Q)[2] = 0.  Found 2026-11-08 (decision 8); #128's model is not used.
c := AssociativeArray();
c["D"] := 15; c["N"] := 46;
c["W"] := [1, 5, 6, 23, 30, 115, 138, 690];
c["shape"] := "C";
c["E"] := "345f1";
c["divisors"] := [ [<-1380, 1>, <-60, 1>], [<-40, 1>, <-15, 1>] ];
c["dualtests"] := [ <1380, 4000, "odd", "8260c52f4c625395b86d5a940fa80371"> ];
c["relations"] := [ <1, 2, 1> ];
c["phi"] := [3, 5];
c["tamagawa"] := [3, 5];
c["argument"] := "paper/hyperelliptic-15-38.tex, Sections 1-4 (shape C); paper/hyperelliptic-39-34.tex (the odd variant)";
certificates := [c];
