// Certificates of hyperellipticity for quotients of X_0^4130(1).  Read by ReadHyperellipticCertificates.
//
// X_0^4130(1)/W_full, genus 3, shape A with relations from the dual test: the four fibres
//     P_-35 + P_-532,  P_-163 + P_-595,  P_-315 + P_-520,  P_-340 + P_-420
// (rational CM points, found by the Shimizu screen) have pairwise torsion differences, certified by
// the sharp dual test against the weight-3/2 cusp forms of level 8260 to precision 4000 (W has
// dimension 3), and the Jacobian is torsion-free over Q, so two of them already give a g^1_2.
// Provenance: PR #157, paper/hyperelliptic-4130-1.tex, tests/_offline/Hyperelliptic4130_1.m; the
// basis took 61 h on lovelace (tools/proof-4130-1/basisS.m).
c := AssociativeArray();
c["D"] := 4130; c["N"] := 1;
c["W"] := [1, 2, 5, 7, 10, 14, 35, 59, 70, 118, 295, 413, 590, 826, 2065, 4130];
c["shape"] := "A";
c["divisors"] := [ [<-35, 1>, <-532, 1>], [<-163, 1>, <-595, 1>], [<-315, 1>, <-520, 1>], [<-340, 1>, <-420, 1>] ];
c["dualtests"] := [ <8260, 4000, "sharp", "3d024cd48b1bb7aa1aaa5839e0d8ac49"> ];
c["relations"] := [ <i, j, 1> : i, j in [1..4] | i lt j ];
c["torsion"] := [120, 30];
c["argument"] := "paper/torsion-counting-hyperelliptic.tex, Theorem 1.2; paper/hyperelliptic-4130-1.tex";
certificates := [c];
