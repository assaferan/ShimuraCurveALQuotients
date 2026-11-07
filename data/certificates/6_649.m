// Certificates of hyperellipticity for quotients of X_0^6(649).  Read by ReadHyperellipticCertificates.
//
// X_0^6(649)/W_full, genus 3, shape A with relations from the dual test and the torsion bound lowered by
// a component group: the three fibres D_1 = P_-24 + P_-132, D_2 = P_-40 + P_-187, D_3 = P_-52 + P_-660
// (rational CM points, from the Shimizu screen) have torsion differences, certified by the sharp dual test
// at level 7788, precision 1200 (W of dimension 53).  The rational torsion of the Jacobian has order
// dividing 6; its 3-part injects into Phi_3 = Z/102, where the differences have class zero, so they lie
// in a group of order at most 2, against three divisors.  Provenance: PR #156 and #158,
// paper/hyperelliptic-6-649.tex, paper/component-groups.tex, tests/Hyperelliptic6_649.m and
// tests/_offline/Hyperelliptic6_649.m.
c := AssociativeArray();
c["D"] := 6; c["N"] := 649;
c["W"] := [1, 2, 3, 6, 11, 22, 33, 59, 66, 118, 177, 354, 649, 1298, 1947, 3894];
c["shape"] := "A";
c["divisors"] := [ [<-24, 1>, <-132, 1>], [<-40, 1>, <-187, 1>], [<-52, 1>, <-660, 1>] ];
c["dualtests"] := [ <7788, 1200, "sharp", "6898a417123c25e1e6430b61278c8ae1"> ];
c["relations"] := [ <1, 2, 1>, <1, 3, 1>, <2, 3, 1> ];
c["torsion"] := [120, 30];
c["phi"] := [3];
c["argument"] := "paper/torsion-counting-hyperelliptic.tex, Theorem 1.2; paper/component-groups.tex; paper/hyperelliptic-6-649.tex";
certificates := [c];
