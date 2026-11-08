// Certificates of hyperellipticity for quotients of X_0^39(11).  Read by ReadHyperellipticCertificates.
//
// X_0^39(11)/<w_11, w_39>, genus 3, shape C: a double cover of the star E = 429b1 branched at the place
// Z_-132 of degree two (over Q(sqrt -11)) and the place Z_-52 of degree two (over Q(sqrt 13)).
// Z_-132 - Z_-52 is torsion by the sharp dual test at level 1716, precision 4000, and zero by its class
// in Phi_3 and Phi_13 and the reduction of E(Q)_tors; the Prym of Y -> E is simple (irreducible at 5),
// so Y is the untwisted cover.  Found 2026-11-08 (decision 8); #128's model of this curve
// (paper/heegner-height-models.tex) is not used.
c := AssociativeArray();
c["D"] := 39; c["N"] := 11;
c["W"] := [1, 11, 39, 429];
c["shape"] := "C";
c["E"] := "429b1";
c["divisors"] := [ [<-132, 1>], [<-52, 1>] ];
c["dualtests"] := [ <1716, 4000, "sharp", "4cc75c550f6ce856fbf180a95cd5a959"> ];
c["relations"] := [ <1, 2, 1> ];
c["phi"] := [3, 13];
c["tamagawa"] := [3, 13];
c["prym"] := 5;
c["argument"] := "paper/hyperelliptic-15-38.tex, Sections 1-4 (shape C)";
certificates := [c];
