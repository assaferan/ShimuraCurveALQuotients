// Certificates of hyperellipticity for quotients of X_0^22(51).  Read by ReadHyperellipticCertificates.
//
// X_0^22(51)/<w_6, w_17, w_22>, genus 3, shape D: a double cover of the star E = 1122f1 branched at the
// places Z_-4488 and Z_-264, both of degree two over Q(sqrt -2).  The certified Heegner heights (the rank-one
// shift recognised with denominators dividing phi(22) psi(51) = 720) give k^2 = 1 at both places and b = 1,
// 1; one configuration is consistent with the trace formula, and it is hyperelliptic.  Shape C certified
// the pairing Z_-4488 | Z_-264 torsion and zero but could not settle the twist (no simple Prym, and no odd
// prime of D with an even cyclic component group).  Found 2026-11-08.
c := AssociativeArray();
c["D"] := 22; c["N"] := 51;
c["W"] := [1, 6, 17, 22, 33, 102, 374, 561];
c["shape"] := "D";
c["E"] := "1122f1";
c["phi"] := [];
c["tamagawa"] := [2, 11];
c["argument"] := "paper/heegner-height-models.tex (certified); paper/classification-certificates.tex";
certificates := [c];
