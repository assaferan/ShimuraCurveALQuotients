// Certificates of hyperellipticity for quotients of X_0^6(187).  Read by ReadHyperellipticCertificates.
//
// X_0^6(187)/<w_3, w_11, w_17>, genus 3, shape D: a double cover of the star E = 1122i1 branched at the
// places Z_-4488 and Z_-264, both of degree two over Q(sqrt -2).  The certified Heegner heights (Cai-Shu-Tian,
// with L(1) certified or exactly zero by the twisted winding element and L'(1) certified; the rank-one
// shift recognised with denominators dividing phi(6) psi(187) = 432) give k^2 = 1 at both places and b = 7,
// 1; the configurations consistent with the trace formula are two twists of one cover, isomorphic, and
// hyperelliptic.  Shape C could not decide this curve: its Galois-stable pairing is a mixed one, and
// neither twist test applies.  #128's conic model of this curve (paper/heegner-height-models.tex) is the
// survivor.  Found 2026-11-08.
c := AssociativeArray();
c["D"] := 6; c["N"] := 187;
c["W"] := [1, 3, 11, 17, 33, 51, 187, 561];
c["shape"] := "D";
c["E"] := "1122i1";
c["phi"] := [];
c["tamagawa"] := [2, 3];
c["argument"] := "paper/heegner-height-models.tex (certified); paper/classification-certificates.tex";
certificates := [c];
