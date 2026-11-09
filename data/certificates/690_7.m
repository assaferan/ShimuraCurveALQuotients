// Certificates of hyperellipticity for quotients of X_0^690(7).  Read by ReadHyperellipticCertificates.
//
// X_0^690(7)/<w_6, w_10, w_14, w_23>, genus 3, shape D: a double cover of the star E = 4830bj1 branched at the
// rational CM points P_-3, P_-115, P_-280, P_-420; the certified Heegner heights and the
// trace formula leave only hyperelliptic covers (paper/heegner-height-models.tex, certified mode).
c := AssociativeArray();
c["D"] := 690; c["N"] := 7;
c["W"] := [1, 6, 10, 14, 15, 21, 23, 35, 138, 210, 230, 322, 345, 483, 805, 4830];
c["shape"] := "D";
c["E"] := "4830bj1";
c["phi"] := [];
c["tamagawa"] := [2, 3, 5, 23];
c["argument"] := "paper/heegner-height-models.tex (certified); paper/classification-certificates.tex";
certificates := [c];
