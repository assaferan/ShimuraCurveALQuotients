// Certificates of hyperellipticity for quotients of X_0^6(455).  Read by ReadHyperellipticCertificates.
//
// X_0^6(455)/<w_6, w_7, w_10, w_13>, genus 3, shape D: a double cover of the star E = 2730y1 branched at the
// places Z_-1560 (over Q(sqrt 10)) and Z_-10920 (over Q(sqrt -26)); the certified Heegner heights and the
// trace formula leave only hyperelliptic covers (paper/heegner-height-models.tex, certified mode).
c := AssociativeArray();
c["D"] := 6; c["N"] := 455;
c["W"] := [1, 6, 7, 10, 13, 15, 42, 70, 78, 91, 105, 130, 195, 546, 910, 1365];
c["shape"] := "D";
c["E"] := "2730y1";
c["phi"] := [];
c["tamagawa"] := [2, 3];
c["argument"] := "paper/heegner-height-models.tex (certified); paper/classification-certificates.tex";
certificates := [c];
