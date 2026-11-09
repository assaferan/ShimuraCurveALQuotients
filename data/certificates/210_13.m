// Certificates of hyperellipticity for quotients of X_0^210(13).  Read by ReadHyperellipticCertificates.
//
// X_0^210(13)/<w_2, w_21, w_65>, genus 3, shape D over an intermediate quotient: its Atkin-Lehner group has
// index four, and Y is a double cover of the genus-one quotient E = X_0^210(13)/W' = 2730u1, W' below
// (identified by its traces up to the Sturm bound and its component groups), branched at the places Z_-420
// and Z_-1092 over Q(i); the certified Heegner heights and the trace formula leave only hyperelliptic covers
// (paper/heegner-height-models.tex, certified mode).
c := AssociativeArray();
c["D"] := 210; c["N"] := 13;
c["W"] := [1, 2, 21, 42, 65, 130, 1365, 2730];
c["Wp"] := [1, 2, 5, 10, 13, 21, 26, 42, 65, 105, 130, 210, 273, 546, 1365, 2730];
c["shape"] := "D";
c["E"] := "2730u1";
c["phi"] := [];
c["tamagawa"] := [2, 3, 5, 7];
c["argument"] := "paper/heegner-height-models.tex (certified); paper/classification-certificates.tex";
certificates := [c];
