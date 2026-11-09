// Certificates of hyperellipticity for quotients of X_0^1326(1).  Read by ReadHyperellipticCertificates.
//
// X_0^1326(1)/<w_17, w_78>, genus 3, shape D over an intermediate quotient: its Atkin-Lehner group has index
// four, and Y is a double cover of the genus-one quotient E = X_0^1326(1)/{1, 6, 13, 17, 78, 102, 221, 1326} =
// 1326b2 (identified by its traces up to the Sturm bound and its component groups), branched at the places
// Z_-24 and Z_-408 over Q(sqrt 2); the certified Heegner heights and the trace formula leave only
// hyperelliptic covers (paper/heegner-height-models.tex, certified mode).
c := AssociativeArray();
c["D"] := 1326; c["N"] := 1;
c["W"] := [1, 17, 78, 1326];
c["Wp"] := [1, 6, 13, 17, 78, 102, 221, 1326];
c["shape"] := "D";
c["E"] := "1326b2";
c["phi"] := [];
c["tamagawa"] := [2, 3, 13, 17];
c["argument"] := "paper/heegner-height-models.tex (certified); paper/classification-certificates.tex";
certificates := [c];
