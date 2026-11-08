// Certificates of hyperellipticity for quotients of X_0^6(323).  Read by ReadHyperellipticCertificates.
//
// Y = X_0^6(323)/<w_6, w_19, w_34>, genus 3, shape C4: a double cover of the star E = 1938h1 branched at the
// place of degree four of discriminant -7752 (the fixed points of w_1938).  Delta = P + s^2 P - s P - s^3 P
// is the CM cycle twisted by the genus character (-8, 969); E^(-8) has root number -1 and E^(969) root
// number +1 with L(E^(969), 1) = 0 exactly (twisted modular symbol of E at level 1938), so L(E/k, chi, s) vanishes to order at least three and Delta is
// torsion (Yuan-Zhang-Zhang); the minus part of the torsion over Q(sqrt -2) and Q(sqrt 969) is {0, T}, and
// at 2 the four points lie two each on two components whose doubled difference is neither 0 nor the image
// of T in Phi_2 over either field, so Delta = 0; the Prym is simple (irreducible at 7).  Provenance: PR #154,
// paper/hyperelliptic-6-323.tex, tests/Hyperelliptic6_323.m.  Decision 2 (2026-11-05): the L-value is
// checked exactly here, where the test used the proved rank.
c := AssociativeArray();
c["D"] := 6; c["N"] := 323;
c["W"] := [1, 6, 19, 34, 51, 114, 646, 969];
c["shape"] := "C4";
c["E"] := "1938h1";
c["d"] := -7752;
c["chi"] := <-8, 969>;
c["tamagawa"] := [2, 3];
c["phi"] := [2];
c["prym"] := 7;
c["argument"] := "paper/hyperelliptic-15-38.tex, Section 1; paper/hyperelliptic-6-323.tex";
certificates := [c];
