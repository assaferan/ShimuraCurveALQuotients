// Certificates for quotients of X_0^6(35).  Read by ReadHyperellipticCertificates.
//
// Both genus-3 quotients X_0^6(35)/<w_3, w_5> and X_0^6(35)/<w_5, w_7> are NOT hyperelliptic, by shape V. Each is
// a Klein-four cover of the genus-zero star whose three intermediate quotients have genus one, and a double cover
// of E = X_0^6(35)/<w_3, w_5, w_7> (Jacobian 30a6) branched at the places Z_-420 and Z_-84, of degree two over
// Q(sqrt -5) and Q(sqrt 7).  So Z_-420 | Z_-84 is the only Galois-stable pairing, both places are fibres of E over
// the star, and a hyperelliptic cover would contain a double cover of the star of genus zero.  The two curves
// differ by a rational 2-torsion twist, which the argument does not use.  The stored models, which rest on the
// m = 0 term at a composite level, are plane quartics; these certificates do not use them.  Found 2026-11-09.
//
// The first curve has a second certificate, shape T, kept as an independent proof: the hyperelliptic cover would
// be E(sqrt(c x)) with div x = Z_-420 - Z_-84, and at 3 the dual graph of Y has one singular point over the
// self-intersection of thickness two of E, where E(sqrt(c x)) has two.  It does not decide the second curve,
// whose dual graph at 3 agrees with E(sqrt(c x)).  The stage uses the first certificate it accepts.
c := AssociativeArray();
c["D"] := 6; c["N"] := 35;
c["W"] := [1, 3, 5, 15];
c["shape"] := "V";
c["over"] := 7;
c["argument"] := "paper/classification-certificates.tex, Section 7 (shape V)";
t := AssociativeArray();
t["D"] := 6; t["N"] := 35;
t["W"] := [1, 3, 5, 15];
t["shape"] := "T";
t["over"] := 7;
t["p"] := 3;
t["argument"] := "paper/classification-certificates.tex, Section 7 (shape T)";
c2 := AssociativeArray();
c2["D"] := 6; c2["N"] := 35;
c2["W"] := [1, 5, 7, 35];
c2["shape"] := "V";
c2["over"] := 3;
c2["argument"] := "paper/classification-certificates.tex, Section 7 (shape V)";
certificates := [c, t, c2];
