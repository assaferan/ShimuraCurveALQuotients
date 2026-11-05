// Subhyperelliptic cover models for X_0(142,1)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
//
// The star curve has genus 1; the model is its minimal Weierstrass model, 142a1, the curve of the
// newform of level 142 with Atkin-Lehner signs -1 at 2 and 71 (a class with one curve), built by
// the genus-one star driver (tools/genus1star.m) from Borcherds forms with prescribed divisors
// and their values by Schofer's formula; data/models/PROVENANCE.md records the run.
// Guo-Yang, Compositio 153 (2017), Example 42. E(Q) is torsion-free, so the divisors of the
// Borcherds forms are principal and no assumption on their characters is used.
//
// Rational CM points of fundamental discriminant on the model, with Q = (1, 1) a generator of
// E(Q) = Z (the origin is the CM point of discriminant -8), the sign as in Guo-Yang's table:
//   -3: -Q = (1, -3)      -4: -2Q = (5, -11)     -8: O             -19: Q = (1, 1)
//   -20: 2Q = (5, 5)      -24: 3Q = (-3, 5)      -40: 4Q = (3, -3) -43: -3Q = (-3, -3)
//   -148: -4Q = (3, -1)   -232: -6Q = (13, 37)
// The two conjugate points of discriminant -187 sum to -2Q, those of discriminant -120 to 2Q.
// The table is Guo-Yang's, entry for entry. The positions of -4 and -8, which fix the branch loci
// of the three double covers, are confirmed independently of the Borcherds computation by
// tests/PrymBranchPositions.m: of the possible P_-4 - P_-8 only 2Q gives the Prym of X/<w_142>
// the Atkin-Lehner signs that Jacquet-Langlands requires.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,71,142]] := [* <1, P![ 15, -12, -1, 1 ], P![ 1, 1 ]> *];
