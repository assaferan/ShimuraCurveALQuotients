// Subhyperelliptic cover models for X_0(158,1)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
//
// The star curve has genus 1; the model is its minimal Weierstrass model, 158a1, the curve of the
// newform of level 158 with Atkin-Lehner signs -1 at 2 and 79 (a class with one curve), built by
// the genus-one star driver (tools/genus1star.m) from Borcherds forms with prescribed divisors
// and their values by Schofer's formula. We know of no published equation. E(Q) is torsion-free,
// so the divisors of the Borcherds forms are principal and no assumption on their characters is
// used.
//
// Rational CM points of fundamental discriminant on the model, with Q = (-1, -4) a generator of
// E(Q) = Z (the origin is the CM point of discriminant -4):
//   -4: O                 -8: 6Q = (11, 28)      -11: 4Q = (7, -20)    -19: 2Q = (3, 0)
//   -20: 5Q = (-3, 2)     -40: 3Q = (1, 0)       -52: Q = (-1, -4)     -67: 10Q = (139, -1704)
//   -88: 8Q = (11/4, -23/8)                      -163: -4Q = (7, 12)
// The two conjugate points of discriminant -115 sum to 12Q, those of -123 to 4Q, of -184 to 10Q
// and of -51 to 6Q.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,79,158]] := [* <1, P![ 9, -9, -1, 1 ], P![ 1, 1 ]> *];
