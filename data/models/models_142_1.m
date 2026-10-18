// Subhyperelliptic cover models for X_0(142,1)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
//
// The star curve has genus 1 (the first entry in this directory whose star is not rational), so
// there is no Hauptmodul; the model is the minimal Weierstrass model of the star, built exactly
// by the genus-one star driver (tools/genus1star.m): two coordinates of degree 2 from Borcherds
// forms with prescribed divisors on the rational CM points, their values by Schofer's formula,
// the bidegree-(2,2) relation 4x^2y^2 + 5xy - x - y + 1 = 0, and its minimal model. The curve is
// 142a1, the Jacobian of the newform of level 142 with Atkin-Lehner signs -1 at 2 and 71; the
// isogeny class has one curve. This is Guo-Yang, Compositio 153 (2017), Example 42.
//
// Rational CM points on the model, with Q = (1, 1) a generator of E(Q) = Z (the origin is the CM
// point of discriminant -8); Guo-Yang's table lists the same points with Q replaced by -Q:
//   -3: Q = (1, 1)        -4: 2Q = (5, 5)        -8: O             -19: -Q = (1, -3)
//   -20: -2Q = (5, -11)   -24: -3Q = (-3, -3)    -40: -4Q = (3, -1) -43: 3Q = (-3, 5)
//   -148: 4Q = (3, -3)    -232: 6Q = (13, -51)
// The two conjugate points of discriminant -187 sum to 2Q, those of discriminant -120 to -2Q.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,71,142]] := [* <1, P![ 15, -12, -1, 1 ], P![ 1, 1 ]> *];
