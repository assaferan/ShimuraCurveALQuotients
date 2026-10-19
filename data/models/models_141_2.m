// Subhyperelliptic cover models for X_0(141,2) -- the quotient by W = {1, 3, 94, 282} only
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
//
// Built from the Heegner heights of the branch points of X/W -> E = X/<W, w_2> = 141a1, by
// tools/heegner_model.m with the degree of the map doubled; see PROVENANCE.md. The quaternionic form
// here is old at 2, so E has conductor 141 rather than 282. The branch locus is the rational points
// of discriminants -4 and -24 and a place of degree 2 of discriminant -564 over Q(sqrt -47). The
// residual Atkin-Lehner involution acts by x -> -x.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,3,94,282]] := [* <3, P![ -47, 0, -260, 0, -26, 0, 252, 0, 81 ], P![]> *];
