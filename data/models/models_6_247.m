// Subhyperelliptic cover models for X_0(6,247) -- the quotient by W = {1, 2, 19, 38, 39, 78, 741, 1482} only
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
//
// Built from the Heegner heights of the branch points of X/W -> E = X/<W, w_3> = 1482i1, by
// tools/heegner_model.m; see PROVENANCE.md. The branch locus consists of two rational points, of
// discriminants -52 and -3, and one place of degree 2 of discriminant -456 over Q(sqrt 57). The
// residual Atkin-Lehner involution acts by x -> -x.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,19,38,39,78,741,1482]] := [* <3, P![ 513, 0, 2540, 0, 1654, 0, -4788, 0, 81 ], P![]> *];
