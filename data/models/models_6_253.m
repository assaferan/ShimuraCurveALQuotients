// Subhyperelliptic cover models for X_0(6,253) -- the quotient by W = {1, 2, 23, 33, 46, 66, 759, 1518} only
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
//
// Built from the Heegner heights of the branch points of X/W -> E = X/<W, w_3> = 1518p1, by
// tools/heegner_model.m; see PROVENANCE.md. The branch locus consists of one place of degree 2 of
// discriminant -552 over Q(sqrt -23) and two rational points, of discriminants -88 and -1012. The
// residual Atkin-Lehner involution acts by x -> -x.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,23,33,46,66,759,1518]] := [* <3, P![ -1863, 0, -3060, 0, 3526, 0, 1388, 0, 9 ], P![]> *];
