// Subhyperelliptic cover models for X_0(1326,1) -- the quotient by W = {1, 17, 78, 1326} only
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
//
// Built from the Heegner heights of the branch points of X/W -> E = X/<W, w_6> = 1326b2 (branch
// points of discriminants -24 and -408 over Q(sqrt 2)), by tools/heegner_model.m; see PROVENANCE.md.
// The Klein four group of Atkin-Lehner involutions acts by x -> -x (w_2, quotient 1326d1),
// x -> 1/(2x) (w_6, quotient 1326b2) and x -> -1/(2x) (w_3, quotient 1326c1).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,17,78,1326]] := [* <3, P![ 17, 0, -18, 0, -55, 0, -72, 0, 272 ], P![]> *];
