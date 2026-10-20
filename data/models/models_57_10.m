// Subhyperelliptic cover models for X_0(57,10) -- the quotient by W = {1, 5, 6, 19, 30, 95, 114, 570} only
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
//
// Built from the Heegner heights of the branch points of X/W -> E = X/<W, w_2> = 285a1, by
// tools/heegner_model.m; see PROVENANCE.md. The quaternionic form is old at 2, so E has conductor
// 285 rather than 570, and the degree entering the formula is the one of the combination of the two
// degeneracy maps. The residual Atkin-Lehner involution acts by x -> -x.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,5,6,19,30,95,114,570]] := [* <3, P![ -15, 0, 68, 0, -2, 0, 500, 0, 729 ], P![]> *];
