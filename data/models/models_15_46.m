// Subhyperelliptic cover models for X_0(15,46) -- the quotient by W = {1, 5, 6, 23, 30, 115, 138, 690} only
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
//
// Built from the Heegner heights of the branch points of X/W -> E = X/<W, w_2> = 345f1, by
// tools/heegner_model.m; see PROVENANCE.md. The quaternionic form is old at 2, so E has conductor
// 345 rather than 690, and the degree entering the formula is the one of the combination of the two
// degeneracy maps. The residual Atkin-Lehner involution acts by x -> -x.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,5,6,23,30,115,138,690]] := [* <3, P![ 81, 0, 380, 0, 358, 0, -900, 0, 81 ], P![]> *];
