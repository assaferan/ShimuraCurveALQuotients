// Subhyperelliptic cover models for X_0(21,34) -- the quotient by W = {1, 3, 14, 17, 42, 51, 238, 714} only
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
//
// Built from the Heegner heights of the branch points of X/W -> E = X/<W, w_2> = 357d1, by
// tools/heegner_model.m with the degree of the map doubled; see PROVENANCE.md. The quaternionic form
// here is old at 2, so E has conductor 357 rather than 714. The residual Atkin-Lehner involution
// acts by x -> -x.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,3,14,17,42,51,238,714]] := [* <3, P![ 49, 0, 220, 0, 262, 0, -612, 0, 81 ], P![]> *];
