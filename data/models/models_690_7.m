// Subhyperelliptic cover models for X_0(690,7) -- the quotient by
// W = {1, 6, 10, 14, 15, 21, 23, 35, 138, 210, 230, 322, 345, 483, 805, 4830} only
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
//
// Built from the Heegner heights of the branch points of X/W -> E = X/<W, w_2> = 4830bj1, by
// tools/heegner_model.m; see PROVENANCE.md. All four branch points are rational, of discriminants
// -3, -280, -420 and -115. The residual Atkin-Lehner involution acts by x -> -x.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,6,10,14,15,21,23,35,138,210,230,322,345,483,805,4830]] := [* <3, P![ 144, 0, -520, 0, -767, 0, -130, 0, 9 ], P![]> *];
