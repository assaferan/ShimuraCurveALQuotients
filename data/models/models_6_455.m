// Subhyperelliptic cover models for X_0(6,455) -- the quotient by
// W = {1, 6, 7, 10, 13, 15, 42, 70, 78, 91, 105, 130, 195, 546, 910, 1365} only
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
//
// Built from the Heegner heights of the branch points of X/W -> E = X/<W, w_30> = 2730y1 (branch
// points of discriminants -10920 over Q(sqrt -26) and -1560 over Q(sqrt 10)), by
// tools/heegner_model.m; see PROVENANCE.md. The residual Atkin-Lehner involution acts by x -> -x.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,6,7,10,13,15,42,70,78,91,105,130,195,546,910,1365]] := [* <3, P![ 40, 0, -251, 0, 582, 0, -515, 0, -416 ], P![]> *];
