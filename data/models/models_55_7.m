// Subhyperelliptic cover models for X_0(55,7) -- the quotient by W = {1, 35, 55, 77} only
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
//
// Built from the Heegner heights of the branch points of X/W -> E = X/<W, w_5> = 385a1 (branch
// points of discriminants -1540 over Q(sqrt -35) and -20 over Q(sqrt 5)), by tools/heegner_model.m;
// see PROVENANCE.md. The residual Atkin-Lehner involution acts by x -> -x.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,35,55,77]] := [* <3, P![ 5, 0, 68, 0, 198, 0, -2028, 0, -35 ], P![]> *];
