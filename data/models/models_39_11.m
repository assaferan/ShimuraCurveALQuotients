// Subhyperelliptic cover models for X_0(39,11) -- the quotient by W = {1, 11, 39, 429} only
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
//
// Built from the Heegner heights of the branch points of X/W -> E = X/<W, w_3> = 429b1 (branch
// points of discriminants -132 over Q(sqrt -11) and -52 over Q(sqrt 13)), by tools/heegner_model.m;
// see PROVENANCE.md. The residual Atkin-Lehner involution acts by x -> -x.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,11,39,429]] := [* <3, P![ 13, 0, -124, 0, 214, 0, 1044, 0, -891 ], P![]> *];
