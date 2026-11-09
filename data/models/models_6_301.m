// Subhyperelliptic cover models for X_0(6,301) -- the quotient by W = {1, 3, 14, 42, 43, 129, 602, 1806} only
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
//
// Built from the Heegner heights of the branch points of X/W -> E = X/<W, w_903> = 1806l1 (branch
// points of discriminants -1204 over Q(sqrt -7) and -1032 over Q(sqrt 129)), by tools/heegner_model.m; see
// PROVENANCE.md. The residual Atkin-Lehner involution acts by x -> -x.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,3,14,42,43,129,602,1806]] := [* <3, P![ -63, 0, 1180, 0, 3334, 0, -2276, 0, 129 ], P![]> *];
