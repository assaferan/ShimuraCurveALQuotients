// Subhyperelliptic cover models for X_0(6,185) -- the quotient by W = {1, 2, 15, 30, 111, 185, 222, 370} only
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
//
// Built from the Heegner heights of the branch points of X/W -> E = X/<W, w_3> = 1110m1, by
// tools/heegner_model.m; see PROVENANCE.md. The branch locus consists of one rational point of
// discriminant -555, one place of degree 2 of discriminant -4440 over Q(sqrt -111), and one rational
// point of discriminant -40. The residual Atkin-Lehner involution acts by x -> -x.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,15,30,111,185,222,370]] := [* <3, P![ -999, 0, -5620, 0, -122, 0, 4716, 0, 2025 ], P![]> *];
