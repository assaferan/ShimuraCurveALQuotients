// Subhyperelliptic cover models for X_0(714,5) -- the quotient by
// W = {1, 3, 10, 14, 17, 30, 35, 42, 51, 105, 170, 238, 510, 595, 714, 1785} only
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
//
// Built from the Heegner heights of the branch points of X/W -> E = X/<W, w_2> = 3570x2, by
// tools/heegner_model.m; see PROVENANCE.md. The branch locus consists of one place of degree 2 of
// discriminant -14280 and the rational points of discriminants -840 and -340. The residual
// Atkin-Lehner involution acts by x -> -x.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,3,10,14,17,30,35,42,51,105,170,238,510,595,714,1785]] := [* <3, P![ 16, 0, -155, 0, 598, 0, 2061, 0, 9720 ], P![]> *];
