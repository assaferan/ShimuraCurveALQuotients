// Subhyperelliptic cover models for X_0(6,385) -- the quotient by
// W = {1, 5, 6, 11, 14, 21, 30, 55, 66, 70, 105, 154, 231, 330, 770, 1155} only
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
//
// Built from the Heegner heights of the branch points of X/W -> E = X/<W, w_2> = 2310q1, by
// tools/heegner_model.m; see PROVENANCE.md. The branch locus consists of the rational points of
// discriminants -40 and -1540 and one place of degree 2 of discriminant -9240. The residual
// Atkin-Lehner involution acts by x -> -x.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,5,6,11,14,21,30,55,66,70,105,154,231,330,770,1155]] := [* <3, P![ 3465, 0, 60092, 0, -9482, 0, -580, 0, 9 ], P![]> *];
