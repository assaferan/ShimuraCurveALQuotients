// Subhyperelliptic cover models for X_0(15,38) -- the quotient by W = {1, 2, 15, 19, 30, 38, 285, 570} only
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
//
// Built from the Heegner heights of the branch points of X/W -> E = X/<W, w_2> = 570e4, by
// tools/heegner_model.m; see PROVENANCE.md. All four branch points are rational, of discriminants
// -228, -12, -760 and -40; the place of discriminant -12 is the one whose order is not maximal, and
// its position was searched rather than predicted. The residual Atkin-Lehner involution acts by
// x -> -x.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,15,19,30,38,285,570]] := [* <3, P![ 361, 0, -1660, 0, 2062, 0, 180, 0, 81 ], P![]> *];
