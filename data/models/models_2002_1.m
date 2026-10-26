// Models for X_0(2002,1): the Atkin-Lehner star (genus two, from the differentials, tools/genus2model.m)
// and its double cover by the quotient by <w_2, w_13, w_77> (genus four, from the star: see
// PROVENANCE.md and paper/star-2002-proof.tex, Section 3).
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,7,11,13,14,22,26,77,91,143,154,182,286,1001,2002]] :=
    [* <2, P![ 4, 64, 196, -456, 180, 4, 16 ], P![ 0, 1, 1 ]> *];
models[[Integers()|1,2,13,26,77,154,1001,2002]] :=
    [* <4, P![ 65536, 0, -16384, 0, 185856, 0, -279360, 0, 98289, 0, -1001 ], P![ 0 ]> *];
