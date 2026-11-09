// Models for X_0(2030,1): the Atkin-Lehner star (genus two, from the differentials, tools/genus2model.m;
// a CANDIDATE, see PROVENANCE.md) and its double cover by the quotient by <w_10, w_14, w_58> (genus
// four, derived from the star: tools/models-2030_1/, PROVENANCE.md; hyperellipticity proved by the
// Borcherds certificate of tests/Hyperelliptic2030_1.m).
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,5,7,10,14,29,35,58,70,145,203,290,406,1015,2030]] :=
    [* <2, P![ 39, 149, 159, 9, 13, 29, 1 ], P![ 0, 1, 1 ]> *];
models[[Integers()|1,10,14,35,58,145,203,2030]] :=
    [* <4, P![ 625, 0, -1100, 0, -186, 0, 2796, 0, 81, 0, -1960 ], P![ 0 ]> *];
