// Models for X_0(2418,1): the Atkin-Lehner star (genus two, from the differentials, tools/genus2model.m;
// a CANDIDATE, see PROVENANCE.md) and its double cover by the quotient by <w_2, w_13, w_93> (genus
// four, derived from the star: tools/models-2418_1/, PROVENANCE.md; hyperellipticity proved by the
// Borcherds certificate of tests/Hyperelliptic2418_1.m).
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,3,6,13,26,31,39,62,78,93,186,403,806,1209,2418]] :=
    [* <2, P![ 4, 30, -6, -5, -9, -3, 9 ], P![ 0, 1, 1 ]> *];
models[[Integers()|1,2,13,26,93,186,1209,2418]] :=
    [* <4, P![ 169, 0, 172, 0, 5190, 0, 5068, 0, 3481, 0, 2304 ], P![ 0 ]> *];
