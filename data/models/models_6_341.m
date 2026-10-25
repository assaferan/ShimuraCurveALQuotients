// Models for X_0(6,341): the Atkin-Lehner star (genus two, from the differentials, tools/genus2model.m;
// a CANDIDATE, see PROVENANCE.md) and its double cover by the quotient by <w_3, w_11, w_62> (genus
// four, derived from the star: tools/models-6_341/, PROVENANCE.md; hyperellipticity of the cover is
// proved GIVEN the star, by the Borcherds certificate of tests/Hyperelliptic6_341.m and the star's
// two-torsion).
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,3,6,11,22,31,33,62,66,93,186,341,682,1023,2046]] :=
    [* <2, P![ 9, -33, 9, 29, -12, 6 ], P![ 0, 1, 1 ]> *];
models[[Integers()|1,3,11,33,62,186,682,2046]] :=
    [* <4, P![ 81, 0, 468, 0, 1126, 0, 1420, 0, -47, 0, -744 ], P![ 0 ]> *];
