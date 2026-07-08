// Subhyperelliptic cover models for X_0(34,7)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,7,14]] := [* <2, P![ -5/16, 9/8, 35/64, -455/64, 625/64, -173/64, -27/16 ], P![]> *];
models[[Integers()|1,2,17,34]] := [* <1, P![ -5, 28, -249/4, 259/4, -27 ], P![]> *];
models[[Integers()|1,14]] := [*  *];
models[[Integers()|1,2]] := [*  *];
models[[Integers()|1,34]] := [*  *];
models[[Integers()|1,7]] := [*  *];
models[[Integers()|1,2,119,238]] := [* <2, P![ 1/4, -3/4, -11/16, 39/8, -95/16, 3/2, 1 ], P![]> *];
models[[Integers()|1,238]] := [*  *];
models[[Integers()|1,7,17,119]] := [* <2, P![ 5, 1, 2, -5, -8, -5, -3 ], P![ 0, 1, 1 ]> *];
models[[Integers()|1,17]] := [* <3, P![ -187/81, 323/81, -1199/324, 613/162, -127/36, 173/81, -89/81, 44/81, -4/27 ], P![]> *];
models[[Integers()|1,14,17,238]] := [* <0, P![ -1/4, -1/4, 1 ], P![]> *];
models[[Integers()|1,119]] := [*  *];
