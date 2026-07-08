// Subhyperelliptic cover models for X_0(26,5)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,10]] := [*  *];
models[[Integers()|1,13]] := [*  *];
models[[Integers()|1,5,13,65]] := [* <1, P![ -76, 164, -107, 60, -44 ], P![]> *];
models[[Integers()|1,2,65,130]] := [* <1, P![ 1, -3, 9/4, -1, 1 ], P![]> *];
models[[Integers()|1]] := [*  *];
models[[Integers()|1,10,13,130]] := [* <2, P![ 15, 1, -3, -5, -1 ], P![ 0, 1, 0, 1 ]> *];
models[[Integers()|1,2]] := [*  *];
models[[Integers()|1,2,5,10]] := [* <2, P![ -125, 45, 44, 39, -15, -8, -5 ], P![ 0, 1, 0, 1 ]> *];
models[[Integers()|1,5]] := [* <3, P![ -3596, 488, 9748, 5560, -6064, -8488, -4044, -888, -76 ], P![]> *];
models[[Integers()|1,2,13,26]] := [* <1, P![ -19/4, 15, -63/4, 15, -11 ], P![]> *];
models[[Integers()|1,10,26,65]] := [* <0, P![ -19/4, 15, -11 ], P![]> *];
models[[Integers()|1,5,26,130]] := [* <0, P![ 4, 0, 4 ], P![]> *];
models[[Integers()|1,26]] := [* <1, P![ -155/4, 40, 63/2, -8, -19/4 ], P![]>, <1, P![ 61/121, 454/121, 1299/121, 1690/121, 845/121 ], P![]> *];
models[[Integers()|1,130]] := [* <3, P![ 145, 130, -199, -346, -156, 14, 33, 10, 1 ], P![]> *];
models[[Integers()|1,65]] := [*  *];
