// Subhyperelliptic cover models for X_0(77,1)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,11]] := [* <1, P![ -11, -10, -31, -20, -16 ], P![]> *];
models[[Integers()|1,7]] := [* <3, P![ 0, 11, 10, 53, 149/4, 151/2, 129/4, 27, -4 ], P![]> *];
models[[Integers()|1]] := [*  *];
models[[Integers()|1,77]] := [* <1, P![ 0, -4, 0, -8, 1 ], P![]> *];
