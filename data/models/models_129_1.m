// Subhyperelliptic cover models for X_0(129,1)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,3]] := [* <4, P![ -5/27, 82/27, -1691/81, 6386/81, -14692/81, 21746/81, -21295/81, 13694/81, -5591/81, 148/9, -16/9 ], P![]> *];
models[[Integers()|1,129]] := [* <1, P![ 5/9, -32/9, 62/9, -4, 1 ], P![]> *];
models[[Integers()|1,43]] := [* <2, P![ -1/3, 10/3, -109/9, 62/3, -175/9, 28/3, -16/9 ], P![]> *];
models[[Integers()|1]] := [*  *];
