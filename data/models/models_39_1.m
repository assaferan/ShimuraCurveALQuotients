// Subhyperelliptic cover models for X_0(39,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,3]] := [* <2, P![ -5/3, 8/3, 4/9, 14/9, -32/9, -4/3, -7/9 ], P![]> *];
models[[Integers()|1,39]] := [* <0, P![ 5, 2, 1 ], P![]> *];
models[[Integers()|1]] := [* <3, P![ -7/9, -10/3, -7/3, 16/3, 37/9, -16/3, -7/3, 10/3, -7/9 ], P![]> *];
models[[Integers()|1,13]] := [* <1, P![ -1/3, 2/3, -1/9, 2/9, -7/9 ], P![]> *];
