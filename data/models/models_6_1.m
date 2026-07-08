// Subhyperelliptic cover models for X_0(6,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,3]] := [* <0, P![ 0, 1/34992, -1/34992 ], P![]> *];
models[[Integers()|1,2]] := [* <0, P![ -1/16, 1/16 ], P![]> *];
models[[Integers()|1]] := [* <0, P![ -1/2187, 0, -16/2187 ], P![]>, <0, P![ -1/16, 0, -2187/16 ], P![]>, <0, P![ -1/2704, -1/1352, -1/208 ], P![]> *];
models[[Integers()|1,6]] := [* <0, P![ 0, -1/2187 ], P![]> *];
