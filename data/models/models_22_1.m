// Subhyperelliptic cover models for X_0(22,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,22]] := [* <0, P![ 0, -11/16 ], P![]> *];
models[[Integers()|1,2]] := [* <0, P![ 0, 121/16, -121/16 ], P![]> *];
models[[Integers()|1,11]] := [* <0, P![ -11, 11 ], P![]> *];
models[[Integers()|1]] := [* <0, P![ -11/16, 0, -11/16 ], P![]>, <0, P![ -11, 0, -16 ], P![]>, <0, P![ -11/16, 0, -1/16 ], P![]> *];
