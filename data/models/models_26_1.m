// Subhyperelliptic cover models for X_0(26,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2]] := [* <1, P![ 0, -169/16, -3/2, 19/16, -1/8 ], P![]> *];
models[[Integers()|1]] := [* <2, P![ -169, 0, -24, 0, 19, 0, -2 ], P![]> *];
models[[Integers()|1,13]] := [* <1, P![ -169, -24, 19, -2 ], P![]> *];
models[[Integers()|1,26]] := [* <0, P![ 0, 1 ], P![]> *];
