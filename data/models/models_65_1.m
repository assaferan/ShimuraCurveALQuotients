// Subhyperelliptic cover models for X_0(65,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,65]] := [* <1, P![ 625, 1250, -3125, -1250, 625 ], P![]> *];
models[[Integers()|1,5]] := [* <2, P![ -125000, -62500, 1078125, -656250, -1078125, 718750, -109375 ], P![]> *];
models[[Integers()|1]] := [*  *];
models[[Integers()|1,13]] := [* <2, P![ -125000, 437500, -171875, -718750, 484375, 281250, -109375 ], P![]> *];
