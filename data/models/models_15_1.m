// Subhyperelliptic cover models for X_0(15,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,3]] := [* <0, P![ -1/48, 41/144, -3/64 ], P![]> *];
models[[Integers()|1,5]] := [* <1, P![ 0, 3/32, -41/32, 27/128 ], P![]> *];
models[[Integers()|1]] := [* <1, P![ -1/48, 0, -41/72, 0, -3/16 ], P![]> *];
models[[Integers()|1,15]] := [* <0, P![ 0, -1/2 ], P![]> *];
