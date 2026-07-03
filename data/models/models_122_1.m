// Subhyperelliptic cover models for X_0(122,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,122]] := [* <1, P![ 16, -48, 52, -32 ], P![]> *];
models[[Integers()|1,2]] := [* <3, P![ -11264, 73728, -211456, 360960, -408512, 313216, -161728, 51712, -8192 ], P![]> *];
models[[Integers()|1]] := [*  *];
models[[Integers()|1,61]] := [* <2, P![ -704, 2496, -3440, 2720, -1200, 256 ], P![]> *];
