// Subhyperelliptic cover models for X_0(206,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2]] := [* <4, P![ -1/512, 13/4096, 21/2048, 331/4096, 55/1024, -733/4096, -3323/2048, -19883/4096, -3605/512, -1139/256, -1 ], P![]> *];
models[[Integers()|1,206]] := [* <0, P![ 0, 1 ], P![]> *];
models[[Integers()|1]] := [* <9, P![ -1/512, 0, 13/4096, 0, 21/2048, 0, 331/4096, 0, 55/1024, 0, -733/4096, 0, -3323/2048, 0, -19883/4096, 0, -3605/512, 0, -1139/256, 0, -1 ], P![]> *];
models[[Integers()|1,103]] := [* <5, P![ 0, 1024, -4556, 7210, -4971, 1661, -184, -56, 82, -11, 3, 2 ], P![ 0, 0, 1, 1, 1, 1 ]> *];
