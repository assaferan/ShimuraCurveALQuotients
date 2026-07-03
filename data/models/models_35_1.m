// Subhyperelliptic cover models for X_0(35,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,35]] := [* <0, P![ 0, 25 ], P![]> *];
models[[Integers()|1,5]] := [* <2, P![ 0, -7/256, -345/64, -277/128, -25/64, -7/256 ], P![]> *];
models[[Integers()|1,7]] := [* <1, P![ -7/256, -345/64, -277/128, -25/64, -7/256 ], P![]> *];
models[[Integers()|1]] := [* <3, P![ -7/256, 0, -69/320, 0, -277/80000, 0, -1/40000, 0, -7/100000000 ], P![]> *];
