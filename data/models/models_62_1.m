// Subhyperelliptic cover models for X_0(62,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2]] := [* <1, P![ -1/2, 43/32, -45/32, 99/128, -1/4 ], P![]> *];
models[[Integers()|1,31]] := [* <2, P![ 0, 1/4, -43/64, 45/64, -99/256, 1/8 ], P![]> *];
models[[Integers()|1,62]] := [* <0, P![ 0, -1/8 ], P![]> *];
models[[Integers()|1]] := [* <3, P![ -1/2, 0, -43/4, 0, -90, 0, -396, 0, -1024 ], P![]> *];
