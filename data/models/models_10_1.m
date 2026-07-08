// Subhyperelliptic cover models for X_0(10,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2]] := [* <0, P![ -1/2, 29/4, -27/4 ], P![]> *];
models[[Integers()|1,5]] := [* <0, P![ 0, 2/25, -27/25 ], P![]> *];
models[[Integers()|1,10]] := [* <0, P![ 0, -1/100, 1/100 ], P![]> *];
models[[Integers()|1]] := [* <0, P![ -1/2916, 0, -1/1458 ], P![]>, <0, P![ -1/2, 0, -25/4 ], P![]>, <0, P![ -1/392, -1/392, -13/784 ], P![]> *];
