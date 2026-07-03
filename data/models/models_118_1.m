// Subhyperelliptic cover models for X_0(118,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,118]] := [* <1, P![ 16, -82, 629/4, -267/2, 169/4 ], P![]> *];
models[[Integers()|1,2]] := [* <2, P![ -6912, 56096, -189920, 686641/2, -5591559/16, 1519313/8, -688675/16 ], P![]> *];
models[[Integers()|1]] := [*  *];
models[[Integers()|1,59]] := [* <1, P![ -432, 2156, -16075/4, 6627/2, -4075/4 ], P![]> *];
