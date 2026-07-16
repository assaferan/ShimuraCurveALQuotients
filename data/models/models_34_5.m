// Subhyperelliptic cover models for X_0(34,5)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,34]] := [*  *];
models[[Integers()|1,2]] := [*  *];
models[[Integers()|1,17]] := [*  *];
models[[Integers()|1,2,17,34]] := [* <2, P![ -1/6, -8/9, -23/16, -13/72, 173/144, 13/36, -11/36 ], P![]> *];
models[[Integers()|1]] := [*  *];
