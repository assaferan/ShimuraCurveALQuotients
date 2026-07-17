// Subhyperelliptic cover models for X_0(38,7)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2]] := [*  *];
models[[Integers()|1,2,19,38]] := [* <2, P![ -5, 81/2, -2327/16, 2371/8, -5775/16, 249, -76 ], P![]> *];
models[[Integers()|1,38]] := [*  *];
models[[Integers()|1,19]] := [*  *];
models[[Integers()|1,133]] := [*  *];
models[[Integers()|1,2,133,266]] := [* <1, P![ 4, -12, 17, -8 ], P![]> *];
models[[Integers()|1,266]] := [*  *];
