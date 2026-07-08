// Subhyperelliptic cover models for X_0(34,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,34]] := [* <0, P![ 1/4, 0, 1 ], P![]> *];
models[[Integers()|1,2]] := [* <1, P![ -11/16, 17/8, -71/16, 17/2, -27/4 ], P![]> *];
models[[Integers()|1,17]] := [* <0, P![ -11/16, 17/8, -27/16 ], P![]> *];
models[[Integers()|1]] := [* <1, P![ -3/16, 1/8, -65/16, -27/8, -11/16 ], P![]> *];
