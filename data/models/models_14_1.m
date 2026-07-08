// Subhyperelliptic cover models for X_0(14,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2]] := [* <0, P![ -49, 931/16, -539/16 ], P![]> *];
models[[Integers()|1,7]] := [* <1, P![ 0, 2401, -84035/16, 36015/8, -26411/16 ], P![]> *];
models[[Integers()|1]] := [* <1, P![ -49, 0, 637/16, 0, -49/2 ], P![]> *];
models[[Integers()|1,14]] := [* <0, P![ 0, -1/4, 1/4 ], P![]> *];
