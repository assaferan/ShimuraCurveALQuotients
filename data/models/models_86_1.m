// Subhyperelliptic cover models for X_0(86,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2]] := [* <2, P![ 0, -43/65536, -185/16384, -753/32768, -189/16384, 245/65536, -1/4096 ], P![]> *];
models[[Integers()|1,86]] := [* <0, P![ 0, 1 ], P![]> *];
models[[Integers()|1,43]] := [* <2, P![ -43/65536, -185/16384, -753/32768, -189/16384, 245/65536, -1/4096 ], P![]> *];
models[[Integers()|1]] := [* <4, P![ -43/65536, 0, -185/16384, 0, -753/32768, 0, -189/16384, 0, 245/65536, 0, -1/4096 ], P![]> *];
