// Subhyperelliptic cover models for X_0(10,3)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,10]] := [*  *];   // degenerate CRV removed 2026-09-06: stored its parent conic twice
models[[Integers()|1,2,5,10]] := [* <0, P![ -2, 43/20, -7/20 ], P![]> *];
models[[Integers()|1,2,3,6]] := [* <0, P![ 1, -15/8, 7/8 ], P![]> *];
models[[Integers()|1,30]] := [*  *];
models[[Integers()|1,2,15,30]] := [* <0, P![ 1, -6/5, 1/5 ], P![]> *];
models[[Integers()|1,15]] := [*  *];   // degenerate CRV removed 2026-09-06: stored its parent conic twice
models[[Integers()|1]] := [*  *];   // degenerate CRV removed 2026-09-06: stored its parent conic twice
models[[Integers()|1,3]] := [*  *];
models[[Integers()|1,2]] := [* <0, P![ 8/845, 16/845, -20/169 ], P![]>, <0, P![ 1, -11/8, 1 ], P![]>, <0, P![ 27/2450, 0, 27/1960 ], P![]> *];
models[[Integers()|1,5]] := [*  *];
models[[Integers()|1,6]] := [*  *];   // degenerate CRV removed 2026-09-06: stored its parent conic twice