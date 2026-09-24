// Subhyperelliptic cover models for X_0(62,3)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,62]] := [* <2, P![ -9/4096, 0, -1/512, 0, -29/4096, 0, -1/2048 ], P![]> *];
models[[Integers()|1,6,31,186]] := [* <2, P![ 0, -155/65536, 29175/16384, -17625/32768, -625/16384, 3125/65536 ], P![]> *];
models[[Integers()|1,31]] := [*  *];
models[[Integers()|1,93]] := [*  *];
models[[Integers()|1]] := [*  *];
models[[Integers()|1,186]] := [* <3, P![ -31/65536, 0, 1167/16384, 0, -141/32768, 0, -1/16384, 0, 1/65536 ], P![]> *];
models[[Integers()|1,2,93,186]] := [* <1, P![ -31/65536, 5835/16384, -3525/32768, -125/16384, 625/65536 ], P![]> *];
models[[Integers()|1,6,62,93]] := [* <1, P![ -9/4096, -5/512, -725/4096, -125/2048 ], P![]> *];
models[[Integers()|1,3]] := [*  *];
models[[Integers()|1,2]] := [*  *];
models[[Integers()|1,3,62,186]] := [* <0, P![ 0, 5 ], P![]> *];
models[[Integers()|1,6]] := [*  *];
models[[Integers()|1,2,31,62]] := [* <1, P![ 0, -45/4096, -25/512, -3625/4096, -625/2048 ], P![]> *];
