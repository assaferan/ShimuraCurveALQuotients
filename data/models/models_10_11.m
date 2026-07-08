// Subhyperelliptic cover models for X_0(10,11)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,10]] := [* <2, P![ -11, 0, -538, 0, -131, 0, -8 ], P![]> *];
models[[Integers()|1,110]] := [* <0, P![ 8, 0, 1 ], P![]>, <0, P![ -2, -2 ], P![]>, <0, P![ -8, 0, 1 ], P![]> *];
models[[Integers()|1,55]] := [*  *];
models[[Integers()|1,2,11,22]] := [* <1, P![ 0, 5/4, 11/2, 61/4, -2 ], P![]> *];
models[[Integers()|1]] := [* <5, P![ -688, 816, 3052, 1740, -4935, -11814, -13625, -10296, -5505, -2110, -563, -96, -8 ], P![]> *];
models[[Integers()|1,2]] := [*  *];
models[[Integers()|1,2,5,10]] := [* <1, P![ -1/10, -171/400, -233/200, 5/16, -1/50 ], P![]> *];
models[[Integers()|1,5]] := [* <3, P![ -1/10, 0, -171/400, 0, -233/200, 0, 5/16, 0, -1/50 ], P![]> *];
models[[Integers()|1,5,22,110]] := [* <0, P![ 0, 1 ], P![]> *];
models[[Integers()|1,5,11,55]] := [* <2, P![ 0, -5/32, -171/256, -233/128, 125/256, -1/32 ], P![]> *];
models[[Integers()|1,10,22,55]] := [* <1, P![ 5, 22, 61, -8 ], P![]> *];
models[[Integers()|1,2,55,110]] := [* <0, P![ 0, -8/25, 1/25 ], P![]> *];
models[[Integers()|1,10,11,110]] := [* <0, P![ -8, 1 ], P![]> *];
models[[Integers()|1,11]] := [* <3, P![ -22, 0, -4315/4, 0, -793/2, 0, -195/4, 0, -2 ], P![]> *];
models[[Integers()|1,22]] := [* <2, P![ 5, 0, 22, 0, 61, 0, -8 ], P![]> *];
