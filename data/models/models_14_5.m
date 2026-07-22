// Subhyperelliptic cover models for X_0(14,5)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,10]] := [*  *];
models[[Integers()|1,2,7,14]] := [* <0, P![ -11, 216, -1024 ], P![]> *];
models[[Integers()|1,14]] := [* <0, P![ 5/17, 150/289, 5/17 ], P![]>, <0, P![ 5, 0, -16 ], P![]>, <0, P![ 5/16, 0, -1/16 ], P![]> *];
models[[Integers()|1,70]] := [* <1, P![ -7/16, 0, 11/8, 0, 1/16 ], P![]> *];
models[[Integers()|1]] := [* <3, P![ -14375/4096, -22875/1024, -124175/2048, -23565/256, -352037/4096, -3249/64, -38075/2048, -3963/1024, -1439/4096 ], P![]> *];
models[[Integers()|1,2]] := [* <1, P![ -8/289, 0, 9/289, 0, -8/289 ], P![]>, <1, P![ -11, 86, -267, 404, -256 ], P![]> *];
models[[Integers()|1,2,5,10]] := [* <1, P![ -100/289, 0, 57/289, 0, -8/289 ], P![]> *];
models[[Integers()|1,5]] := [*  *];
models[[Integers()|1,35]] := [* <1, P![ -7/4096, 0, -181/2048, 0, 1/4096 ], P![]> *];
models[[Integers()|1,7]] := [*  *];
models[[Integers()|1,2,35,70]] := [* <0, P![ 1, -12, 4 ], P![]> *];
models[[Integers()|1,5,14,70]] := [* <0, P![ 1, -8 ], P![]> *];
models[[Integers()|1,10,14,35]] := [* <0, P![ -11, 128 ], P![]> *];
