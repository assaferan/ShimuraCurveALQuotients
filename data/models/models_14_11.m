// Subhyperelliptic cover models for X_0(14,11)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,7,14]] := [* <2, P![ -80/7, 128/7, 296/49, -1248/49, -79/49, 340/49, -32/49 ], P![]> *];
models[[Integers()|1,2,11,22]] := [* <1, P![ -4/7, 4/7, 33/49, -44/49, -32/49 ], P![]> *];
models[[Integers()|1,14]] := [*  *];
models[[Integers()|1,154]] := [* <2, P![ 0, -32/49, -332/49, -1304/49, -2412/49, -2112/49, -704/49 ], P![]> *];
models[[Integers()|1]] := [*  *];
models[[Integers()|1,2]] := [* <3, P![ -512/49, -10592/49, -13220/7, -9080, -1285468/49, -2284864/49, -2443712/49, -1441792/49, -360448/49 ], P![]> *];
models[[Integers()|1,7]] := [*  *];
models[[Integers()|1,2,77,154]] := [* <0, P![ 20, -12, 1 ], P![]> *];
models[[Integers()|1,7,11,77]] := [* <1, P![ -20/7, -4, -59/49, 8/49 ], P![]> *];
models[[Integers()|1,7,22,154]] := [* <1, P![ 4/49, -12/49, 13/49, -4/49 ], P![]> *];
models[[Integers()|1,11,14,154]] := [* <1, P![ 20/49, -52/49, 45/49, -4/49 ], P![]> *];
models[[Integers()|1,77]] := [* <2, P![ 0, 64/49, 788/49, 3496/49, 6868/49, 6144/49, 2048/49 ], P![]> *];
models[[Integers()|1,14,22,77]] := [* <1, P![ -4/7, -4/7, 5/49, 8/49 ], P![]> *];
models[[Integers()|1,11]] := [*  *];
models[[Integers()|1,22]] := [*  *];
