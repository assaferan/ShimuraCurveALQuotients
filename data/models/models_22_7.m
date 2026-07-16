// Subhyperelliptic cover models for X_0(22,7)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,7,14]] := [* <1, P![ -20, -134, -1095/4, -164, -27 ], P![]> *];
models[[Integers()|1,2,11,22]] := [* <2, P![ -5, -97/2, -2723/16, -4075/16, -9399/64, -61/2, -27/16 ], P![]> *];
models[[Integers()|1,14]] := [* <1, P![ -3/1024, 5/512, 39/1024, -59/512, -291/1024 ], P![]> *];
models[[Integers()|1,154]] := [* <1, P![ 5, 22, 21, 0, 16 ], P![]>, <1, P![ 1/4096, -1/2048, -9/4096, 21/2048, 89/4096 ], P![]> *];
models[[Integers()|1]] := [*  *];
models[[Integers()|1,2]] := [* <3, P![ -20, -188, -679, -1232, -1422, -1460, -1063, -416, -432 ], P![]> *];
models[[Integers()|1,7]] := [*  *];
models[[Integers()|1,2,77,154]] := [* <0, P![ 1/4, 3/4, 1/16 ], P![]> *];
models[[Integers()|1,7,11,77]] := [* <1, P![ 12, 116, 209, 82, -27 ], P![]> *];
models[[Integers()|1,7,22,154]] := [* <1, P![ 44, 44, -75, -64, 16 ], P![]> *];
models[[Integers()|1,11,14,154]] := [* <0, P![ 5, 16, 4 ], P![]> *];
models[[Integers()|1,77]] := [* <1, P![ -4, -20, -31, -26, -27 ], P![]> *];
models[[Integers()|1,14,22,77]] := [* <0, P![ -4, -14, -27/4 ], P![]> *];
models[[Integers()|1,22]] := [*  *];
models[[Integers()|1,11]] := [*  *];
