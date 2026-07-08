// Subhyperelliptic cover models for X_0(94,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2]] := [* <1, P![ -27, 313/4, -1341/16, 2491/64, -423/64 ], P![]> *];
models[[Integers()|1,94]] := [* <0, P![ 1, -9/4, 5/4 ], P![]> *];
models[[Integers()|1]] := [* <3, P![ -1/2500, 0, 381/40000, 0, -117/1250, 0, 276/625, 0, -512/625 ], P![]> *];
models[[Integers()|1,47]] := [* <2, P![ -1/3, 139/81, -29/8, 1735/432, -5659/2304, 8131/10368, -235/2304 ], P![]> *];
