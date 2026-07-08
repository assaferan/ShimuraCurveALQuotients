// Subhyperelliptic cover models for X_0(106,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2]] := [* <2, P![ -19/4096, 39/256, -3663/2048, 4131/512, -19683/4096, -17253/512, -19683/1024 ], P![]> *];
models[[Integers()|1]] := [*  *];
models[[Integers()|1,106]] := [* <1, P![ 1/16, -9/8, 45/16, 27, 81/4 ], P![]> *];
models[[Integers()|1,53]] := [* <1, P![ -19/256, 21/32, 63/128, -297/32, -2187/256 ], P![]> *];
