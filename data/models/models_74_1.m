// Subhyperelliptic cover models for X_0(74,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2]] := [* <2, P![ 0, -1369/65536, -2079/32768, 473/32768, -41/8192, 47/65536, -1/32768 ], P![]> *];
models[[Integers()|1,37]] := [* <2, P![ -1369/65536, -2079/32768, 473/32768, -41/8192, 47/65536, -1/32768 ], P![]> *];
models[[Integers()|1]] := [* <4, P![ -1369/65536, 0, -2079/32768, 0, 473/32768, 0, -41/8192, 0, 47/65536, 0, -1/32768 ], P![]> *];
models[[Integers()|1,74]] := [* <0, P![ 0, 1 ], P![]> *];
