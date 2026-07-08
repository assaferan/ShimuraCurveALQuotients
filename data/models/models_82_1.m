// Subhyperelliptic cover models for X_0(82,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,41]] := [* <0, P![ -108, 164, -67 ], P![]> *];
models[[Integers()|1,2]] := [* <2, P![ -27648, 138752, -294128, 336544, -218968, 76760, -11323 ], P![]> *];
models[[Integers()|1]] := [* <3, "CRV", [ Strings() | "y^2 - 169/1024*s^4 + 183/256*s^3*z - 301/256*s^2*z^2 + 7/8*s*z^3 - 1/4*z^4", "x^2 + 67*s^2 - 164*s*z + 108*z^2" ]> *];
models[[Integers()|1,82]] := [* <1, P![ 1/4, -7/8, 301/256, -183/256, 169/1024 ], P![]> *];
