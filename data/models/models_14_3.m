// Subhyperelliptic cover models for X_0(14,3)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2]] := [* <1, "CRV", [ Strings() | "y^2 + 128*s^2 - 104*s*z + 11*z^2", "x^2 + 128*s^2 - 104*s*z + 11*z^2" ]> *];
models[[Integers()|1,7]] := [*  *];
models[[Integers()|1]] := [*  *];
models[[Integers()|1,2,7,14]] := [* <0, P![ -11, 104, -128 ], P![]> *];
models[[Integers()|1,14]] := [*  *];
