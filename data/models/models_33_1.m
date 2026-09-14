// Subhyperelliptic cover models for X_0(33,1)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,3]] := [* <0, P![ -16/1089, 31/1089, -27/1936 ], P![]> *];
models[[Integers()|1,11]] := [* <1, P![ -16/121, 47/121, -739/1936, 243/1936 ], P![]> *];
models[[Integers()|1]] := [* <1, P![ -1/5808, 0, -5/8712, 0, -27/1936 ], P![]> *];
models[[Integers()|1,33]] := [* <0, P![ 1, -1 ], P![]> *];
