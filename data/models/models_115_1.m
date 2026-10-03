// Subhyperelliptic cover models for X_0(115,1)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1]] := [*  *];
models[[Integers()|1,23]] := [* <1, P![ -8/25, 44/25, -96/25, 92/25, -27/25 ], P![]> *];
