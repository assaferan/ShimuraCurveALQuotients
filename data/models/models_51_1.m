// Subhyperelliptic cover models for X_0(51,1)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,3]] := [* <1, P![ -2187/256, -5589/64, -27297/128, 6507/64, -2187/256 ], P![]> *];
models[[Integers()|1,17]] := [* <2, P![ 0, 6561/256, 16767/64, 81891/128, -19521/64, 6561/256 ], P![]> *];
models[[Integers()|1,51]] := [* <0, P![ 0, -3 ], P![]> *];
models[[Integers()|1]] := [* <3, P![ -2187/256, 0, 1863/64, 0, -3033/128, 0, -241/64, 0, -27/256 ], P![]> *];
