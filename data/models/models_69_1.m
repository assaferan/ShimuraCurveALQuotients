// Subhyperelliptic cover models for X_0(69,1)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,3]] := [* <1, P![ -4/6561, 1/8748, 35/139968, -37/839808, -1/27648 ], P![]> *];
models[[Integers()|1,69]] := [* <0, P![ 9, -9 ], P![]> *];
models[[Integers()|1]] := [* <3, P![ -1/3072, 0, -7/186624, 0, -37/30233088, 0, 317/1224440064, 0, -1/181398528 ], P![]> *];
models[[Integers()|1,23]] := [* <2, P![ -4/81, 19/324, 19/1728, -247/10368, 53/82944, 3/1024 ], P![]> *];
