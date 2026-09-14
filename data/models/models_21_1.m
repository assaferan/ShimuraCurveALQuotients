// Subhyperelliptic cover models for X_0(21,1)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,3]] := [* <1, P![ 0, 734832, 7663248, 18475776, -26873856 ], P![]> *];
models[[Integers()|1,7]] := [* <0, P![ -7, -80, -256 ], P![]> *];
models[[Integers()|1]] := [* <1, P![ -7, 0, 94, 0, -343 ], P![]> *];
models[[Integers()|1,21]] := [* <0, P![ 0, -16/9, 16/9 ], P![]> *];
