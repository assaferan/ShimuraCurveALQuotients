// Subhyperelliptic cover models for X_0(55,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,5]] := [* <1, P![ -1/45375, -2/45375, -21/15125, -138/15125, -243/3025 ], P![]> *];
models[[Integers()|1,11]] := [* <2, P![ -1/5671875, 4/5671875, -32/1890625, -42/1890625, -1332/1890625, 216/378125, -2187/75625 ], P![]> *];
models[[Integers()|1]] := [* <3, P![ -1/5808, -1/8712, 1/5808, 1/2178, -1/5808, -1/2178, 1/5808, 1/8712, -1/5808 ], P![]> *];
models[[Integers()|1,55]] := [* <0, P![ 1/5, -6/5, 9 ], P![]> *];
