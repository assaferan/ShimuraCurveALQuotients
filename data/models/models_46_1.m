// Subhyperelliptic cover models for X_0(46,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2]] := [* <0, P![ -27, 299/4, -207/4 ], P![]> *];
models[[Integers()|1]] := [* <1, P![ -4/25, 0, 9/20, 0, -8/25 ], P![]> *];
models[[Integers()|1,46]] := [* <0, P![ 1, -9/4, 5/4 ], P![]> *];
models[[Integers()|1,23]] := [* <1, P![ -243, 2439/2, -36531/16, 15111/8, -9315/16 ], P![]> *];
