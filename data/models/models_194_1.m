// Subhyperelliptic cover models for X_0(194,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,97]] := [* <4, P![ -1, 3, -55/4, 47/2, -243/4, 86, -114, 152, -108, 96, -64 ], P![]> *];
models[[Integers()|1,2]] := [* <5, P![ -13, 35, -747/4, 621/2, -3883/4, 1345, -2353, 3240, -3076, 3856, -2608, 1664, -1280 ], P![]> *];
models[[Integers()|1]] := [* <9, P![ -19, -92, -286, -592, -921, -1016, -872, 460, 1545, 1752, 34, -1752, 1545, -460, -872, 1016, -921, 592, -286, 92, -19 ], P![]> *];
models[[Integers()|1,194]] := [* <0, P![ 13, 4, 20 ], P![]> *];
