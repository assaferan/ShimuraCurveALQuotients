// Subhyperelliptic cover models for X_0(58,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,29]] := [* <0, P![ -1/27, 1/27 ], P![]> *];
models[[Integers()|1,2]] := [* <1, P![ -4/19683, -5/6561, -1519/1679616, -841/10077696 ], P![]> *];
models[[Integers()|1]] := [* <2, P![ -1/512, 0, -39/512, 0, -431/512, 0, -841/512 ], P![]> *];
models[[Integers()|1,58]] := [* <1, P![ 4/531441, 11/531441, 239/45349632, -8273/272097792, -841/272097792 ], P![]> *];
