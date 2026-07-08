// Subhyperelliptic cover models for X_0(38,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2]] := [* <1, P![ 0, -171/256, -3321/128, -43011/256, -6561/16 ], P![]> *];
models[[Integers()|1,38]] := [* <0, P![ 0, 9 ], P![]> *];
models[[Integers()|1,19]] := [* <1, P![ -19/256, -369/128, -4779/256, -729/16 ], P![]> *];
models[[Integers()|1]] := [* <2, P![ -19/256, 0, -41/128, 0, -59/256, 0, -1/16 ], P![]> *];
