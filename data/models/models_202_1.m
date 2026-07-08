// Subhyperelliptic cover models for X_0(202,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2]] := [* <4, P![ -1728, 18752, -89343, 248029, -1790021/4, 550989, -7522581/16, 4398465/16, -6748619/64, 766647/32, -156643/64 ], P![]> *];
models[[Integers()|1,101]] := [* <1, P![ -27, 131, -879/4, 313/2, -163/4 ], P![]> *];
models[[Integers()|1]] := [*  *];
