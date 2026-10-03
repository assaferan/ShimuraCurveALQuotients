// Subhyperelliptic cover models for X_0(161,1)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,23]] := [* <3, P![ -8/49, 76/49, -319/49, 830/49, -1661/49, 2770/49, -3403/49, 356/7, -16 ], P![]> *];
