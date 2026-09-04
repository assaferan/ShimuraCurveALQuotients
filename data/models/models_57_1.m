// Subhyperelliptic cover models for X_0(57,1)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,3]] := [* <2, P![ -7/9, 16/3, -110/9, 104/9, -95/9, 20, -16 ], P![]> *];
models[[Integers()|1]] := [* <3, "CRV", [ Strings() | "y^2 - 9*s^4 - 2*s^2*z^2 + 4*s*z^3 - z^4", "x^2 + 144*s^2 - 180*s*z + 63*z^2" ]> *];
models[[Integers()|1,19]] := [* <0, P![ -63, 180, -144 ], P![]> *];
models[[Integers()|1,57]] := [* <1, P![ 1, -4, 2, 0, 9 ], P![]> *];
