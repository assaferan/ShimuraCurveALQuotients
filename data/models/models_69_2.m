// Subhyperelliptic cover models for X_0(69,2)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,23,46]] := [* <2, P![ 0, 3/529, 29/2116, 29/1058, 9/2116, -11/529, -16/529 ], P![]> *];
models[[Integers()|1,3,46,138]] := [* <1, P![ -4/3, -8/9, 4/9, 16/9 ], P![]> *];
models[[Integers()|1,69]] := [*  *];
models[[Integers()|1,6,23,138]] := [* <2, P![ 0, 1/3, 53/36, 49/18, 89/36, 1 ], P![]> *];
models[[Integers()|1,2,69,138]] := [* <1, P![ 0, -1/9, -7/36, 1/18, 1/4 ], P![]> *];
models[[Integers()|1,46]] := [*  *];
models[[Integers()|1,3]] := [*  *];
models[[Integers()|1,2]] := [*  *];
models[[Integers()|1,6]] := [*  *];
models[[Integers()|1,2,3,6]] := [* <2, P![ -4/14283, -74/42849, -955/171396, -35/3174, -2335/171396, -419/42849, -16/4761 ], P![]> *];
models[[Integers()|1,138]] := [*  *];
models[[Integers()|1,3,23,69]] := [* <2, P![ 4/42849, 14/42849, 35/57132, 31/85698, -95/171396, -4/4761 ], P![]> *];
models[[Integers()|1,23]] := [*  *];
models[[Integers()|1,6,46,69]] := [* <1, P![ 0, -4/4761, -7/4761, -16/4761 ], P![]> *];
