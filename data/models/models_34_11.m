// Subhyperelliptic cover models for X_0(34,11)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,11]] := [*  *];
models[[Integers()|1,187]] := [*  *];
models[[Integers()|1,17]] := [*  *];
models[[Integers()|1,34]] := [*  *];
models[[Integers()|1,2]] := [*  *];
models[[Integers()|1,22]] := [*  *];
models[[Integers()|1,2,187,374]] := [* <1, P![ 1, -1, -11/4, -1/2, 1/4 ], P![]> *];
models[[Integers()|1,11,34,374]] := [* <1, P![ 5, -13, 9/4, 7/2, -3/4 ], P![]> *];
models[[Integers()|1,17,22,374]] := [* <0, P![ 5/121, 2/121, -3/121 ], P![]> *];
models[[Integers()|1,374]] := [* <2, P![ -12/14641, -152/14641, -740/14641, -1880/14641, -2700/14641, -192/1331, -64/1331 ], P![]> *];
