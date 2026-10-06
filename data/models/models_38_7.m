// Subhyperelliptic cover models for X_0(38,7)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,133,266]] := [* <1, P![ 64, -192, 272, -128 ], P![]> *];
models[[Integers()|1,14,19,266]] := [* <1, P![ 80, -368, 724, -704, 256 ], P![]> *];
models[[Integers()|1,133]] := [*  *];
models[[Integers()|1,2,7,14]] := [* <4, P![ -5/4, 111/8, -4611/64, 14637/64, -124747/256, 92233/128, -188783/256, 16065/32, -821/4, 38 ], P![]> *];
models[[Integers()|1,266]] := [* <2, P![ 19, 0, 1/4, 0, 1/2, 0, 1/4 ], P![]> *];
models[[Integers()|1,7,38,266]] := [* <0, P![ 5, -8 ], P![]> *];
models[[Integers()|1,14]] := [*  *];
models[[Integers()|1,19]] := [*  *];
models[[Integers()|1,2]] := [*  *];
models[[Integers()|1,2,19,38]] := [* <2, P![ -5, 81/2, -2327/16, 2371/8, -5775/16, 249, -76 ], P![]> *];
models[[Integers()|1,38]] := [* <4, P![ -76, 0, 9/2, 0, -225/512, 0, 141/16384, 0, -1763/4194304, 0, -49/67108864 ], P![]> *];
models[[Integers()|1,7]] := [*  *];
models[[Integers()|1,14,38,133]] := [* <2, P![ -4, 2, -5, -2, 3, -4, 2 ], P![ 0, 1, 1 ]> *];
