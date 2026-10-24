// Model for the Atkin-Lehner star X_0(2002,1)* -- from the differentials (tools/genus2model.m)
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,7,11,13,14,22,26,77,91,143,154,182,286,1001,2002]] :=
    [* <2, P![ 4, 64, 196, -456, 180, 4, 16 ], P![ 0, 1, 1 ]> *];
