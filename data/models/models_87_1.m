// Subhyperelliptic cover models for X_0(87,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,29]] := [* <3, P![ -5, 14, 23, -81, -36, 93, 70, 18, 3 ], P![ 1, 0, 1, 1 ]> *];
models[[Integers()|1,3]] := [* <2, P![ -129140163/3444736, 1190959281/1722368, -15635525661/3444736, 10581521751/861184, -34231709133/3444736, -8451506223/1722368, -10460353203/3444736 ], P![]> *];
models[[Integers()|1,87]] := [* <0, P![ 0, -27 ], P![]> *];
models[[Integers()|1]] := [* <5, P![ -129140163/3444736, 0, -44109603/1722368, 0, -21447909/3444736, 0, -537597/861184, 0, -64413/3444736, 0, 589/1722368, 0, -27/3444736 ], P![]> *];
