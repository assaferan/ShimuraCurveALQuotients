// Subhyperelliptic cover models for X_0(22,5)* -- Guo-Yang /  AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually  0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[ 1, 2, 5, 10 ]] := [* <1, P![ -1, 4755/1024, -8267/1024, 797/128,  -115/64 ], P![]> *];
models[[ 1, 2, 55, 110 ]] := [* <0, P![ 1, -9/4, 5/4 ], P![]> *];
models[[ 1, 2, 11, 22 ]] := [* <1, P![ -4096, 20044, -36799, 30040, -9200 ],  P![]> *];
