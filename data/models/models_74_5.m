// Subhyperelliptic cover models for X_0(74,5)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,5,37,185]] := [* <3, P![ -4864/25, -1792/25, 11584/25, -20224/25, -20352/25, 9216/25, -24768/25, 6144/25, -4096/25 ], P![]> *];
models[[Integers()|1,2,5,10]] := [* <5, P![ -311296/625, -1048576/625, -303104/625, 516096/625, -3963904/625, -6002688/625, -2466816/625, -1060864/125, -3982336/625, 43008/625, -2374656/625, 360448/625, -65536/125 ], P![]> *];
models[[Integers()|1,10]] := [*  *];
models[[Integers()|1,185]] := [*  *];
models[[Integers()|1,74]] := [*  *];
models[[Integers()|1,5,74,370]] := [* <1, P![ 1600, 4800, 3600, 800, 2000 ], P![]> *];
models[[Integers()|1,10,74,185]] := [* <1, P![ -76/25, -104/25, 4/5, -32/25 ], P![]> *];
models[[Integers()|1,37]] := [*  *];
models[[Integers()|1,2]] := [*  *];
models[[Integers()|1,2,185,370]] := [* <2, P![ 64, -64, -48, 288, -112, 128 ], P![]> *];
models[[Integers()|1,5]] := [*  *];
models[[Integers()|1,370]] := [*  *];
