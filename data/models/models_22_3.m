// Subhyperelliptic cover models for X_0(22,3)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
// Non-hyperelliptic covers are stored as <genus, "CRV", [defining-poly strings]>.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[ 1 ]] := [* <3, P![ -251/64, -377/24, -8179/288, -715/24, -11281/576, -589/72, -607/288, -11/36, -11/576 ], P![]> *];
models[[ 1, 2 ]] := [* <2, P![ -11/64, 0, -35/32, 0, -131/64, 0, -27/16 ], P![]> *];
models[[ 1, 6 ]] := [* <2, P![ 11/4096, 0, -49/36864, 0, -23/36864, 0, -3/4096 ], P![]> *];
models[[ 1, 66 ]] := [* <0, P![ 1, 0, 4 ], P![]>, <0, P![ 1/20, 0, -1/4 ], P![]>, <0, P![ -1/4, 0, 1/4 ], P![]> *];
models[[ 1, 3 ]] := [*  *];   // degenerate CRV removed 2026-09-06: stored its parent conic twice
models[[ 1, 33 ]] := [* <1, P![ -11/576, 0, -13/288, 0, -3/64 ], P![]> *];
models[[ 1, 11 ]] := [* <1, P![ -11/1024, 0, -25/4608, 0, -3/1024 ], P![]> *];
models[[ 1, 2, 11, 22 ]] := [* <1, P![ -11/64, 35/128, -131/1024, 27/1024 ], P![]> *];
models[[ 1, 6, 11, 66 ]] := [* <0, P![ 1, -1 ], P![]> *];
models[[ 1, 3, 11, 33 ]] := [* <0, P![ -11/576, 13/1152, -3/1024 ], P![]> *];
models[[ 1, 2, 3, 6 ]] := [* <1, P![ -81, -12, 1, 1 ], P![ 1, 1 ]> *];
models[[ 1, 2, 33, 66 ]] := [* <0, P![ 0, -1/4 ], P![]> *];
models[[ 1, 3, 22, 66 ]] := [* <0, P![ 0, -1/20, 1/20 ], P![]> *];
models[[ 1, 6, 22, 33 ]] := [* <1, P![ 0, 11/2304, -13/4608, 3/4096 ], P![]> *];
