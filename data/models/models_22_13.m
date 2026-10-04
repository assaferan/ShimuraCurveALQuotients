// Subhyperelliptic cover models for X_0(22,13)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,11]] := [* <5, P![ -1936, 0, -286891/16, 0, 175571/16, 0, -405713/128, 0, 102971/256, 0, -58443/4096, 0, -1331/512 ], P![]> *];
models[[Integers()|1,26]] := [* <2, P![ -1331, 0, -2176, 0, 249856/121, 0, -16777216/14641 ], P![]> *];
models[[Integers()|1,2,13,26]] := [* <1, P![ -1, 963/256, -2723/512, 13739/4096, -3267/4096 ], P![]> *];
models[[Integers()|1,13]] := [*  *];
models[[Integers()|1,286]] := [* <2, P![ 121/86528, 0, -55/43264, 0, 19/40898, 0, 16/2474329 ], P![]> *];
models[[Integers()|1,2,11,22]] := [* <2, P![ -1, 3907/256, -3561/64, 187561/2048, -320771/4096, 2237411/65536, -395307/65536 ], P![]> *];
models[[Integers()|1,22]] := [*  *];
models[[Integers()|1,2,143,286]] := [* <1, P![ 1, -27/2, 497/16, -209/8, 121/16 ], P![]> *];
models[[Integers()|1,143]] := [*  *];

// Built 2026-10-03 as fibre products over the star line of committed double covers
// (FibreProductCovers.m); each verified against the trace formula at three primes, F_p and F_{p^2}.
// Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
models[[Integers()|1,2]] := [* <4, "CRV", [ Strings() | "-121/16*s^4 + 209/8*s^3*z - 497/16*s^2*z^2 + 27/2*s*z^3 - z^4 + y1^2", "395307/65536*s^6 - 2237411/65536*s^5*z + 320771/4096*s^4*z^2 - 187561/2048*s^3*z^3 + 3561/64*s^2*z^4 - 3907/256*s*z^5 + z^6 + y2^2" ]> *];
models[[Integers()|1,22,26,143]] := [* <1, P![ 0, 495616, -1368752, 1267112, -395307 ], P![]> *];
models[[Integers()|1,13,22,286]] := [* <1, P![ 0, -121/43264, 3025/86528, -35937/692224, 14641/692224 ], P![]> *];
models[[Integers()|1,11,26,286]] := [* <0, P![ 0, -121/16, 121/16 ], P![]> *];
models[[Integers()|1,11,13,143]] := [* <3, P![ 0, -16, -593, 1451, -1677, 851, -121, -88 ], P![ 0, 1, 0, 1 ]> *];
