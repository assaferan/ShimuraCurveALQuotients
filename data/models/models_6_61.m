// Subhyperelliptic cover models for X_0(6,61)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,6,122,183]] := [* <3, P![ 0, -1/36864, 4589/21233664, -3109/5308416, 1781/3538944, 1511/5308416, -9283/21233664, -3/16384 ], P![]> *];
models[[Integers()|1,2,61,122]] := [* <1, P![ -1/4, 325/256, -161/128, -243/256 ], P![]> *];
models[[Integers()|1,3,61,183]] := [* <1, P![ 0, -1/36, 325/2304, -161/1152, -27/256 ], P![]> *];
models[[Integers()|1,6,61,366]] := [* <0, P![ 0, 1/36 ], P![]> *];
models[[Integers()|1,6]] := [* <5, P![ -1/4096, 0, 4589/65536, 0, -27981/4096, 0, 432783/2048, 0, 1101519/256, 0, -60905763/256, 0, -14348907/4 ], P![]> *];
models[[Integers()|1,2,3,6]] := [* <2, P![ -1/4096, 4589/2359296, -3109/589824, 1781/393216, 1511/589824, -9283/2359296, -27/16384 ], P![]> *];
models[[Integers()|1,2,183,366]] := [* <1, P![ 81/1024, -117/512, 153/1024, 9/64 ], P![]> *];
models[[Integers()|1,3,122,366]] := [* <1, P![ 0, 9/1024, -13/512, 17/1024, 1/64 ], P![]> *];
models[[Integers()|1,61]] := [* <2, P![ -1/4, 0, 2925/64, 0, -13041/8, 0, -177147/4 ], P![]> *];
models[[Integers()|1,366]] := [* <2, P![ 81/1024, 0, -1053/128, 0, 12393/64, 0, 6561 ], P![]> *];

// Built 2026-10-03 as fibre products over the star line of committed double covers
// (FibreProductCovers.m); each verified against the trace formula at three primes, F_p and F_{p^2}.
// Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
models[[Integers()|1,2]] := [* <4, "CRV", [ Strings() | "243/256*s^3*z + 161/128*s^2*z^2 - 325/256*s*z^3 + 1/4*z^4 + y1^2", "-9/64*s^3*z - 153/1024*s^2*z^2 + 117/512*s*z^3 - 81/1024*z^4 + y2^2" ]> *];
models[[Integers()|1,122]] := [* <5, "CRV", [ Strings() | "-1/64*s^4 - 17/1024*s^3*z + 13/512*s^2*z^2 - 9/1024*s*z^3 + y1^2", "243/256*s^3*z + 161/128*s^2*z^2 - 325/256*s*z^3 + 1/4*z^4 + y2^2" ]> *];
models[[Integers()|1,3]] := [* <4, "CRV", [ Strings() | "-1/64*s^4 - 17/1024*s^3*z + 13/512*s^2*z^2 - 9/1024*s*z^3 + y1^2", "27/16384*s^6 + 9283/2359296*s^5*z - 1511/589824*s^4*z^2 - 1781/393216*s^3*z^3 + 3109/589824*s^2*z^4 - 4589/2359296*s*z^5 + 1/4096*z^6 + y2^2" ]> *];
models[[Integers()|1,183]] := [* <5, "CRV", [ Strings() | "3/16384*s^7*z + 9283/21233664*s^6*z^2 - 1511/5308416*s^5*z^3 - 1781/3538944*s^4*z^4 + 3109/5308416*s^3*z^5 - 4589/21233664*s^2*z^6 + 1/36864*s*z^7 + y1^2", "27/256*s^4 + 161/1152*s^3*z - 325/2304*s^2*z^2 + 1/36*s*z^3 + y2^2" ]> *];
