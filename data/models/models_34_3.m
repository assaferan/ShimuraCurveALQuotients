// Subhyperelliptic cover models for X_0(34,3)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[ 1, 102 ]] := [* <2, P![ -3/16, 0, 1/8, 0, 13/16, 0, 1/4 ], P![]> *];
models[[ 1, 51 ]] := [* <2, P![ 288/13845841, -864/13845841, 23760/13845841, -46080/13845841, 347760/13845841, -324864/13845841, -144/226981 ], P![]> *];
models[[ 1, 17 ]] := [* <1, P![ -1/816, 0, -49/20808, 0, -3/4624 ], P![]> *];
models[[ 1, 2, 17, 34 ]] := [* <1, P![ 0, -36, 1305/16, -729/16 ], P![]> *];
models[[ 1, 6, 17, 102 ]] := [* <0, P![ -3, 3 ], P![]> *];
models[[ 1, 3, 17, 51 ]] := [* <0, P![ 0, 4/867, -27/4624 ], P![]> *];
models[[ 1, 2, 3, 6 ]] := [* <1, P![ -1/6, 283/384, -4967/4096, 3603/4096, -243/1024 ], P![]> *];
models[[ 1, 2, 51, 102 ]] := [* <1, P![ 0, 6, -207/16, 27/4 ], P![]> *];
models[[ 1, 3, 34, 102 ]] := [* <1, P![ 0, -2, 1, 1 ], P![ 0, 1 ]> *];
models[[ 1, 6, 34, 51 ]] := [* <1, P![ 1550, -256, 0, 1 ], P![ 1, 1 ]> *];

// Built 2026-10-03 as fibre products over the star line of committed double covers
// (FibreProductCovers.m); each verified against the trace formula at three primes, F_p and F_{p^2}.
// Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
models[[Integers()|1]] := [* <5, "CRV", [ Strings() | "243/1024*s^4 - 3603/4096*s^3*z + 4967/4096*s^2*z^2 - 283/384*s*z^3 + 1/6*z^4 + y1^2", "-3*s*z + 3*z^2 + y2^2", "-27/4*s^3*z + 207/16*s^2*z^2 - 6*s*z^3 + y3^2" ]> *];
models[[Integers()|1,2]] := [* <3, "CRV", [ Strings() | "729/16*s^3*z - 1305/16*s^2*z^2 + 36*s*z^3 + y1^2", "-27/4*s^3*z + 207/16*s^2*z^2 - 6*s*z^3 + y2^2" ]> *];
models[[Integers()|1,6]] := [* <2, "CRV", [ Strings() | "243/1024*s^4 - 3603/4096*s^3*z + 4967/4096*s^2*z^2 - 283/384*s*z^3 + 1/6*z^4 + y1^2", "-3*s*z + 3*z^2 + y2^2" ]> *];
models[[Integers()|1,3]] := [* <2, "CRV", [ Strings() | "243/1024*s^4 - 3603/4096*s^3*z + 4967/4096*s^2*z^2 - 283/384*s*z^3 + 1/6*z^4 + y1^2", "27/4624*s^2 - 4/867*s*z + y2^2" ]> *];
