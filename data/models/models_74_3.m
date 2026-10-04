// Subhyperelliptic cover models for X_0(74,3)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,111]] := [* <7, "CRV", [ Strings() | "-1/8*s^3*z + 11/256*s^2*z^2 - 1/128*s*z^3 + 3/256*z^4 + y1^2", "-256*s^6 + 896*s^5*z - 1424*s^4*z^2 + 1120*s^3*z^3 - 528*s^2*z^4 + 128*s*z^5 + y2^2" ]> *];
models[[Integers()|1,222]] := [* <4, P![ 32, 0, 33, 0, 35/2, 0, 89/16, 0, 7/8, 0, 1/16 ], P![]> *];
models[[Integers()|1,3]] := [* <7, P![ -3/2048, 0, -115/65536, 0, -91/65536, 0, -513/524288, 0, -543/1048576, 0, -3379/16777216, 0, -437/8388608, 0, -123/16777216, 0, -1/2097152 ], P![]> *];
models[[Integers()|1,6,37,222]] := [* <2, P![ 32, -132, 280, -356, 224, -64 ], P![]> *];
models[[Integers()|1,2,37,74]] := [* <1, P![ 0, 3/256, -1/128, 11/256, -1/8 ], P![]> *];
models[[Integers()|1,6,74,111]] := [* <1, P![ -3/256, 1/128, -11/256, 1/8 ], P![]> *];
models[[Integers()|1,2,3,6]] := [* <3, P![ -3/2048, 115/16384, -91/4096, 513/8192, -543/4096, 3379/16384, -437/2048, 123/1024, -1/32 ], P![]> *];
models[[Integers()|1,3,74,222]] := [* <0, P![ 0, -4 ], P![]> *];
models[[Integers()|1,2,111,222]] := [* <2, P![ 0, -128, 528, -1120, 1424, -896, 256 ], P![]> *];
models[[Integers()|1,3,37,111]] := [* <4, P![ 20, 24, -34, -126, -49, 61, 79, 77, 25, 6 ], P![ 0, 0, 1, 0, 1 ]> *];
models[[Integers()|1,74]] := [* <2, P![ -3/256, 0, -1/512, 0, -11/4096, 0, -1/512 ], P![]> *];

// Built 2026-10-03 as fibre products over the star line of committed double covers
// (FibreProductCovers.m); each verified against the trace formula at three primes, F_p and F_{p^2}.
// Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
models[[Integers()|1,2]] := [* <6, "CRV", [ Strings() | "1/32*s^8 - 123/1024*s^7*z + 437/2048*s^6*z^2 - 3379/16384*s^5*z^3 + 543/4096*s^4*z^4 - 513/8192*s^3*z^5 + 91/4096*s^2*z^6 - 115/16384*s*z^7 + 3/2048*z^8 + y1^2", "-256*s^6 + 896*s^5*z - 1424*s^4*z^2 + 1120*s^3*z^3 - 528*s^2*z^4 + 128*s*z^5 + y2^2" ]> *];
models[[Integers()|1,6]] := [* <6, "CRV", [ Strings() | "1/32*s^8 - 123/1024*s^7*z + 437/2048*s^6*z^2 - 3379/16384*s^5*z^3 + 543/4096*s^4*z^4 - 513/8192*s^3*z^5 + 91/4096*s^2*z^6 - 115/16384*s*z^7 + 3/2048*z^8 + y1^2", "64*s^5*z - 224*s^4*z^2 + 356*s^3*z^3 - 280*s^2*z^4 + 132*s*z^5 - 32*z^6 + y2^2" ]> *];
models[[Integers()|1,37]] := [* <7, "CRV", [ Strings() | "64*s^5*z - 224*s^4*z^2 + 356*s^3*z^3 - 280*s^2*z^4 + 132*s*z^5 - 32*z^6 + y1^2", "1/8*s^4 - 11/256*s^3*z + 1/128*s^2*z^2 - 3/256*s*z^3 + y2^2" ]> *];
