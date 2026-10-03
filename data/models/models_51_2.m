// Subhyperelliptic cover models for X_0(51,2)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,3,17,51]] := [* <1, P![ 0, 2, -47/16, 3/8, 9/16 ], P![]> *];
models[[Integers()|1,2,51,102]] := [* <1, P![ 0, -24, 45/4, 27/4 ], P![]> *];
models[[Integers()|1,2,17,34]] := [* <2, P![ 0, 9, -63/4, 189/16, -567/256, -729/256 ], P![]> *];
models[[Integers()|1,34]] := [* <3, P![ -867/256, 0, -1413/64, 0, -9801/128, 0, -9477/64, 0, -19683/256 ], P![]> *];
models[[Integers()|1,6]] := [* <4, P![ 289/32, 0, -20475/256, 0, -24057/64, 0, -204849/128, 0, -124659/64, 0, -177147/256 ], P![]> *];
models[[Integers()|1,2,3,6]] := [* <2, P![ -1536, 3408, -2844, 567, 14013/16, -2673/8, -2187/16 ], P![]> *];
models[[Integers()|1,6,34,51]] := [* <0, P![ -1/3, 1/3 ], P![]> *];
models[[Integers()|1,51]] := [* <2, P![ -6, 0, 225/4, 0, 567/2, 0, 729/4 ], P![]> *];
models[[Integers()|1,3,34,102]] := [* <1, P![ 0, -3, 9/4, -27/16, -243/256 ], P![]> *];
models[[Integers()|1,6,17,102]] := [* <2, P![ 32, -39, 81/4, 135/16, -2511/256, -729/256 ], P![]> *];

// Built 2026-10-03 as fibre products over the star line of committed double covers
// (FibreProductCovers.m); each verified against the trace formula at three primes, F_p and F_{p^2}.
// Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
models[[Integers()|1]] := [* <9, "CRV", [ Strings() | "2187/16*s^6 + 2673/8*s^5*z - 14013/16*s^4*z^2 - 567*s^3*z^3 + 2844*s^2*z^4 - 3408*s*z^5 + 1536*z^6 + y1^2", "729/256*s^5*z + 2511/256*s^4*z^2 - 135/16*s^3*z^3 - 81/4*s^2*z^4 + 39*s*z^5 - 32*z^6 + y2^2", "-27/4*s^3*z - 45/4*s^2*z^2 + 24*s*z^3 + y3^2" ]> *];
models[[Integers()|1,2]] := [* <5, "CRV", [ Strings() | "729/256*s^5*z + 567/256*s^4*z^2 - 189/16*s^3*z^3 + 63/4*s^2*z^4 - 9*s*z^5 + y1^2", "-27/4*s^3*z - 45/4*s^2*z^2 + 24*s*z^3 + y2^2" ]> *];
models[[Integers()|1,102]] := [* <4, "CRV", [ Strings() | "729/256*s^5*z + 2511/256*s^4*z^2 - 135/16*s^3*z^3 - 81/4*s^2*z^4 + 39*s*z^5 - 32*z^6 + y1^2", "-27/4*s^3*z - 45/4*s^2*z^2 + 24*s*z^3 + y2^2" ]> *];
models[[Integers()|1,3]] := [* <4, "CRV", [ Strings() | "2187/16*s^6 + 2673/8*s^5*z - 14013/16*s^4*z^2 - 567*s^3*z^3 + 2844*s^2*z^4 - 3408*s*z^5 + 1536*z^6 + y1^2", "-9/16*s^4 - 3/8*s^3*z + 47/16*s^2*z^2 - 2*s*z^3 + y2^2" ]> *];
models[[Integers()|1,17]] := [* <5, "CRV", [ Strings() | "729/256*s^5*z + 567/256*s^4*z^2 - 189/16*s^3*z^3 + 63/4*s^2*z^4 - 9*s*z^5 + y1^2", "729/256*s^5*z + 2511/256*s^4*z^2 - 135/16*s^3*z^3 - 81/4*s^2*z^4 + 39*s*z^5 - 32*z^6 + y2^2" ]> *];
