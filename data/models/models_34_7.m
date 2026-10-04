// Subhyperelliptic cover models for X_0(34,7)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,7,14]] := [* <2, P![ -5/16, 9/8, 35/64, -455/64, 625/64, -173/64, -27/16 ], P![]> *];
models[[Integers()|1,2,17,34]] := [* <1, P![ -5, 28, -249/4, 259/4, -27 ], P![]> *];
models[[Integers()|1,34]] := [*  *];
models[[Integers()|1,7]] := [*  *];
models[[Integers()|1,2,119,238]] := [* <2, P![ 1/4, -3/4, -11/16, 39/8, -95/16, 3/2, 1 ], P![]> *];
models[[Integers()|1,7,17,119]] := [* <2, P![ 5, 1, 2, -5, -8, -5, -3 ], P![ 0, 1, 1 ]> *];
models[[Integers()|1,17]] := [* <3, P![ -187/81, 323/81, -1199/324, 613/162, -127/36, 173/81, -89/81, 44/81, -4/27 ], P![]> *];
models[[Integers()|1,14,17,238]] := [* <0, P![ -1/4, -1/4, 1 ], P![]> *];
models[[Integers()|1,119]] := [*  *];

// Built 2026-10-03 as fibre products over the star line of committed double covers
// (FibreProductCovers.m); each verified against the trace formula at three primes, F_p and F_{p^2}.
// Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
models[[Integers()|1,2]] := [* <5, "CRV", [ Strings() | "27*s^4 - 259/4*s^3*z + 249/4*s^2*z^2 - 28*s*z^3 + 5*z^4 + y1^2", "27/16*s^6 + 173/64*s^5*z - 625/64*s^4*z^2 + 455/64*s^3*z^3 - 35/64*s^2*z^4 - 9/8*s*z^5 + 5/16*z^6 + y2^2" ]> *];
models[[Integers()|1,14]] := [* <3, "CRV", [ Strings() | "-s^2 + 1/4*s*z + 1/4*z^2 + y1^2", "27/16*s^6 + 173/64*s^5*z - 625/64*s^4*z^2 + 455/64*s^3*z^3 - 35/64*s^2*z^4 - 9/8*s*z^5 + 5/16*z^6 + y2^2" ]> *];
models[[Integers()|1,238]] := [* <3, "CRV", [ Strings() | "-s^2 + 1/4*s*z + 1/4*z^2 + y1^2", "-s^6 - 3/2*s^5*z + 95/16*s^4*z^2 - 39/8*s^3*z^3 + 11/16*s^2*z^4 + 3/4*s*z^5 - 1/4*z^6 + y2^2" ]> *];
models[[Integers()|1,14,34,119]] := [* <1, P![ 5/16, -23/16, 137/64, -25/32, -27/64 ], P![]> *];
models[[Integers()|1,7,34,238]] := [* <1, P![ -1/16, 1/4, -21/64, 7/64, 1/16 ], P![]> *];
