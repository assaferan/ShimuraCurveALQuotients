// Subhyperelliptic cover models for X_0(122,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,122]] := [* <1, P![ 16, -48, 52, -32 ], P![]> *];
models[[Integers()|1,2]] := [* <3, P![ -11264, 73728, -211456, 360960, -408512, 313216, -161728, 51712, -8192 ], P![]> *];
models[[Integers()|1,61]] := [* <2, P![ -704, 2496, -3440, 2720, -1200, 256 ], P![]> *];

// Built 2026-10-03 as fibre products over the star line of committed double covers
// (FibreProductCovers.m); each verified against the trace formula at three primes, F_p and F_{p^2}.
// Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
models[[Integers()|1]] := [* <6, "CRV", [ Strings() | "32*s^3*z - 52*s^2*z^2 + 48*s*z^3 - 16*z^4 + y1^2", "8192*s^8 - 51712*s^7*z + 161728*s^6*z^2 - 313216*s^5*z^3 + 408512*s^4*z^4 - 360960*s^3*z^5 + 211456*s^2*z^6 - 73728*s*z^7 + 11264*z^8 + y2^2" ]> *];
