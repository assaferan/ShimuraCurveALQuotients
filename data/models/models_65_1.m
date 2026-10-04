// Subhyperelliptic cover models for X_0(65,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,65]] := [* <1, P![ 625, 1250, -3125, -1250, 625 ], P![]> *];
models[[Integers()|1,5]] := [* <2, P![ -125000, -62500, 1078125, -656250, -1078125, 718750, -109375 ], P![]> *];
models[[Integers()|1,13]] := [* <2, P![ -125000, 437500, -171875, -718750, 484375, 281250, -109375 ], P![]> *];

// Built 2026-10-03 as fibre products over the star line of committed double covers
// (FibreProductCovers.m); each verified against the trace formula at three primes, F_p and F_{p^2}.
// Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
models[[Integers()|1]] := [* <5, "CRV", [ Strings() | "109375*s^6 - 281250*s^5*z - 484375*s^4*z^2 + 718750*s^3*z^3 + 171875*s^2*z^4 - 437500*s*z^5 + 125000*z^6 + y1^2", "109375*s^6 - 718750*s^5*z + 1078125*s^4*z^2 + 656250*s^3*z^3 - 1078125*s^2*z^4 + 62500*s*z^5 + 125000*z^6 + y2^2" ]> *];
