// Subhyperelliptic cover models for X_0(77,1)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,11]] := [* <1, P![ -11, -10, -31, -20, -16 ], P![]> *];
models[[Integers()|1,7]] := [* <3, P![ 0, 11, 10, 53, 149/4, 151/2, 129/4, 27, -4 ], P![]> *];
models[[Integers()|1,77]] := [* <1, P![ 0, -4, 0, -8, 1 ], P![]> *];

// Built 2026-10-03 as fibre products over the star line of committed double covers
// (FibreProductCovers.m); each verified against the trace formula at three primes, F_p and F_{p^2}.
// Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
models[[Integers()|1]] := [* <5, "CRV", [ Strings() | "16*s^4 + 20*s^3*z + 31*s^2*z^2 + 10*s*z^3 + 11*z^4 + y1^2", "4*s^8 - 27*s^7*z - 129/4*s^6*z^2 - 151/2*s^5*z^3 - 149/4*s^4*z^4 - 53*s^3*z^5 - 10*s^2*z^6 - 11*s*z^7 + y2^2" ]> *];
