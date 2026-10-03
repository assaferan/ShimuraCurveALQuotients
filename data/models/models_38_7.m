// Subhyperelliptic cover models for X_0(38,7)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,19,38]] := [* <2, P![ -5, 81/2, -2327/16, 2371/8, -5775/16, 249, -76 ], P![]> *];
models[[Integers()|1,38]] := [*  *];
models[[Integers()|1,19]] := [*  *];
models[[Integers()|1,133]] := [*  *];
models[[Integers()|1,2,133,266]] := [* <1, P![ 4, -12, 17, -8 ], P![]> *];
models[[Integers()|1,266]] := [*  *];

// Built 2026-10-03 as fibre products over the star line of committed double covers
// (FibreProductCovers.m); each verified against the trace formula at three primes, F_p and F_{p^2}.
// Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
models[[Integers()|1,2]] := [* <7, "CRV", [ Strings() | "76*s^6 - 249*s^5*z + 5775/16*s^4*z^2 - 2371/8*s^3*z^3 + 2327/16*s^2*z^4 - 81/2*s*z^5 + 5*z^6 + y1^2", "8*s^3*z - 17*s^2*z^2 + 12*s*z^3 - 4*z^4 + y2^2" ]> *];
