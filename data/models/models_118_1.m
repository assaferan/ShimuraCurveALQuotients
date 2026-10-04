// Subhyperelliptic cover models for X_0(118,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,118]] := [* <1, P![ 16, -82, 629/4, -267/2, 169/4 ], P![]> *];
models[[Integers()|1,2]] := [* <2, P![ -6912, 56096, -189920, 686641/2, -5591559/16, 1519313/8, -688675/16 ], P![]> *];
models[[Integers()|1,59]] := [* <1, P![ -432, 2156, -16075/4, 6627/2, -4075/4 ], P![]> *];

// Built 2026-10-03 as fibre products over the star line of committed double covers
// (FibreProductCovers.m); each verified against the trace formula at three primes, F_p and F_{p^2}.
// Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
models[[Integers()|1]] := [* <4, "CRV", [ Strings() | "688675/16*s^6 - 1519313/8*s^5*z + 5591559/16*s^4*z^2 - 686641/2*s^3*z^3 + 189920*s^2*z^4 - 56096*s*z^5 + 6912*z^6 + y1^2", "4075/4*s^4 - 6627/2*s^3*z + 16075/4*s^2*z^2 - 2156*s*z^3 + 432*z^4 + y2^2" ]> *];
