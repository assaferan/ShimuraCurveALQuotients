// Subhyperelliptic cover models for X_0(202,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2]] := [* <4, P![ -1728, 18752, -89343, 248029, -1790021/4, 550989, -7522581/16, 4398465/16, -6748619/64, 766647/32, -156643/64 ], P![]> *];
models[[Integers()|1,101]] := [* <1, P![ -27, 131, -879/4, 313/2, -163/4 ], P![]> *];

// Built 2026-10-03 as fibre products over the star line of committed double covers
// (FibreProductCovers.m); each verified against the trace formula at three primes, F_p and F_{p^2}.
// Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
models[[Integers()|1]] := [* <8, "CRV", [ Strings() | "163/4*s^4 - 313/2*s^3*z + 879/4*s^2*z^2 - 131*s*z^3 + 27*z^4 + y1^2", "156643/64*s^10 - 766647/32*s^9*z + 6748619/64*s^8*z^2 - 4398465/16*s^7*z^3 + 7522581/16*s^6*z^4 - 550989*s^5*z^5 + 1790021/4*s^4*z^6 - 248029*s^3*z^7 + 89343*s^2*z^8 - 18752*s*z^9 + 1728*z^10 + y2^2" ]> *];
