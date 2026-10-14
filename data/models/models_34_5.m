// Subhyperelliptic cover models for X_0(34,5)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
// Fibre-product entries (FibreProductCovers.m, 2026-10-14): built over the star line from the double covers of
// this file by tools/writefp.m with the coordinate check of the three covers above each pair of factors;
// each checked against the trace formula at the first three good primes, F_p and F_{p^2}, and again at the
// 4th and 5th.  Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,10]] := [* <4, "CRV", [ Strings() | "11/16*s^4 - 1/8*s^3*z - 37/16*s^2*z^2 - 2*s*z^3 - 1/2*z^4 + y1^2", "-1/4*s^4 + 3/4*s^3*z - 1/16*s^2*z^2 - 7/8*s*z^3 + 3/16*z^4 + y2^2" ]> *];
models[[Integers()|1,5,34,170]] := [* <0, P![ -3, -4, 4 ], P![]> *];
models[[Integers()|1,10,34,85]] := [* <1, P![ 1/2, 2, 37/16, 1/8, -11/16 ], P![]> *];
models[[Integers()|1,2,17,34]] := [* <2, P![ -1/6, -8/9, -23/16, -13/72, 173/144, 13/36, -11/36 ], P![]> *];
models[[Integers()|1]] := [* <9, "CRV", [ Strings() | "11/16*s^4 - 1/8*s^3*z - 37/16*s^2*z^2 - 2*s*z^3 - 1/2*z^4 + y1^2", "11/36*s^6 - 13/36*s^5*z - 173/144*s^4*z^2 + 13/72*s^3*z^3 + 23/16*s^2*z^4 + 8/9*s*z^5 + 1/6*z^6 + y2^2", "-1/4*s^4 + 3/4*s^3*z - 1/16*s^2*z^2 - 7/8*s*z^3 + 3/16*z^4 + y3^2" ]> *];
models[[Integers()|1,2]] := [* <5, "CRV", [ Strings() | "11/16*s^6 - 57/16*s^5*z + 299/64*s^4*z^2 + 49/32*s^3*z^3 - 241/64*s^2*z^4 - s*z^5 + 3/8*z^6 + y1^2", "-s^4 + s^3*z + 11/4*s^2*z^2 + 1/2*s*z^3 - 1/4*z^4 + y2^2" ]> *];
models[[Integers()|1,2,5,10]] := [* <2, P![ -3/8, 1, 241/64, -49/32, -299/64, 57/16, -11/16 ], P![]> *];
models[[Integers()|1,5]] := [* <4, P![ -289/4096, 0, 957/8192, 0, -1311/32768, 0, -239/32768, 0, 103/1048576, 0, -5/2097152 ], P![]> *];
models[[Integers()|1,34]] := [* <3, P![ -11/256, -5/32, 25/128, 25/16, 647/256, 25/16, 25/128, -5/32, -11/256 ], P![]> *];
models[[Integers()|1,5,17,85]] := [* <2, P![ 1/8, 0, -123/64, -105/32, 25/64, 35/16, -11/16 ], P![]> *];
models[[Integers()|1,85]] := [* <4, "CRV", [ Strings() | "11/16*s^6 - 35/16*s^5*z - 25/64*s^4*z^2 + 105/32*s^3*z^3 + 123/64*s^2*z^4 - 1/8*z^6 + y1^2", "11/16*s^4 - 1/8*s^3*z - 37/16*s^2*z^2 - 2*s*z^3 - 1/2*z^4 + y2^2" ]> *];
models[[Integers()|1,170]] := [* <2, P![ 4/625, 24/625, 28/625, -72/625, -236/625, -48/125, -16/125 ], P![]> *];
models[[Integers()|1,17]] := [* <5, "CRV", [ Strings() | "11/16*s^6 - 35/16*s^5*z - 25/64*s^4*z^2 + 105/32*s^3*z^3 + 123/64*s^2*z^4 - 1/8*z^6 + y1^2", "-1/4*s^4 + 3/4*s^3*z - 1/16*s^2*z^2 - 7/8*s*z^3 + 3/16*z^4 + y2^2" ]> *];
models[[Integers()|1,2,85,170]] := [* <1, P![ 1/4, -1/2, -11/4, -1, 1 ], P![]> *];
models[[Integers()|1,10,17,170]] := [* <1, P![ -3/16, 7/8, 1/16, -3/4, 1/4 ], P![]> *];
