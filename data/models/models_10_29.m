// Subhyperelliptic cover models for X_0(10,29)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,10,29,290]] := [* <1, P![ -700, 100, -75, 1000 ], P![]> *];
models[[Integers()|1,2,145,290]] := [* <0, P![ -7, 8 ], P![]> *];
models[[Integers()|1,10,58,145]] := [* <2, P![ 5, -5, -7/4, 3/2, 5/4, -2 ], P![]> *];
models[[Integers()|1,290]] := [* <1, P![ 725/64, 0, 51/32, 0, 5/64 ], P![]> *];
models[[Integers()|1,2]] := [* <6, P![ -18125/524288, 0, -10590825/1048576, 0, -313069/131072, 0, -434755/1048576, 0, -29777/524288, 0, -5295/1048576, 0, -63/262144, 0, -5/1048576 ], P![]> *];
models[[Integers()|1,2,5,10]] := [* <3, P![ -140, 160, 14, 166, -895/4, -27/2, 193/4, 56, -80 ], P![]> *];
models[[Integers()|1,145]] := [* <4, P![ -25/8192, 0, -14601/16384, 0, -175/2048, 0, -151/8192, 0, -15/8192, 0, -1/16384 ], P![]> *];
models[[Integers()|1,5]] := [* <5, P![ -35, -360, -1176, -648, 2968, 2264, -5354, -2264, 2968, 648, -1176, 360, -35 ], P![]> *];
models[[Integers()|1,2,29,58]] := [* <3, P![ 20, 0, -2, -26, 9/4, 9/2, -7/4, -10 ], P![]> *];
models[[Integers()|1,5,58,290]] := [* <0, P![ 4, 4, 5 ], P![]> *];
models[[Integers()|1,5,29,145]] := [* <2, P![ -35, 75, -111/4, -49/2, 13/4, 24, -16 ], P![]> *];

// Built 2026-10-03 as fibre products over the star line of committed double covers
// (FibreProductCovers.m); each verified against the trace formula at three primes, F_p and F_{p^2}.
// Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
models[[Integers()|1]] := [* <11, "CRV", [ Strings() | "80*s^8 - 56*s^7*z - 193/4*s^6*z^2 + 27/2*s^5*z^3 + 895/4*s^4*z^4 - 166*s^3*z^5 - 14*s^2*z^6 - 160*s*z^7 + 140*z^8 + y1^2", "10*s^7*z + 7/4*s^6*z^2 - 9/2*s^5*z^3 - 9/4*s^4*z^4 + 26*s^3*z^5 + 2*s^2*z^6 - 20*z^8 + y2^2", "-1000*s^3*z + 75*s^2*z^2 - 100*s*z^3 + 700*z^4 + y3^2" ]> *];
models[[Integers()|1,58]] := [* <5, "CRV", [ Strings() | "2*s^5*z - 5/4*s^4*z^2 - 3/2*s^3*z^3 + 7/4*s^2*z^4 + 5*s*z^5 - 5*z^6 + y1^2", "10*s^7*z + 7/4*s^6*z^2 - 9/2*s^5*z^3 - 9/4*s^4*z^4 + 26*s^3*z^5 + 2*s^2*z^6 - 20*z^8 + y2^2" ]> *];
models[[Integers()|1,10]] := [* <6, "CRV", [ Strings() | "2*s^5*z - 5/4*s^4*z^2 - 3/2*s^3*z^3 + 7/4*s^2*z^4 + 5*s*z^5 - 5*z^6 + y1^2", "-1000*s^3*z + 75*s^2*z^2 - 100*s*z^3 + 700*z^4 + y2^2" ]> *];
models[[Integers()|1,29]] := [* <6, "CRV", [ Strings() | "10*s^7*z + 7/4*s^6*z^2 - 9/2*s^5*z^3 - 9/4*s^4*z^4 + 26*s^3*z^5 + 2*s^2*z^6 - 20*z^8 + y1^2", "16*s^6 - 24*s^5*z - 13/4*s^4*z^2 + 49/2*s^3*z^3 + 111/4*s^2*z^4 - 75*s*z^5 + 35*z^6 + y2^2" ]> *];
