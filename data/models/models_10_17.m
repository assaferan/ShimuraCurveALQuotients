// Subhyperelliptic cover models for X_0(10,17)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,10]] := [* <4, P![ -10625/1048576, 0, -3823/524288, 0, 15169/13107200, 0, -307/3276800, 0, 19/5242880, 0, -1/13107200 ], P![]> *];
models[[Integers()|1,5,34,170]] := [* <1, P![ 25/16, 75/16, -275/64, 125/16 ], P![]> *];
models[[Integers()|1,10,34,85]] := [* <2, P![ -7/256, -3/128, 41/512, -165/512, 1125/4096, -125/512 ], P![]> *];
models[[Integers()|1,2,17,34]] := [* <1, P![ -175/64, -1025/64, -3575/256, 375/64, -625/8 ], P![]> *];
models[[Integers()|1,2]] := [* <3, P![ -26875/256, -55625/256, -94375/1024, 71875/256, 161875/512, -71875/256, -94375/1024, 55625/256, -26875/256 ], P![]> *];
models[[Integers()|1,2,5,10]] := [* <2, P![ -35/256, -85/128, -35/512, -5/512, -20775/4096, 4375/1024, -625/128 ], P![]> *];
models[[Integers()|1,5,17,85]] := [* <1, P![ -35/64, -65/64, 325/256, -125/32 ], P![]> *];
models[[Integers()|1,170]] := [* <1, P![ 17/1024, 0, -13/12800, 0, 1/25600 ], P![]> *];
models[[Integers()|1,17]] := [* <2, P![ -625/4096, 0, -61/512, 0, 43/4096, 0, -1/2048 ], P![]> *];
models[[Integers()|1,2,85,170]] := [* <0, P![ 1/80, -1/80, 1/64 ], P![]> *];
models[[Integers()|1,10,17,170]] := [* <0, P![ 5, 20 ], P![]> *];

// Built 2026-10-03 as fibre products over the star line of committed double covers
// (FibreProductCovers.m); each verified against the trace formula at three primes, F_p and F_{p^2}.
// Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
models[[Integers()|1]] := [* <7, "CRV", [ Strings() | "125/512*s^5*z - 1125/4096*s^4*z^2 + 165/512*s^3*z^3 - 41/512*s^2*z^4 + 3/128*s*z^5 + 7/256*z^6 + y1^2", "625/8*s^4 - 375/64*s^3*z + 3575/256*s^2*z^2 + 1025/64*s*z^3 + 175/64*z^4 + y2^2", "-20*s*z - 5*z^2 + y3^2" ]> *];
models[[Integers()|1,34]] := [* <4, "CRV", [ Strings() | "-125/16*s^3*z + 275/64*s^2*z^2 - 75/16*s*z^3 - 25/16*z^4 + y1^2", "625/8*s^4 - 375/64*s^3*z + 3575/256*s^2*z^2 + 1025/64*s*z^3 + 175/64*z^4 + y2^2" ]> *];
models[[Integers()|1,5]] := [* <4, "CRV", [ Strings() | "125/32*s^3*z - 325/256*s^2*z^2 + 65/64*s*z^3 + 35/64*z^4 + y1^2", "625/128*s^6 - 4375/1024*s^5*z + 20775/4096*s^4*z^2 + 5/512*s^3*z^3 + 35/512*s^2*z^4 + 85/128*s*z^5 + 35/256*z^6 + y2^2" ]> *];
models[[Integers()|1,85]] := [* <3, "CRV", [ Strings() | "125/32*s^3*z - 325/256*s^2*z^2 + 65/64*s*z^3 + 35/64*z^4 + y1^2", "125/512*s^5*z - 1125/4096*s^4*z^2 + 165/512*s^3*z^3 - 41/512*s^2*z^4 + 3/128*s*z^5 + 7/256*z^6 + y2^2" ]> *];
