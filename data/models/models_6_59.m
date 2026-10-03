// Subhyperelliptic cover models for X_0(6,59)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,354]] := [* <2, P![ -1/16, 0, 47/36, 0, -37/144, 0, 1/72 ], P![]> *];
models[[Integers()|1,3,118,354]] := [* <0, P![ 5, -4 ], P![]> *];
models[[Integers()|1,3]] := [* <6, P![ 845883/4096, 0, -551079/128, 0, 2757793/4096, 0, -147635/2048, 0, 67773/4096, 0, -2101/1024, 0, 391/4096, 0, -3/2048 ], P![]> *];
models[[Integers()|1,6,118,177]] := [* <1, P![ -512, 128, 20, -20, -3 ], P![]> *];
models[[Integers()|1,2,3,6]] := [* <3, P![ -40960, 22528, 32832, -4512, -19020, 5268, 1081, -604, -96 ], P![]> *];
models[[Integers()|1,2,177,354]] := [* <1, P![ 16/9, 8/9, -7/9, -8/9 ], P![]> *];
models[[Integers()|1,6,59,354]] := [* <1, P![ 80, -24, -67, -12, 32 ], P![]> *];
models[[Integers()|1,2,59,118]] := [* <2, P![ -2560, 2688, -412, -180, 65, 12 ], P![]> *];
models[[Integers()|1,118]] := [* <3, P![ -93987/256, 0, -973/64, 0, -665/128, 0, 35/64, 0, -3/256 ], P![]> *];
models[[Integers()|1,3,59,177]] := [* <3, P![ -8192, -2048, 4928, 3040, -1372, -44, 181, 24 ], P![]> *];

// Built 2026-10-03 as fibre products over the star line of committed double covers
// (FibreProductCovers.m); each verified against the trace formula at three primes, F_p and F_{p^2}.
// Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
models[[Integers()|1]] := [* <11, "CRV", [ Strings() | "3*s^4 + 20*s^3*z - 20*s^2*z^2 - 128*s*z^3 + 512*z^4 + y1^2", "8/9*s^3*z + 7/9*s^2*z^2 - 8/9*s*z^3 - 16/9*z^4 + y2^2", "4*s*z - 5*z^2 + y3^2" ]> *];
models[[Integers()|1,2]] := [* <6, "CRV", [ Strings() | "-12*s^5*z - 65*s^4*z^2 + 180*s^3*z^3 + 412*s^2*z^4 - 2688*s*z^5 + 2560*z^6 + y1^2", "8/9*s^3*z + 7/9*s^2*z^2 - 8/9*s*z^3 - 16/9*z^4 + y2^2" ]> *];
models[[Integers()|1,6]] := [* <5, "CRV", [ Strings() | "-32*s^4 + 12*s^3*z + 67*s^2*z^2 + 24*s*z^3 - 80*z^4 + y1^2", "3*s^4 + 20*s^3*z - 20*s^2*z^2 - 128*s*z^3 + 512*z^4 + y2^2" ]> *];
models[[Integers()|1,177]] := [* <5, "CRV", [ Strings() | "8/9*s^3*z + 7/9*s^2*z^2 - 8/9*s*z^3 - 16/9*z^4 + y1^2", "-24*s^7*z - 181*s^6*z^2 + 44*s^5*z^3 + 1372*s^4*z^4 - 3040*s^3*z^5 - 4928*s^2*z^6 + 2048*s*z^7 + 8192*z^8 + y2^2" ]> *];
models[[Integers()|1,59]] := [* <6, "CRV", [ Strings() | "-12*s^5*z - 65*s^4*z^2 + 180*s^3*z^3 + 412*s^2*z^4 - 2688*s*z^5 + 2560*z^6 + y1^2", "-32*s^4 + 12*s^3*z + 67*s^2*z^2 + 24*s*z^3 - 80*z^4 + y2^2" ]> *];
