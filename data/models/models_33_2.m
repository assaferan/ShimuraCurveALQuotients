// Subhyperelliptic cover models for X_0(33,2)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
// Fibre-product entries (FibreProductCovers.m, 2026-10-14): built over the star line from the double covers of
// this file by tools/writefp.m with the coordinate check of the three covers above each pair of factors;
// each checked against the trace formula at the first three good primes, F_p and F_{p^2}, and again at the
// 4th and 5th.  Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,3,22,66]] := [* <0, P![ -1/3, 16/81 ], P![]> *];
models[[Integers()|1,66]] := [* <1, P![ 11/3, 0, 38, 0, 81 ], P![]> *];
models[[Integers()|1,2,11,22]] := [* <1, P![ 0, 64/27, -4288/729, 29696/6561, -65536/59049 ], P![]> *];
models[[Integers()|1]] := [* <5, "CRV", [ Strings() | "65536/6561*s^4 - 332800/6561*s^3*z + 68288/729*s^2*z^2 - 6016/81*s*z^3 + 64/3*z^4 + y1^2", "-16/81*s*z + 1/3*z^2 + y2^2", "-256/81*s^2 + 256/81*s*z + y3^2" ]> *];
models[[Integers()|1,6,11,66]] := [* <1, P![ 0, 64/3, -2752/81, 1024/81 ], P![]> *];
models[[Integers()|1,3]] := [* <2, P![ -11/9, 0, -86, 0, -963, 0, -2916 ], P![]> *];
models[[Integers()|1,33]] := [* <2, P![ 256/9, 0, -2560/81, 0, 256/729, 0, 2048/729 ], P![]> *];
models[[Integers()|1,3,11,33]] := [* <1, P![ 256/9, -6656/81, 55552/729, -16384/729 ], P![]> *];
models[[Integers()|1,2]] := [* <2, P![ -5632/59049, 0, 13/324, 0, 63/512, 0, -19683/262144 ], P![]> *];
models[[Integers()|1,2,33,66]] := [* <0, P![ 0, -256/81, 256/81 ], P![]> *];
models[[Integers()|1,6]] := [* <3, "CRV", [ Strings() | "16384/6561*s^3*z - 4352/729*s^2*z^2 + 256/81*s*z^3 + y1^2", "65536/6561*s^4 - 332800/6561*s^3*z + 68288/729*s^2*z^2 - 6016/81*s*z^3 + 64/3*z^4 + y2^2" ]> *];
models[[Integers()|1,6,22,33]] := [* <1, P![ 0, -256/81, 4352/729, -16384/6561 ], P![]> *];
models[[Integers()|1,2,3,6]] := [* <1, P![ -64/3, 6016/81, -68288/729, 332800/6561, -65536/6561 ], P![]> *];
models[[Integers()|1,22]] := [* <2, P![ -1/3, 0, -22, 0, -171, 0, -324 ], P![]> *];
models[[Integers()|1,11]] := [* <3, "CRV", [ Strings() | "16384/729*s^3*z - 55552/729*s^2*z^2 + 6656/81*s*z^3 - 256/9*z^4 + y1^2", "-1024/81*s^3*z + 2752/81*s^2*z^2 - 64/3*s*z^3 + y2^2" ]> *];
