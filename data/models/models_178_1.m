// Subhyperelliptic cover models for X_0(178,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2]] := [* <4, P![ -110592, 2129920/3, -17779712/9, 29249536/9, -294268928/81, 709578752/243, -139225088/81, 1604882432/2187, -478117888/2187, 810070016/19683, -228502528/59049 ], P![]> *];
models[[Integers()|1,89]] := [* <1, P![ -432, 1408/3, -3616/9, 1408/9, -2608/81 ], P![]> *];
models[[Integers()|1,178]] := [* <2, P![ 256, -4096/3, 25664/9, -82688/27, 155264/81, -168704/243, 87616/729 ], P![]> *];

// Built 2026-10-03 as fibre products over the star line of committed double covers
// (FibreProductCovers.m); each verified against the trace formula at three primes, F_p and F_{p^2}.
// Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
models[[Integers()|1]] := [* <7, "CRV", [ Strings() | "-87616/729*s^6 + 168704/243*s^5*z - 155264/81*s^4*z^2 + 82688/27*s^3*z^3 - 25664/9*s^2*z^4 + 4096/3*s*z^5 - 256*z^6 + y1^2", "2608/81*s^4 - 1408/9*s^3*z + 3616/9*s^2*z^2 - 1408/3*s*z^3 + 432*z^4 + y2^2" ]> *];
