// Subhyperelliptic cover models for X_0(6,47)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,6,47,282]] := [* <2, P![ -9/8, 3/8, 59/16, -5/2, -19/16, 3/4, 1/4 ], P![]> *];
models[[Integers()|1,3,94,282]] := [* <0, P![ -8, 8, 4 ], P![]> *];
models[[Integers()|1,3]] := [* <5, P![ -43/1296, -121/648, -605/1296, -3497/648, -57329/1296, -56441/324, -264863/648, -206717/324, -885401/1296, -106855/216, -32981/144, -429/8, -67/16 ], P![]> *];
models[[Integers()|1,2,47,94]] := [* <1, P![ -10, 38/3, 173/9, -164/9, -76/9 ], P![]> *];
models[[Integers()|1,282]] := [* <3, P![ 1/4, -1/2, 9/4, 57/2, 77, 213/2, 329/4, 51/2, 9/4 ], P![]> *];
models[[Integers()|1,2,141,282]] := [* <1, P![ 9/16, 3/8, -19/16, 1/4, 1/4 ], P![]> *];
models[[Integers()|1,6]] := [* <5, "CRV", [ Strings() | "y^2 - 1/4*s^6 - 3/4*s^5*z + 19/16*s^4*z^2 + 5/2*s^3*z^3 - 59/16*s^2*z^4 - 3/8*s*z^5 + 9/8*z^6", "x^2 + 19*s^2 + 3*s*z - 45/4*z^2" ]> *];
models[[Integers()|1,6,94,141]] := [* <0, P![ 45/4, -3, -19 ], P![]> *];
models[[Integers()|1,94]] := [* <1, P![ -43/4, -82, -437/2, -222, -603/4 ], P![]> *];
models[[Integers()|1,2,3,6]] := [* <3, P![ -5/72, 1/24, 439/1296, -247/972, -4859/11664, 238/729, 185/1458, -20/243, -19/729 ], P![]> *];
models[[Integers()|1,3,47,141]] := [* <2, P![ 5/64, 1/32, -179/576, -1/108, 197/648, -11/162, -19/324 ], P![]> *];
models[[Integers()|1,141]] := [* <3, "CRV", [ Strings() | "y^2 - 1/4*s^4 - 1/4*s^3*z + 19/16*s^2*z^2 - 3/8*s*z^3 - 9/16*z^4", "x^2 + 19*s^2 + 3*s*z - 45/4*z^2" ]> *];

// Built 2026-10-03 as fibre products over the star line of committed double covers
// (FibreProductCovers.m); each verified against the trace formula at three primes, F_p and F_{p^2}.
// Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
models[[Integers()|1]] := [* <9, "CRV", [ Strings() | "76/9*s^4 + 164/9*s^3*z - 173/9*s^2*z^2 - 38/3*s*z^3 + 10*z^4 + y1^2", "19*s^2 + 3*s*z - 45/4*z^2 + y2^2", "-1/4*s^4 - 1/4*s^3*z + 19/16*s^2*z^2 - 3/8*s*z^3 - 9/16*z^4 + y3^2" ]> *];
models[[Integers()|1,2]] := [* <5, "CRV", [ Strings() | "19/729*s^8 + 20/243*s^7*z - 185/1458*s^6*z^2 - 238/729*s^5*z^3 + 4859/11664*s^4*z^4 + 247/972*s^3*z^5 - 439/1296*s^2*z^6 - 1/24*s*z^7 + 5/72*z^8 + y1^2", "-1/4*s^4 - 1/4*s^3*z + 19/16*s^2*z^2 - 3/8*s*z^3 - 9/16*z^4 + y2^2" ]> *];
models[[Integers()|1,47]] := [* <5, "CRV", [ Strings() | "76/9*s^4 + 164/9*s^3*z - 173/9*s^2*z^2 - 38/3*s*z^3 + 10*z^4 + y1^2", "19/324*s^6 + 11/162*s^5*z - 197/648*s^4*z^2 + 1/108*s^3*z^3 + 179/576*s^2*z^4 - 1/32*s*z^5 - 5/64*z^6 + y2^2" ]> *];
