// Subhyperelliptic cover models for X_0(57,2)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,57,114]] := [* <1, P![ 0, 27, -81/4, -243/8, 6561/256 ], P![]> *];
models[[Integers()|1,3,19,57]] := [* <1, P![ 0, 864, 81, -972 ], P![]> *];
models[[Integers()|1,3,38,114]] := [* <2, P![ 0, -1, 7/4, 3/8, -531/256, 243/256 ], P![]> *];
models[[Integers()|1,6]] := [* <4, P![ -1539/256, 0, -4293/32, 0, 7047/128, 0, -21951/64, 0, -32643/256, 0, -729/64 ], P![]> *];
models[[Integers()|1,2,3,6]] := [* <2, P![ -7776, 12879, 51759/4, -124659/4, 662661/256, 4822335/256, -531441/64 ], P![]> *];
models[[Integers()|1,6,19,114]] := [* <0, P![ -3, 3 ], P![]> *];
models[[Integers()|1,6,38,57]] := [* <2, P![ 2592, -1701, -24057/4, 4374, 898857/256, -177147/64 ], P![]> *];
models[[Integers()|1,114]] := [* <3, P![ 513/256, 0, -45/64, 0, 603/128, 0, 171/64, 0, 81/256 ], P![]> *];
models[[Integers()|1,19]] := [* <2, P![ -27, 0, -630, 0, -315, 0, -36 ], P![]> *];
models[[Integers()|1,2,19,38]] := [* <1, P![ 0, -5832, 21141/4, 28431/4, -6561 ], P![]> *];

// Built 2026-10-03 as fibre products over the star line of committed double covers
// (FibreProductCovers.m); each verified against the trace formula at three primes, F_p and F_{p^2}.
// Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
models[[Integers()|1]] := [* <9, "CRV", [ Strings() | "972*s^3*z - 81*s^2*z^2 - 864*s*z^3 + y1^2", "6561*s^4 - 28431/4*s^3*z - 21141/4*s^2*z^2 + 5832*s*z^3 + y2^2", "-6561/256*s^4 + 243/8*s^3*z + 81/4*s^2*z^2 - 27*s*z^3 + y3^2" ]> *];
models[[Integers()|1,2]] := [* <4, "CRV", [ Strings() | "6561*s^4 - 28431/4*s^3*z - 21141/4*s^2*z^2 + 5832*s*z^3 + y1^2", "531441/64*s^6 - 4822335/256*s^5*z - 662661/256*s^4*z^2 + 124659/4*s^3*z^3 - 51759/4*s^2*z^4 - 12879*s*z^5 + 7776*z^6 + y2^2" ]> *];
models[[Integers()|1,38]] := [* <5, "CRV", [ Strings() | "6561*s^4 - 28431/4*s^3*z - 21141/4*s^2*z^2 + 5832*s*z^3 + y1^2", "-243/256*s^5*z + 531/256*s^4*z^2 - 3/8*s^3*z^3 - 7/4*s^2*z^4 + s*z^5 + y2^2" ]> *];
models[[Integers()|1,3]] := [* <5, "CRV", [ Strings() | "972*s^3*z - 81*s^2*z^2 - 864*s*z^3 + y1^2", "-243/256*s^5*z + 531/256*s^4*z^2 - 3/8*s^3*z^3 - 7/4*s^2*z^4 + s*z^5 + y2^2" ]> *];
models[[Integers()|1,57]] := [* <4, "CRV", [ Strings() | "972*s^3*z - 81*s^2*z^2 - 864*s*z^3 + y1^2", "-6561/256*s^4 + 243/8*s^3*z + 81/4*s^2*z^2 - 27*s*z^3 + y2^2" ]> *];
