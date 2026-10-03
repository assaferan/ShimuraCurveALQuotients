// Subhyperelliptic cover models for X_0(62,3)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,62]] := [* <2, P![ -9/4096, 0, -1/512, 0, -29/4096, 0, -1/2048 ], P![]> *];
models[[Integers()|1,6,31,186]] := [* <2, P![ 0, -155/65536, 29175/16384, -17625/32768, -625/16384, 3125/65536 ], P![]> *];
models[[Integers()|1,186]] := [* <3, P![ -31/65536, 0, 1167/16384, 0, -141/32768, 0, -1/16384, 0, 1/65536 ], P![]> *];
models[[Integers()|1,2,93,186]] := [* <1, P![ -31/65536, 5835/16384, -3525/32768, -125/16384, 625/65536 ], P![]> *];
models[[Integers()|1,6,62,93]] := [* <1, P![ -9/4096, -5/512, -725/4096, -125/2048 ], P![]> *];
models[[Integers()|1,3]] := [*  *];
models[[Integers()|1,3,62,186]] := [* <0, P![ 0, 5 ], P![]> *];
models[[Integers()|1,2,31,62]] := [* <1, P![ 0, -45/4096, -25/512, -3625/4096, -625/2048 ], P![]> *];

// Built 2026-10-03 as fibre products over the star line of committed double covers
// (FibreProductCovers.m); each verified against the trace formula at three primes, F_p and F_{p^2}.
// Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
models[[Integers()|1]] := [* <11, "CRV", [ Strings() | "-625/65536*s^4 + 125/16384*s^3*z + 3525/32768*s^2*z^2 - 5835/16384*s*z^3 + 31/65536*z^4 + y1^2", "125/2048*s^3*z + 725/4096*s^2*z^2 + 5/512*s*z^3 + 9/4096*z^4 + y2^2", "-5*s*z + y3^2" ]> *];
models[[Integers()|1,2]] := [* <5, "CRV", [ Strings() | "-625/65536*s^4 + 125/16384*s^3*z + 3525/32768*s^2*z^2 - 5835/16384*s*z^3 + 31/65536*z^4 + y1^2", "625/2048*s^4 + 3625/4096*s^3*z + 25/512*s^2*z^2 + 45/4096*s*z^3 + y2^2" ]> *];
models[[Integers()|1,6]] := [* <6, "CRV", [ Strings() | "125/2048*s^3*z + 725/4096*s^2*z^2 + 5/512*s*z^3 + 9/4096*z^4 + y1^2", "-3125/65536*s^5*z + 625/16384*s^4*z^2 + 17625/32768*s^3*z^3 - 29175/16384*s^2*z^4 + 155/65536*s*z^5 + y2^2" ]> *];
models[[Integers()|1,93]] := [* <5, "CRV", [ Strings() | "-625/65536*s^4 + 125/16384*s^3*z + 3525/32768*s^2*z^2 - 5835/16384*s*z^3 + 31/65536*z^4 + y1^2", "125/2048*s^3*z + 725/4096*s^2*z^2 + 5/512*s*z^3 + 9/4096*z^4 + y2^2" ]> *];
models[[Integers()|1,31]] := [* <6, "CRV", [ Strings() | "625/2048*s^4 + 3625/4096*s^3*z + 25/512*s^2*z^2 + 45/4096*s*z^3 + y1^2", "-3125/65536*s^5*z + 625/16384*s^4*z^2 + 17625/32768*s^3*z^3 - 29175/16384*s^2*z^4 + 155/65536*s*z^5 + y2^2" ]> *];
