// Subhyperelliptic cover models for X_0(6,41)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,3,6]] := [* <2, P![ -405/2, -1539/2, -74439/64, -111213/64, -735399/256, -156411/64, -729 ], P![]> *];
models[[Integers()|1,6,82,123]] := [* <1, P![ -8, -32, -163/4, -16 ], P![]> *];
models[[Integers()|1,2,123,246]] := [* <0, P![ 1/16, -1/16, 9/64 ], P![]> *];
models[[Integers()|1,246]] := [* <1, P![ 369/1024, 0, -53/512, 0, 9/1024 ], P![]> *];
models[[Integers()|1,2,41,82]] := [* <1, P![ -3240, -15552, -107487/4, -19683, -5184 ], P![]> *];
models[[Integers()|1,3]] := [* <4, P![ -807003/65536, 0, -554769/16384, 0, 1492911/32768, 0, -299781/16384, 0, 193509/65536, 0, -729/4096 ], P![]> *];
models[[Integers()|1,2]] := [* <3, P![ -1128492, 1320948, 525609, -872856, -106758, 207684, 23409, -18144, -3240 ], P![]> *];
models[[Integers()|1,3,82,246]] := [* <0, P![ 5, 4 ], P![]> *];
models[[Integers()|1,82]] := [* <2, P![ -27/64, 0, -41/32, 0, 77/64, 0, -1/4 ], P![]> *];
models[[Integers()|1,6,41,246]] := [* <1, P![ 45/16, -9/16, 261/64, 81/16 ], P![]> *];
models[[Integers()|1,3,41,123]] := [* <2, P![ -81/2, -243/2, -8667/64, -15309/64, -98091/256, -729/4 ], P![]> *];

// Built 2026-10-03 as fibre products over the star line of committed double covers
// (FibreProductCovers.m); each verified against the trace formula at three primes, F_p and F_{p^2}.
// Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
models[[Integers()|1]] := [* <7, "CRV", [ Strings() | "729*s^6 + 156411/64*s^5*z + 735399/256*s^4*z^2 + 111213/64*s^3*z^3 + 74439/64*s^2*z^4 + 1539/2*s*z^5 + 405/2*z^6 + y1^2", "16*s^3*z + 163/4*s^2*z^2 + 32*s*z^3 + 8*z^4 + y2^2", "5184*s^4 + 19683*s^3*z + 107487/4*s^2*z^2 + 15552*s*z^3 + 3240*z^4 + y3^2" ]> *];
models[[Integers()|1,6]] := [* <4, "CRV", [ Strings() | "16*s^3*z + 163/4*s^2*z^2 + 32*s*z^3 + 8*z^4 + y1^2", "-81/16*s^3*z - 261/64*s^2*z^2 + 9/16*s*z^3 - 45/16*z^4 + y2^2" ]> *];
models[[Integers()|1,123]] := [* <3, "CRV", [ Strings() | "16*s^3*z + 163/4*s^2*z^2 + 32*s*z^3 + 8*z^4 + y1^2", "729/4*s^5*z + 98091/256*s^4*z^2 + 15309/64*s^3*z^3 + 8667/64*s^2*z^4 + 243/2*s*z^5 + 81/2*z^6 + y2^2" ]> *];
models[[Integers()|1,41]] := [* <4, "CRV", [ Strings() | "-81/16*s^3*z - 261/64*s^2*z^2 + 9/16*s*z^3 - 45/16*z^4 + y1^2", "5184*s^4 + 19683*s^3*z + 107487/4*s^2*z^2 + 15552*s*z^3 + 3240*z^4 + y2^2" ]> *];
