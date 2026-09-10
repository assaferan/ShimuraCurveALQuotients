// Subhyperelliptic cover models for X_0(10,13)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
// ✅ REGENERATED 2026-09-09: previously EMPTY cover keys are now filled, unlocked by
// EquationsByRebase (EquationsCovers.m). No flag needed. Verified against the Guo-Yang quotient
// oracle before installing; see data/models/PROVENANCE.md.

models := AssociativeArray();
models[[Integers()|1,10]] := [* <1, P![ 325/3136, 325/784, 1275/784, 475/196, 725/196 ], P![]>, <1, P![ -387/1024, -631/512, -717/1024, 259/512, -27/1024 ], P![]> *];
models[[Integers()|1,13]] := [* <1, P![ -43/256, -51/128, 11/256, 43/128, -27/256 ], P![]>, <1, P![ 5/2916, 0, 199/1458, 0, 8125/2916 ], P![]> *];
models[[Integers()|1,5,13,65]] := [* <0, P![ 0, -1/2, -27/16 ], P![]> *];
models[[Integers()|1,2,65,130]] := [* <0, P![ 0, 1/16, 1/64 ], P![]> *];
models[[Integers()|1]] := [* <3, "CRV", [ Strings() | "y^2 - 725/196*s^4 - 475/196*s^3*z - 1275/784*s^2*z^2 - 325/784*s*z^3 - 325/3136*z^4", "x^2 + 75/784*s^2 + 25/392*s*z + 25/784*z^2" ]> *];
models[[Integers()|1,10,13,130]] := [* <0, P![ 5/4, -4, 4 ], P![]> *];
models[[Integers()|1,2]] := [* <2, P![ -5/512, 0, 23/16, 0, 4800, 0, -2080000 ], P![]> *];
models[[Integers()|1,2,5,10]] := [* <1, P![ -5/2, -17/16, 1209/64, -89/4, -27/4 ], P![]> *];
models[[Integers()|1,5]] := [* <2, P![ -1/40000000, 0, -523/200000000, 0, -723/8000000, 0, -13/12800 ], P![]> *];
models[[Integers()|1,2,13,26]] := [* <1, P![ 0, -5/2, -7/16, 19, -27 ], P![]> *];
models[[Integers()|1,10,26,65]] := [* <0, P![ -1/2, -29/16, -27/64 ], P![]> *];
models[[Integers()|1,5,26,130]] := [* <1, P![ 0, 5/16, -59/64, 3/4, 1/4 ], P![]> *];
models[[Integers()|1,26]] := [* <2, P![ -13/12800, 0, -577/10000, 0, -608/625, 0, -3328/625 ], P![]> *];
models[[Integers()|1,130]] := [* <1, P![ 9/1024, 13/512, 7/1024, -9/512, 1/1024 ], P![]>, <1, P![ 5/64, 0, -37/8, 0, 325/4 ], P![]> *];
models[[Integers()|1,65]] := [* <0, P![ -25/784, -25/392, -75/784 ], P![]>, <0, P![ -1/1458, 0, -25/1458 ], P![]>, <0, P![ -1/32, 0, -25/16 ], P![]> *];
