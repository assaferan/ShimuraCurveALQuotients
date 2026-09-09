// Subhyperelliptic cover models for X_0(14,5)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
// ✅ REGENERATED 2026-09-09: previously EMPTY cover keys are now filled, unlocked by
// EquationsByRebase (EquationsCovers.m). No flag needed. Verified against the Guo-Yang quotient
// oracle before installing; see data/models/PROVENANCE.md.

models := AssociativeArray();
models[[Integers()|1,10]] := [* <2, P![ -35/1048576, 0, -1803/1048576, 0, 367/1048576, 0, -1/1048576 ], P![]> *];
models[[Integers()|1,2,7,14]] := [* <0, P![ -11, 216, -1024 ], P![]> *];
models[[Integers()|1,14]] := [* <0, P![ 5/17, 150/289, 5/17 ], P![]>, <0, P![ 5, 0, -16 ], P![]>, <0, P![ 5/16, 0, -1/16 ], P![]> *];
models[[Integers()|1,70]] := [* <1, P![ -7/256, 0, 11/128, 0, 1/256 ], P![]> *];
models[[Integers()|1]] := [* <3, P![ -14375/65536, -22875/16384, -124175/32768, -23565/4096, -352037/65536, -3249/1024, -38075/32768, -3963/16384, -1439/65536 ], P![]> *];
models[[Integers()|1,2]] := [* <1, P![ -1/578, 0, 9/4624, 0, -1/578 ], P![]>, <1, P![ -4, 16, -15, -2, -11 ], P![]> *];
models[[Integers()|1,2,5,10]] := [* <1, P![ -11/16, 87/4, -915/4, 822, -256 ], P![]> *];
models[[Integers()|1,5]] := [* <2, P![ -35/256, 0, 111/128, 0, -347/256, 0, -1/16 ], P![]> *];
models[[Integers()|1,35]] := [* <1, P![ -7/65536, 0, -181/32768, 0, 1/65536 ], P![]> *];
models[[Integers()|1,7]] := [* <2, P![ -5/9826, -75/83521, 5/78608, 675/668168, 5/78608, -75/83521, -5/9826 ], P![]> *];
models[[Integers()|1,5,7,35]] := [* <1, P![ -11/16, 65/4, -395/4, 32 ], P![]> *];
models[[Integers()|1,2,35,70]] := [* <0, P![ 1/16, -3/4, 1/4 ], P![]> *];
models[[Integers()|1,5,14,70]] := [* <0, P![ 1, -8 ], P![]> *];
models[[Integers()|1,10,14,35]] := [* <0, P![ -11, 128 ], P![]> *];
models[[Integers()|1,7,10,70]] := [* <1, P![ 1/16, -5/4, 25/4, -2 ], P![]> *];
