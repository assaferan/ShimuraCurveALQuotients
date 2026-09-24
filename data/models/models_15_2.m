// Subhyperelliptic cover models for X_0(15,2)*
//
// REGENERATED 2026-09-07 with DEFAULT FLAGS -- no flag is needed any more:
//     NORMALIZ_BIN=... magma -b D_s:=15 N_s:=2 OUTDIR:=... genmodels.m < /dev/null
// Populated covers went 12 -> 15.
//
// ⚠ WHAT ACTUALLY UNLOCKED THIS, because it is easy to misattribute: the COPRIME-TO-LEVEL CM
// FILTER becoming OFF BY DEFAULT (same day). NOT Y2TWIST. That flag was written for exactly these
// three bases, but a controlled run -- default vs the selector disabled, on the SAME code -- gives
// IDENTICAL output at all three, and the deferral path logs zero "unpinned y2-scale" messages.
// The selector never fires here now, so Y2TWIST stays off by default; see EquationsCovers.m.
// (An earlier evaluation compared against the COMMITTED files, which predate the coprime flip, and
// so credited the gains to the wrong flag. Compare against a current baseline, not an artifact.)
//
// ⚠ The p | gcd(d,N) local factor still has NO live implementation, so what makes this file
// trustworthy is the published Guo-Yang equation, not regeneration. See data/models/PROVENANCE.md.
//
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,10]] := [* <1, P![ -25/48, 0, -17/8, 0, -27/16 ], P![]> *];
models[[Integers()|1,15]] := [* <0, P![ -8/9, 0, 1/9 ], P![]>, <0, P![ 8, 0, 9 ], P![]>, <0, P![ -8, 0, 8 ], P![]> *];
models[[Integers()|1,3,10,30]] := [* <0, P![ 5/144, 7/9, -4/3 ], P![]> *];
models[[Integers()|1]] := [* <3, P![ -7/9, -13/9, -85/36, -103/36, -371/144, -19/12, -47/72, -1/6, -1/48 ], P![]> *];
models[[Integers()|1,3]] := [* <1, P![ -7/9, -38/9, -53/9, -10/3, -5/3 ], P![]>, <1, P![ 5/144, 0, -61/72, 0, -25/48 ], P![]> *];
models[[Integers()|1,2]] := [* <2, P![ -5/2, 0, -107/16, 0, 19/8, 0, -3/16 ], P![]> *];
models[[Integers()|1,2,5,10]] := [* <1, P![ -5/2, -107/2, 152, -96 ], P![]> *];
models[[Integers()|1,5]] := [* <2, P![ -5/2, 0, 127/2, 0, -47/2, 0, -75/2 ], P![]> *];
models[[Integers()|1,6,10,15]] := [* <0, P![ -8/9, 8/9 ], P![]> *];
models[[Integers()|1,2,15,30]] := [* <0, P![ 0, 8 ], P![]> *];
models[[Integers()|1,6]] := [* <2, P![ -25/6, 0, -347/16, 0, -261/8, 0, -243/16 ], P![]> *];
models[[Integers()|1,2,3,6]] := [* <1, P![ 0, -20, -428, 1216, -768 ], P![]> *];
models[[Integers()|1,30]] := [* <1, P![ 5/144, 0, 7/72, 0, -1/48 ], P![]> *];
models[[Integers()|1,3,5,15]] := [* <0, P![ 0, -64/9, 64/9 ], P![]> *];
models[[Integers()|1,5,6,30]] := [* <1, P![ 0, 5/18, 56/9, -32/3 ], P![]> *];
