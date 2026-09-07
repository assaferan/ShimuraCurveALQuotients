// Subhyperelliptic cover models for X_0(22,3)*
//
// REGENERATED 2026-09-07 with DEFAULT FLAGS -- no flag is needed any more:
//     NORMALIZ_BIN=... magma -b D_s:=22 N_s:=3 OUTDIR:=... genmodels.m < /dev/null
// Populated covers went 13 -> 15.
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
models[[Integers()|1,3,22,66]] := [* <0, P![ 0, -729/4, 729/4 ], P![]> *];
models[[Integers()|1,66]] := [* <0, P![ 9, 0, 4 ], P![]>, <0, P![ 9/4, 0, -9/4 ], P![]>, <0, P![ -9/4, 0, 1/4 ], P![]> *];
models[[Integers()|1,2,11,22]] := [* <1, P![ -11/186624, 35/373248, -131/2985984, 1/110592 ], P![]> *];
models[[Integers()|1]] := [* <3, P![ -3/4096, 0, -77/9216, 0, -1073/18432, 0, -77/9216, 0, -3/4096 ], P![]> *];
models[[Integers()|1,6,11,66]] := [* <0, P![ 9, -9 ], P![]> *];
models[[Integers()|1,3]] := [* <1, P![ -11/2304, 0, 31/4608, 0, -11/4096 ], P![]> *];
models[[Integers()|1,33]] := [* <1, P![ -11/2304, 0, -13/10368, 0, -1/6912 ], P![]> *];
models[[Integers()|1,3,11,33]] := [* <0, P![ -11/2304, 13/4608, -3/4096 ], P![]> *];
models[[Integers()|1,2]] := [* <2, P![ -11/186624, 0, -35/839808, 0, -131/15116544, 0, -1/1259712 ], P![]> *];
models[[Integers()|1,2,33,66]] := [* <0, P![ 0, -9/4 ], P![]> *];
models[[Integers()|1,6]] := [* <2, P![ 11/16384, 0, -49/1327104, 0, -23/11943936, 0, -1/3981312 ], P![]> *];
models[[Integers()|1,6,22,33]] := [* <1, P![ 0, 11/9216, -13/18432, 3/16384 ], P![]> *];
models[[Integers()|1,2,3,6]] := [* <1, P![ 0, 11/746496, -35/1492992, 131/11943936, -1/442368 ], P![]> *];
models[[Integers()|1,22]] := [* <2, P![ -11/9216, 0, 53/18432, 0, -347/147456, 0, 11/16384 ], P![]> *];
models[[Integers()|1,11]] := [* <1, P![ -11/4096, 0, -25/165888, 0, -1/110592 ], P![]> *];
