// Subhyperelliptic cover models for X_0(22,5)*
//
// REGENERATED 2026-09-07 with DEFAULT FLAGS -- no flag is needed any more:
//     NORMALIZ_BIN=... magma -b D_s:=22 N_s:=5 OUTDIR:=... genmodels.m < /dev/null
// Populated covers went 3 -> 11.
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
// ✅ REGENERATED 2026-09-09: the four previously-EMPTY keys ({1}, {1,2}, {1,5}, {1,11}) are now
// filled. What unlocked them is EquationsByRebase (EquationsCovers.m), the last stage of
// AllEquationsAboveCovers: this base's genus-2 quotients need a degree-3 equation over a shared
// base and the degrees produced were 1,2,4,6,7,8, so the fibre-product construction had nothing to
// work with. Changing the Hauptmodul on the star base (t -> r + 1/u at a rational root r) fixes
// the degree profile. No flag needed; no previously-present entry changed.
//
// ⚠ The W={1} entry below is the r=0 model. It is ISOMORPHIC to Guo-Yang's published degree-12
// curve but not equal to it; the r=4/5 sweep reproduces their polynomial VERBATIM, which
// tests/_offline/FullCurve_22_5.m pins. Both are correct models of the same curve.
// ⚠ All the new entries are checked against Guo-Yang by tests/GuoYangQuotients_22_5.m.

models := AssociativeArray();
models[[Integers()|1,10]] := [* <3, P![ -1024/625, -4096/625, -6803/625, -6073/625, -3147/625, -951/625, -157/625, -11/625 ], P![]> *];
models[[Integers()|1,110]] := [* <0, P![ 1/4, 1/4 ], P![]>, <0, P![ 1/25, 0, -4/25 ], P![]>, <0, P![ 4, 0, 1 ], P![]> *];
models[[Integers()|1,55]] := [* <3, P![ -11/390625, 0, 6/78125, 0, 37/390625, 0, 56/390625, 0, 16/78125 ], P![]> *];
models[[Integers()|1,2,11,22]] := [* <1, P![ -4096/625, 20044/625, -36799/625, 6008/125, -368/25 ], P![]> *];
models[[Integers()|1]] := [* <5, P![ -43151/125, -40128/25, -2187824/625, -591552/125, -551976/125, -374424/125, -946538/625, -359424/625, -20352/125, -4192/125, -2984/625, -264/625, -11/625 ], P![]> *];
models[[Integers()|1,2]] := [* <2, P![ -11, 0, -56/25, 0, -304/625, 0, -256/3125 ], P![]> *];
models[[Integers()|1,2,5,10]] := [* <1, P![ -1024/625, 951/125, -8267/625, 6376/625, -368/125 ], P![]> *];
models[[Integers()|1,5]] := [* <2, P![ -11/625, 0, -73/5000, 0, -659/160000, 0, -1/2500 ], P![]> *];
models[[Integers()|1,5,22,110]] := [* <0, P![ 0, -16, 16 ], P![]> *];
models[[Integers()|1,5,11,55]] := [* <1, P![ 0, 1024/625, -3731/625, 4536/625, -368/125 ], P![]> *];
models[[Integers()|1,10,22,55]] := [* <2, P![ 0, 65536/50625, -77248/10125, 909488/50625, -1069424/50625, 13952/1125, -5888/2025 ], P![]> *];
models[[Integers()|1,2,55,110]] := [* <0, P![ 1, -9/4, 5/4 ], P![]> *];
models[[Integers()|1,10,11,110]] := [* <0, P![ 0, -4, 5 ], P![]> *];
models[[Integers()|1,11]] := [* <2, P![ 1/2000, 0, -17/5000, 0, 109/10000, 0, -16/625 ], P![]> *];
models[[Integers()|1,22]] := [* <3, P![ -4096/625, 0, -732/125, 0, -1243/625, 0, -38/125, 0, -11/625 ], P![]> *];
