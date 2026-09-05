// Subhyperelliptic cover models for X_0(39,2)*
//
// ⚠ PRODUCED WITH CMNONCOPRIME=1 -- this file does NOT regenerate under default settings.
// With the coprime-to-level CM filter on, 39_2 sees 3 CM points against demand 19 and dies with
// "Could not find enough points"; with the filter off it sees 24 and builds cleanly (15 keys, 0
// empty). Regenerate with:
//     CMNONCOPRIME=1 NORMALIZ_BIN=... magma -b D_s:=39 N_s:=2 OUTDIR:=... genmodels.m
// tests/_offline/ModelRegen.m therefore lists 39_2 in MR_KNOWN_DRIFT; that is expected, not rot.
//
// ⚠ The flag is NOT known to be safe in general -- the p | gcd(d,N) local factor has no live
// implementation (kappaminuszero is dead code), and at 26_3 two non-coprime discriminants give
// provably wrong values. What makes THIS file trustworthy is not the flag but an INDEPENDENT
// oracle: its W={1} genus-7 curve is IsIsomorphic to Guo-Yang's published equation
//   y^2 = -(x^8+11x^7+52x^6+140x^5+243x^4+280x^3+208x^2+88x+16)(7x^4+24x^3+32x^2+24x+16)
//          (x^4+3x^3+8x^2+12x+7)
// (verified 2026-09-05, tests/GuoYangEquations.m), and ModelChecks validates it structurally.
// Any further base produced this way must clear the same bar before being committed.
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,13]] := [* <3, P![ -27/4096, 0, -1585/173056, 0, -641/346112, 0, -25/173056, 0, -3/692224 ], P![]> *];
models[[Integers()|1]] := [* <7, P![ -63/169, 90/169, -261/169, 900/169, -1161/169, 360/169, -4239/169, -450/169, -360/13, 450/169, -4239/169, -360/169, -1161/169, -900/169, -261/169, -90/169, -63/169 ], P![]> *];
models[[Integers()|1,6,26,39]] := [* <0, P![ 0, -81/2, 81/2 ], P![]> *];
models[[Integers()|1,3]] := [* <4, P![ -243/8192, 0, -61623/1384448, 0, -8939/692224, 0, -1091/692224, 0, -127/1384448, 0, -3/1384448 ], P![]> *];
models[[Integers()|1,2]] := [* <4, P![ -81/3328, 0, -549/21632, 0, 193/21632, 0, -7/5408, 0, 35/43264, 0, -3/21632 ], P![]> *];
models[[Integers()|1,2,13,26]] := [* <2, P![ -81/3328, -4941/43264, 15633/86528, -5103/43264, 229635/692224, -177147/692224 ], P![]> *];
models[[Integers()|1,6]] := [* <3, P![ 9/3328, 9/416, 63/2704, -441/2704, -495/1352, 9/338, 153/676, -45/169, -63/169 ], P![]> *];
models[[Integers()|1,3,13,39]] := [* <0, P![ -9, 9 ], P![]> *];
models[[Integers()|1,2,3,6]] := [* <2, P![ 0, -9477/512, -44469/512, 140697/1024, -45927/512, 2066715/8192, -1594323/8192 ], P![]> *];
models[[Integers()|1,39]] := [* <0, P![ 9/2, 0, 1/2 ], P![]>, <0, P![ -9, 0, 2 ], P![]>, <0, P![ -9, -18, 9 ], P![]> *];
models[[Integers()|1,26]] := [* <4, P![ -81/3328, -405/1664, -24867/43264, 6723/5408, 17415/2704, 13203/2704, -7857/1352, -243/169, 6885/676, 729/169, -567/169 ], P![]> *];
models[[Integers()|1,6,13,78]] := [* <1, P![ 9/3328, 333/21632, -405/86528, 729/86528, -19683/692224 ], P![]> *];
models[[Integers()|1,2,39,78]] := [* <0, P![ 0, 9/2 ], P![]> *];
models[[Integers()|1,3,26,78]] := [* <2, P![ 0, 81/6656, 2997/43264, -3645/173056, 6561/173056, -177147/1384448 ], P![]> *];
models[[Integers()|1,78]] := [* <3, P![ 9/3328, 0, 37/10816, 0, -5/21632, 0, 1/10816, 0, -3/43264 ], P![]> *];
