// Subhyperelliptic cover models for X_0(111,1)*
//
// PROVENANCE. Produced 2026-09-06/07 on lovelace, DEFAULT flags (BFPROGRESS/COVPROGRESS/M0PROGRESS
// are diagnostics only):
//     NORMALIZ_BIN=... magma -b D_s:=111 N_s:=1 OUTDIR:=... genmodels.m < /dev/null
// AllEquationsAboveCovers took 72617 s (20.2 h); 4 cover-keys, all populated.
// ⚠ Pinnable to a single commit -- lovelace's clone was at f87b0ae with a CLEAN working tree for
// the whole run. (Contrast models_10_61.m and models_14_43.m, which cannot be, because that clone
// was pulled while they were running.)
// It required the vx fix, like 93_1: before `n_oo` this base died in the odd-D oo-expansion block.
//
// VALIDATED AGAINST GUO-YANG BY EXACT FULL-CURVE ISOMORPHISM -- the strongest form, and cheap here
// because the W={1} entry is hyperelliptic rather than a CRV pair (0.05 s; contrast the CRV pairs,
// where IsIsomorphic runs for hours):
//     GY:  y^2 = -(19x^8 - 44x^7 - 16x^6 + 55x^5 + 37x^4 - 55x^3 - 16x^2 + 44x + 19)
//                 (x^8 - 3x^5 - x^4 + 3x^3 + 1)                      degree 16, genus 7
//     ours: W={1}, degree 16, genus 7    -->    IsIsomorphic TRUE
// (Guo-Yang, Compositio Math. 153 (2017), "Equations of level one" table.)
// Checked in CI by tests/GuoYangEquations.m.
//
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,3]] := [* <4, P![ -10097379, 41452398, -86388687, 119003418, -118111122, 87589350, -48925377, 20234124, -5918022, 1102248, -98415 ], P![]> *];
models[[Integers()|1,37]] := [* <3, P![ -124659, 428652, -711504, 756702, -558414, 288684, -101331, 21870, -2187 ], P![]> *];
models[[Integers()|1]] := [* <7, P![ -124659, 288684, 104976, 13122, -984150, -616734, 1948617, 1115370, -2171691, -1115370, 1948617, 616734, -984150, -13122, 104976, -288684, -124659 ], P![]> *];
models[[Integers()|1,111]] := [* <0, P![ 81, -54, 45 ], P![]> *];
