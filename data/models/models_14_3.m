// Subhyperelliptic cover models for X_0(14,3)*
//
// ⚠ PRODUCED WITH CMNONCOPRIME=1 -- this file does NOT regenerate under default settings.
// With the coprime-to-level CM filter on, 14_3's covers are under-determined and W={1} comes out
// EMPTY (the previous committed file had 6 keys, 3 of them empty); with the filter off it builds
// 16 keys, 0 empty. Regenerate with:
//     CMNONCOPRIME=1 NORMALIZ_BIN=... magma -b D_s:=14 N_s:=3 OUTDIR:=... genmodels.m
// tests/_offline/ModelRegen.m lists 14_3 in MR_KNOWN_DRIFT for this reason.
//
// VERIFIED AGAINST THE LITERATURE. Its W={1} genus-3 curve is IsIsomorphic to Guo-Yang's
//     z^2 = -9x^2-2 ,  y^2 = -7x^4+22x^2+1
// (exact, 6816 s -- see tests/_offline/GuoYangCurve_14_3.m, which is offline BECAUSE of that
// cost). The flag is not trusted; the oracle is.
//
// ⚠ A TRAP, recorded because it nearly discarded this result. Our W={1} is presented over the
// Klein 4-group {1,2,7,14}; Guo-Yang's over {1,3,14,42}. DIFFERENT V4s of the same curve. So
// comparing our presentation's quotients against theirs gives false on both, and the two degree-4
// branch forms even have different Galois groups (C2^2 vs D4) -- from which I briefly, and
// wrongly, concluded the curves differed. W={1} here has SEVEN involution quotients (genera
// 1,2,2,2,0,1,1) and several valid V4s, so "genus distinguishes the quotients, hence they must
// correspond" is FALSE. Asking instead whether each of THEIR quotients matches ANY of ours:
//     GY genus-1 quotient ~ our [1,42]      GY genus-2 quotient ~ our [1,3]
// which is consistent, and the exact test then confirmed it.
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,7,14]] := [* <0, P![ -11, 104, -128 ], P![]> *];
models[[Integers()|1,42]] := [* <1, P![ -28, 0, 88, 0, 4 ], P![]> *];
models[[Integers()|1,14]] := [* <0, P![ -9, 0, -2 ], P![]>, <0, P![ -9/2, 0, -1/2 ], P![]>, <0, P![ -3, -2, -3 ], P![]> *];
models[[Integers()|1,3,7,21]] := [* <1, P![ -176, 2368, -3776, 1024 ], P![]> *];
models[[Integers()|1]] := [* <3, "CRV", [ Strings() | "y^2 + 256*s^4 + 1792/3*s^3*z + 6592/9*s^2*z^2 + 1792/3*s*z^3 + 256*z^4", "x^2 + 3*s^2 + 2*s*z + 3*z^2" ]> *];
models[[Integers()|1,3]] := [* <2, P![ 63, 0, -184, 0, -53, 0, -2 ], P![]> *];
models[[Integers()|1,2]] := [* <1, P![ -851, -1370, -883, -268, -32 ], P![]>, <1, P![ -256, -1792/3, -6592/9, -1792/3, -256 ], P![]> *];
models[[Integers()|1,7]] := [* <2, P![ 768, 2304, 4160, 45440/9, 4160, 2304, 768 ], P![]> *];
models[[Integers()|1,6]] := [* <2, P![ 3087/2, 0, 577/2, 0, 17/2, 0, -1/2 ], P![]> *];
models[[Integers()|1,2,3,6]] := [* <1, P![ -11, 236, -1420, 1952, -512 ], P![]> *];
models[[Integers()|1,2,21,42]] := [* <0, P![ 64, -768, 256 ], P![]> *];
models[[Integers()|1,6,14,21]] := [* <0, P![ -11, 16 ], P![]> *];
models[[Integers()|1,3,14,42]] := [* <0, P![ 1, -8 ], P![]> *];
models[[Integers()|1,21]] := [* <1, P![ -343, 0, -26, 0, 1 ], P![]> *];
models[[Integers()|1,6,7,42]] := [* <1, P![ 64, -1280, 6400, -2048 ], P![]> *];
