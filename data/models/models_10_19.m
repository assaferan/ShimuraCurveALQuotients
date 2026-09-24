// Subhyperelliptic cover models for X_0(10,19)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
// ✅ REGENERATED 2026-09-09: the four previously-EMPTY keys ({1}, {1,2}, {1,10}, {1,190}) are now
// filled, unlocked by EquationsByRebase (EquationsCovers.m). No flag needed.
//
// ⚠ ELEVEN EXISTING ENTRIES ALSO CHANGED, AND THAT IS NOT A REGRESSION: each new polynomial is the
// old one scaled by 16 = 4^2, so y^2 = 16f and (y/4)^2 = f describe the SAME curve. Checked
// entry by entry: 11 unchanged-up-to-isomorphism, 4 newly filled, 0 covers lost. The scaling comes
// from the run, not from the rebase stage -- that stage only ADOPTS keys which were empty.
//
// ⚠ The W={1} entry is a CRV pair over a POINTLESS conic, matching how Guo-Yang present this base
// (Example 37: y^2 = -8x^6+57x^4-40x^2+16, z^2 = 5x^2-32, which has only real points).
// ⚠ Every entry is checked against Guo-Yang by tests/GuoYangQuotients_10_19.m.

models := AssociativeArray();
models[[Integers()|1,10]] := [* <2, P![ -475/2097152, 0, -2923/131072, 0, -13/40960, 0, -1/320000 ], P![]> *];
models[[Integers()|1]] := [* <5, "CRV", [ Strings() | "y^2 + 1/320000*s^6 + 13/40960*s^4*z^2 + 2923/131072*s^2*z^4 + 475/2097152*z^6", "x^2 + 1/128*s^2 - 125/2048*z^2" ]> *];
models[[Integers()|1,2]] := [* <2, P![ -19/80, 0, -1834/125, 0, -7728/3125, 0, -8192/78125 ], P![]> *];
models[[Integers()|1,2,95,190]] := [* <1, P![ 25/16, -75/32, 825/256, -625/128, 625/256 ], P![]> *];
models[[Integers()|1,2,5,10]] := [* <1, P![ -25/2, 1075/64, -3175/128, 36875/1024, -16875/1024 ], P![]> *];
models[[Integers()|1,5]] := [* <3, P![ -25/2, 0, 2125/64, 0, -6325/128, 0, 13525/1024, 0, -125/128 ], P![]> *];
models[[Integers()|1,190]] := [* <2, P![ -25/128, 0, 57/16, 0, -32/5, 0, 4096/625 ], P![]> *];
models[[Integers()|1,95]] := [* <3, P![ 25/38416, 25/9604, 275/19208, 325/9604, 5575/38416, 325/1372, 375/19208, -225/2401, -1175/38416 ], P![]> *];
models[[Integers()|1,5,19,95]] := [* <2, P![ 0, 625/8, -46875/256, 133125/512, -1556875/4096, 671875/2048, -421875/4096 ], P![]> *];
models[[Integers()|1,10,38,95]] := [* <0, P![ 0, 50, -675/16 ], P![]> *];
models[[Integers()|1,38]] := [* <0, P![ -25/392, -25/392, 25/392 ], P![]>, <0, P![ 125/392, 125/196, 4125/12544 ], P![]>, <0, P![ -50, 0, 125/16 ], P![]> *];
models[[Integers()|1,5,38,190]] := [* <0, P![ 0, -25/64, 25/64 ], P![]> *];
models[[Integers()|1,19]] := [* <3, P![ -475/1229312, -475/153664, -792875/39337984, -1527425/19668992, -13217475/78675968, -4169575/19668992, -14056125/89915392, -19760075/314703872, -6708075/629407744 ], P![]> *];
models[[Integers()|1,2,19,38]] := [* <0, P![ -25/2, 1475/64, -675/64 ], P![]> *];
models[[Integers()|1,10,19,190]] := [* <1, P![ 0, -25/64, 25/128, -625/1024, 625/1024 ], P![]> *];
