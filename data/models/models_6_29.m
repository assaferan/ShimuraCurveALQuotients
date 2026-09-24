// Subhyperelliptic cover models for X_0(6,29)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
// ✅ REGENERATED 2026-09-09: previously EMPTY cover keys are now filled, unlocked by
// EquationsByRebase (EquationsCovers.m). No flag needed. Verified against the Guo-Yang quotient
// oracle before installing; see data/models/PROVENANCE.md.

models := AssociativeArray();
models[[Integers()|1,6,29,174]] := [* <0, P![ 1/2, 9/16 ], P![]> *];
models[[Integers()|1,3,29,87]] := [* <2, P![ 0, -1/96, 85/2304, 83/1152, -467/2304, -1/4 ], P![]> *];
models[[Integers()|1]] := [* <5, P![ -67/104976, -263/34992, -1253/46656, -11333/419904, -6293/559872, 1679/279936, 30343/1679616, 931/69984, 613/559872, -2975/839808, -379/186624, -1/2187, -1/26244 ], P![]> *];
models[[Integers()|1,3]] := [* <3, P![ -1/96, 0, 85/144, 0, 166/9, 0, -7472/9, 0, -16384 ], P![]> *];
models[[Integers()|1,2]] := [* <2, P![ -2187/131072, 0, 1737/131072, 0, -2425/1179648, 0, -24389/95551488 ], P![]> *];
models[[Integers()|1,6]] := [* <2, P![ 24389/1679616, 0, -1979/17496, 0, 527/2187, 0, -1024/6561 ], P![]> *];
models[[Integers()|1,2,3,6]] := [* <1, P![ -1/96, 85/2304, 83/1152, -467/2304, -1/4 ], P![]> *];
models[[Integers()|1,6,58,87]] := [* <1, P![ -1/768, 7/1152, 5/2304, -1/36 ], P![]> *];
models[[Integers()|1,29]] := [* <3, P![ -24389/1889568, 0, 119381/944784, 0, -8174/19683, 0, 33488/59049, 0, -16384/59049 ], P![]> *];
models[[Integers()|1,3,58,174]] := [* <0, P![ 0, 1/16 ], P![]> *];
models[[Integers()|1,2,29,58]] := [* <1, P![ 0, -1/768, 7/1152, 5/2304, -1/36 ], P![]> *];
models[[Integers()|1,58]] := [* <2, P![ -1/768, 0, 7/72, 0, 5/9, 0, -1024/9 ], P![]> *];
models[[Integers()|1,87]] := [* <3, P![ -2187/1048576, 0, 495/262144, 0, -2081/4718592, 0, -641/191102976, 0, 24389/6879707136 ], P![]> *];
models[[Integers()|1,2,87,174]] := [* <0, P![ 0, 8, 9 ], P![]> *];
models[[Integers()|1,174]] := [* <0, P![ -1/18, 0, 1/9 ], P![]>, <0, P![ 1/32, 1/32, -1/16 ], P![]>, <0, P![ 1/2, 0, 9 ], P![]> *];
