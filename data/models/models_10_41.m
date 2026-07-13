// Subhyperelliptic cover models for X_0(10,41)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,10]] := [*  *];
models[[Integers()|1,41]] := [*  *];
models[[Integers()|1,2,41,82]] := [* <3, P![ 0, -512/5, 31808/25, -150528/25, 322944/25, -234496/25, -225728/25, 92672/5, -8192 ], P![]> *];
models[[Integers()|1,2]] := [*  *];
models[[Integers()|1,2,5,10]] := [* <3, P![ -20480, 272384, -1344960, 2701312, -255616, -6617600, 6337600, 3392000, -5120000 ], P![]> *];
models[[Integers()|1,5]] := [*  *];
models[[Integers()|1,205]] := [*  *];
models[[Integers()|1,2,205,410]] := [* <1, P![ 0, 1/2, -23/16, -5/8, 25/16 ], P![]> *];
models[[Integers()|1,10,82,205]] := [* <2, P![ 0, -128/25, 144/5, -864/25, -688/25, 256/5 ], P![]> *];
models[[Integers()|1,82]] := [*  *];
models[[Integers()|1,5,82,410]] := [* <1, P![ 5/16, -17/8, 69/16, -5/2 ], P![]> *];
models[[Integers()|1,10,41,410]] := [* <2, P![ 0, 125/128, -7675/1024, 7875/512, 3125/1024, -3125/128 ], P![]> *];
models[[Integers()|1,5,41,205]] := [* <3, P![ -4096/25, 34816/25, -88768/25, 15872/25, 210816/25, -35328/5, -4800, 5120 ], P![]> *];
models[[Integers()|1,410]] := [*  *];
