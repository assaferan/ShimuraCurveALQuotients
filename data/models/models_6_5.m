// Subhyperelliptic cover models for X_0(6,5)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,10]] := [* <0, P![ -3/1024, 0, -3/64 ], P![]>, <0, P![ -3, 0, -64 ], P![]>, <0, P![ -3/64, 0, -1/64 ], P![]> *];
models[[Integers()|1,15]] := [* <1, P![ -125/12288, 0, -61/18432, 0, 1/36864 ], P![]> *];
models[[Integers()|1,3,10,30]] := [* <0, P![ 0, 1/4 ], P![]> *];
models[[Integers()|1]] := [* <1, P![ -64375/16384, -15625/2048, 29375/8192, 15625/2048, -64375/16384 ], P![]>, <1, P![ -1216, 3712, -2240, -4096, -1024 ], P![]> *];
models[[Integers()|1,3]] := [* <1, P![ -8/3, 0, -524/9, 0, -256/9 ], P![]> *];
models[[Integers()|1,2]] := [* <0, P![ -1/96, 0, -125/768 ], P![]>, <0, P![ -1375/256, -375/128, -375/256 ], P![]>, <0, P![ -24, 0, -1000 ], P![]> *];
models[[Integers()|1,2,5,10]] := [* <0, P![ 0, -3, -16 ], P![]> *];
models[[Integers()|1,5]] := [* <1, P![ -1000, 0, 4048, 0, -4096 ], P![]> *];
models[[Integers()|1,6,10,15]] := [* <0, P![ -3, -16 ], P![]> *];
models[[Integers()|1,2,15,30]] := [* <0, P![ 0, 1/18, 1/144 ], P![]> *];
models[[Integers()|1,6]] := [* <0, P![ 125, 0, -256 ], P![]>, <0, P![ 2125/4096, 125/128, 125/256 ], P![]>, <0, P![ 125/256, 0, -1/256 ], P![]> *];
models[[Integers()|1,2,3,6]] := [* <0, P![ -8/3, -131/9, -16/9 ], P![]> *];
models[[Integers()|1,30]] := [* <0, P![ -2, 0, 4 ], P![]>, <0, P![ 1/2, 0, -1/2 ], P![]>, <0, P![ 1/2, 0, 1/4 ], P![]> *];
models[[Integers()|1,5,6,30]] := [* <0, P![ 1/2, 1/16 ], P![]> *];
// X_0(6,5)/<w_3, w_5>, genus 1.  It is the quotient of Gonzalez-Rotger's model of X_0(6,5)
// (J. Math. Soc. Japan 58 (2006), Table 1, p. 8: y^2 = -x^4 + 61x^2 - 1024, with w_30 = (x,-y),
// w_2 = (-x,y), w_6 = (32/x, 32y/x^2)) first by w_15 = w_2 w_30 and then by w_3 = w_2 w_6, the
// elliptic curve 30a6; tests/X0_6_5.m checks the isomorphism.
models[[Integers()|1,3,5,15]] := [* <1, P![ 0, -8/3, -131/9, -16/9 ], P![]> *];
