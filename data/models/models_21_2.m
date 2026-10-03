// Subhyperelliptic cover models for X_0(21,2)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
// Non-hyperelliptic covers are stored as <genus, "CRV", [defining-poly strings]>.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[ 1 ]] := [* <3, "CRV", [ Strings() | "y^2 + 7/256*s^4 + 25/32*s^2*z^2 + 7/16*z^4", "x^2 + 3*s^2 + 3*z^2" ]> *];
models[[ 1, 2 ]] := [* <1, P![ -4/9, -16/9, -23/9, -14/9, -7/9 ], P![]>, <1, P![ -7/16, 0, -25/32, 0, -7/256 ], P![]> *];
models[[ 1, 42 ]] := [* <1, P![ -7/16, 0, -1/32, 0, 9/256 ], P![]> *];
models[[ 1, 21 ]] := [* <1, P![ -7/256, 0, 31/128, 0, 9/256 ], P![]> *];
models[[ 1, 7 ]] := [* <0, P![ -3, 0, -1 ], P![]>, <0, P![ -3, 0, -3 ], P![]>, <0, P![ -3, 0, -1 ], P![]> *];
models[[ 1, 2, 7, 14 ]] := [* <0, P![ 0, 81/4, -81/4 ], P![]> *];
models[[ 1, 6, 7, 42 ]] := [* <0, P![ 0, -3 ], P![]> *];
models[[ 1, 3, 7, 21 ]] := [* <0, P![ -3, 3 ], P![]> *];
models[[ 1, 2, 3, 6 ]] := [* <1, P![ -4/9, 32/9, -71/9, 28/9 ], P![]> *];
models[[ 1, 2, 21, 42 ]] := [* <0, P![ -7/16, 3/32, 81/256 ], P![]> *];
models[[ 1, 3, 14, 42 ]] := [* <1, P![ 0, 7/4, -1/8, -9/64 ], P![]> *];
models[[ 1, 6, 14, 21 ]] := [* <1, P![ 0, 7/64, 31/32, -9/64 ], P![]> *];
// The three genus-2 quotients, by the Atkin-Lehner involutions without fixed points.  Each is
// isomorphic to the quotient of Guo-Yang's equations for X_0^21(2) (Compositio 153 (2017),
// Appendix A, Table A.2); tests/X0_21_2.m checks that.
models[[ 1, 3 ]] := [* <2, P![ 189/1024, 0, -1611/1024, 0, -801/1024, 0, -81/1024 ], P![]> *];
models[[ 1, 6 ]] := [* <2, P![ 189/16, 0, 153/32, 0, -171/256, 0, -81/256 ], P![]> *];
models[[ 1, 14 ]] := [* <2, P![ 189/64, 0, 1053/128, 0, 5589/1024, 0, 189/1024 ], P![]> *];
