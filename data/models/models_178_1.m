// Subhyperelliptic cover models for X_0(178,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2]] := [* <4, P![ -110592, 2129920/3, -17779712/9, 29249536/9, -294268928/81, 709578752/243, -139225088/81, 1604882432/2187, -478117888/2187, 810070016/19683, -228502528/59049 ], P![]> *];
models[[Integers()|1]] := [*  *];
models[[Integers()|1,89]] := [* <1, P![ -432, 1408/3, -3616/9, 1408/9, -2608/81 ], P![]> *];
models[[Integers()|1,178]] := [* <2, P![ 256, -4096/3, 25664/9, -82688/27, 155264/81, -168704/243, 87616/729 ], P![]> *];
