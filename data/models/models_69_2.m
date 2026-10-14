// Subhyperelliptic cover models for X_0(69,2)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
// Fibre-product entries (FibreProductCovers.m, 2026-10-14): built over the star line from the double covers of
// this file by tools/writefp.m with the coordinate check of the three covers above each pair of factors;
// each checked against the trace formula at the first three good primes, F_p and F_{p^2}, and again at the
// 4th and 5th.  Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,23,46]] := [* <2, P![ 0, 3/529, 29/2116, 29/1058, 9/2116, -11/529, -16/529 ], P![]> *];
models[[Integers()|1,3,46,138]] := [* <1, P![ -4/3, -8/9, 4/9, 16/9 ], P![]> *];
models[[Integers()|1,69]] := [* <4, "CRV", [ Strings() | "-1/4*s^4 - 1/18*s^3*z + 7/36*s^2*z^2 + 1/9*s*z^3 + y1^2", "4/4761*s^5*z + 95/171396*s^4*z^2 - 31/85698*s^3*z^3 - 35/57132*s^2*z^4 - 14/42849*s*z^5 - 4/42849*z^6 + y2^2" ]> *];
models[[Integers()|1,6,23,138]] := [* <2, P![ 0, 1/3, 53/36, 49/18, 89/36, 1 ], P![]> *];
models[[Integers()|1,2,69,138]] := [* <1, P![ 0, -1/9, -7/36, 1/18, 1/4 ], P![]> *];
models[[Integers()|1,46]] := [* <4, "CRV", [ Strings() | "16/529*s^6 + 11/529*s^5*z - 9/2116*s^4*z^2 - 29/1058*s^3*z^3 - 29/2116*s^2*z^4 - 3/529*s*z^5 + y1^2", "-16/9*s^3*z - 4/9*s^2*z^2 + 8/9*s*z^3 + 4/3*z^4 + y2^2" ]> *];
models[[Integers()|1,3]] := [* <5, "CRV", [ Strings() | "16/4761*s^6 + 419/42849*s^5*z + 2335/171396*s^4*z^2 + 35/3174*s^3*z^3 + 955/171396*s^2*z^4 + 74/42849*s*z^5 + 4/14283*z^6 + y1^2", "4/4761*s^5*z + 95/171396*s^4*z^2 - 31/85698*s^3*z^3 - 35/57132*s^2*z^4 - 14/42849*s*z^5 - 4/42849*z^6 + y2^2" ]> *];
models[[Integers()|1,2]] := [* <5, "CRV", [ Strings() | "16/529*s^6 + 11/529*s^5*z - 9/2116*s^4*z^2 - 29/1058*s^3*z^3 - 29/2116*s^2*z^4 - 3/529*s*z^5 + y1^2", "16/4761*s^6 + 419/42849*s^5*z + 2335/171396*s^4*z^2 + 35/3174*s^3*z^3 + 955/171396*s^2*z^4 + 74/42849*s*z^5 + 4/14283*z^6 + y2^2" ]> *];
models[[Integers()|1,6]] := [* <5, "CRV", [ Strings() | "16/4761*s^6 + 419/42849*s^5*z + 2335/171396*s^4*z^2 + 35/3174*s^3*z^3 + 955/171396*s^2*z^4 + 74/42849*s*z^5 + 4/14283*z^6 + y1^2", "-s^5*z - 89/36*s^4*z^2 - 49/18*s^3*z^3 - 53/36*s^2*z^4 - 1/3*s*z^5 + y2^2" ]> *];
models[[Integers()|1,2,3,6]] := [* <2, P![ -4/14283, -74/42849, -955/171396, -35/3174, -2335/171396, -419/42849, -16/4761 ], P![]> *];
models[[Integers()|1,138]] := [* <4, "CRV", [ Strings() | "-16/9*s^3*z - 4/9*s^2*z^2 + 8/9*s*z^3 + 4/3*z^4 + y1^2", "-s^5*z - 89/36*s^4*z^2 - 49/18*s^3*z^3 - 53/36*s^2*z^4 - 1/3*s*z^5 + y2^2" ]> *];
models[[Integers()|1,3,23,69]] := [* <2, P![ 4/42849, 14/42849, 35/57132, 31/85698, -95/171396, -4/4761 ], P![]> *];
models[[Integers()|1,23]] := [* <6, "CRV", [ Strings() | "-s^5*z - 89/36*s^4*z^2 - 49/18*s^3*z^3 - 53/36*s^2*z^4 - 1/3*s*z^5 + y1^2", "4/4761*s^5*z + 95/171396*s^4*z^2 - 31/85698*s^3*z^3 - 35/57132*s^2*z^4 - 14/42849*s*z^5 - 4/42849*z^6 + y2^2" ]> *];
models[[Integers()|1,6,46,69]] := [* <1, P![ 0, -4/4761, -7/4761, -16/4761 ], P![]> *];
