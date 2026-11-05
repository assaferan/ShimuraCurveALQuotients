// Subhyperelliptic cover models for X_0(38,7)*
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
// Fibre-product entries (FibreProductCovers.m): built over the star line from the double covers of
// this file.  The W={1,2} entry was rebuilt on 2026-10-13 from the covers as they now stand, and
// checked against the trace formula at five primes over F_p and F_{p^2}.  Coordinates s, z have
// weight 1; y_i has weight half the degree of its own equation.
// Further entries of 2026-10-14 by the same writer, with the coordinate check of the covers above each pair of
 // factors (FibreProductCovers.m), checked at the first three good primes and again at the 4th and 5th.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,133,266]] := [* <1, P![ 64, -192, 272, -128 ], P![]> *];
models[[Integers()|1,14,19,266]] := [* <1, P![ 80, -368, 724, -704, 256 ], P![]> *];
models[[Integers()|1,133]] := [*  *];
models[[Integers()|1,2,7,14]] := [* <4, P![ -5/4, 111/8, -4611/64, 14637/64, -124747/256, 92233/128, -188783/256, 16065/32, -821/4, 38 ], P![]> *];
models[[Integers()|1,266]] := [* <2, P![ 19, 0, 1/4, 0, 1/2, 0, 1/4 ], P![]> *];
models[[Integers()|1,7,38,266]] := [* <0, P![ 5, -8 ], P![]> *];
models[[Integers()|1,14]] := [* <7, "CRV", [ Strings() | "-38*s^9*z + 821/4*s^8*z^2 - 16065/32*s^7*z^3 + 188783/256*s^6*z^4 - 92233/128*s^5*z^5 + 124747/256*s^4*z^6 - 14637/64*s^3*z^7 + 4611/64*s^2*z^8 - 111/8*s*z^9 + 5/4*z^10 + y1^2", "-256*s^4 + 704*s^3*z - 724*s^2*z^2 + 368*s*z^3 - 80*z^4 + y2^2" ]> *];
models[[Integers()|1,19]] := [* <6, "CRV", [ Strings() | "76*s^6 - 249*s^5*z + 5775/16*s^4*z^2 - 2371/8*s^3*z^3 + 2327/16*s^2*z^4 - 81/2*s*z^5 + 5*z^6 + y1^2", "-256*s^4 + 704*s^3*z - 724*s^2*z^2 + 368*s*z^3 - 80*z^4 + y2^2" ]> *];
models[[Integers()|1,2]] := [* <7, "CRV", [ Strings() | "76*s^6 - 249*s^5*z + 5775/16*s^4*z^2 - 2371/8*s^3*z^3 + 2327/16*s^2*z^4 - 81/2*s*z^5 + 5*z^6 + y1^2", "128*s^3*z - 272*s^2*z^2 + 192*s*z^3 - 64*z^4 + y2^2" ]> *];
models[[Integers()|1,2,19,38]] := [* <2, P![ -5, 81/2, -2327/16, 2371/8, -5775/16, 249, -76 ], P![]> *];
models[[Integers()|1,38]] := [* <4, P![ -76, 0, 9/2, 0, -225/512, 0, 141/16384, 0, -1763/4194304, 0, -49/67108864 ], P![]> *];
models[[Integers()|1,7]] := [* <7, "CRV", [ Strings() | "8*s*z - 5*z^2 + y1^2", "-38*s^9*z + 821/4*s^8*z^2 - 16065/32*s^7*z^3 + 188783/256*s^6*z^4 - 92233/128*s^5*z^5 + 124747/256*s^4*z^6 - 14637/64*s^3*z^7 + 4611/64*s^2*z^8 - 111/8*s*z^9 + 5/4*z^10 + y2^2" ]> *];
models[[Integers()|1,14,38,133]] := [* <2, P![ -4, 2, -5, -2, 3, -4, 2 ], P![ 0, 1, 1 ]> *];
