// Subhyperelliptic cover models for X_0(26,3)*
//
// ⚠ PRODUCED WITH CMNONCOPRIME=1 -- this file does NOT regenerate under default settings.
//     CMNONCOPRIME=1 NORMALIZ_BIN=... magma -b D_s:=26 N_s:=3 OUTDIR:=... genmodels.m < /dev/null
// 189 s, 15 cover-keys (2 empty: [1,2] and [1,13]).
// Under the default coprime-to-level CM filter the base dies at "Could not find enough points":
// the filter admits only 3 of Guo-Yang's 14 discriminants (only -8, -11, -20 are coprime to N=3)
// against demand 9. With the filter off the pool is 21.
//
// VALIDATED AGAINST GUO-YANG -- the full V_4 diagram, exactly:
//     GY y-quotient  y^2 = x^6 - 2x^4 + 9x^2 + 8   genus 2  ==  our [1,78]   IsIsomorphic
//     GY product     y^2 * (-8x^2-3)               genus 3  ==  our [1,3]    IsIsomorphic
//     GY conic       z^2 = -8x^2 - 3               genus 0  ==  our [1,26] entry 2, which is
//                                                              -8x^2-3 COEFFICIENT FOR COEFFICIENT
// (Guo-Yang, Compositio 153 (2017), "Equations of level greater than one" table.)
//
// ⚠⚠ THE KNOWN-BAD CM VALUES DID NOT CORRUPT THE MODEL, and that is the interesting part.
// This base has a recorded s <-> s~ SWAP at discriminants -267 and -708 -- Guo-Yang's `s` sits in
// our `s~` row there -- and that misbehaviour is the entire justification for the coprime filter
// existing. Yet with the filter OFF, so those two bad rows are admitted, the model still comes out
// isomorphic to Guo-Yang's on every quotient. So the two wrong values are not load-bearing for the
// covers; the filter was blunter than the defect it guards against.
// ⇒ See data/models/PROVENANCE.md: this is evidence that the coprime guard may not be needed.
//
// ⚠ NOT a defence of CMNONCOPRIME in general. The `p | gcd(d,N)` local factor still has no live
// implementation (`kappaminuszero` is dead code), so there is no theoretical guarantee; what makes
// THIS file trustworthy is the published equation, not the flag. Same bar as 39_2 and 14_3.
//
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,13]] := [*  *];
models[[Integers()|1]] := [* <5, "CRV", [ Strings() | "y^2 + 1/32768*s^6 + 25/32768*s^4*z^2 + 699/32768*s^2*z^4 - 2197/32768*z^6", "x^2 + 1/8*s^2 + 3/8*z^2" ]> *];
models[[Integers()|1,6,26,39]] := [* <0, P![ -11, 16 ], P![]> *];
models[[Integers()|1,3]] := [* <3, P![ -3/8, 0, -91/64, 0, -33/32, 0, 13/64, 0, -1/8 ], P![]> *];
models[[Integers()|1,2]] := [*  *];
models[[Integers()|1,2,13,26]] := [* <0, P![ -11, 38, -32 ], P![]> *];
models[[Integers()|1,6]] := [* <2, P![ 2197/32768, 0, -699/32768, 0, -25/32768, 0, -1/32768 ], P![]> *];
models[[Integers()|1,3,13,39]] := [* <2, P![ -11/4, 49/4, -291/16, 47/4, -27/4, 4 ], P![]> *];
models[[Integers()|1,2,3,6]] := [* <1, P![ -11/4, 27/4, -75/16, 19/8, -2 ], P![]> *];
models[[Integers()|1,39]] := [* <3, P![ -6591/262144, 0, -25/65536, 0, 387/131072, 0, 7/65536, 0, 1/262144 ], P![]> *];
models[[Integers()|1,26]] := [* <0, P![ -3/64, 0, -3/8 ], P![]>, <0, P![ -3, 0, -8 ], P![]>, <0, P![ -3/8, 0, -1/8 ], P![]> *];
models[[Integers()|1,6,13,78]] := [* <1, P![ 1/4, -1/4, 1/16, -1/8 ], P![]> *];
models[[Integers()|1,2,39,78]] := [* <1, P![ 1/4, -3/4, 9/16, -1/4, 1/4 ], P![]> *];
models[[Integers()|1,3,26,78]] := [* <0, P![ 1, -2 ], P![]> *];
models[[Integers()|1,78]] := [* <2, P![ 1/8, 0, 9/64, 0, -1/32, 0, 1/64 ], P![]> *];
