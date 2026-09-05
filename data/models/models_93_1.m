// Subhyperelliptic cover models for X_0(93,1)*
//
// PROVENANCE. Produced 2026-09-05 on lovelace in 50927 s (14.1 h) with NO non-default flags:
//     NORMALIZ_BIN=... magma -b D_s:=93 N_s:=1 OUTDIR:=... genmodels.m < /dev/null
// It regenerates from committed `main` -- but only SINCE the vx fix (`n_oo`, BorcherdsForms.m
// ~:771/:787/:858). Before that fix this base died in the odd-D `oo`-expansion block, which is why
// it had no model until now. The tree that produced it differs from `main` only by `main` carrying
// MORE diagnostics (BFPROGRESS/COVPROGRESS/the LogSum runaway-coefficient message) and the
// EtaQuotient speedups; none of those touch the computation.
//
// VALIDATED AGAINST GUO-YANG -- and it settles a typo in their table. Their entry reads
//     y^2 = (3s^3 - 7s^2 - 3t - 1)(3s^3 + s^2 - 3s - 9),   x^2 = -4s^2 - 6s - 9
// where `t` appears nowhere else in the row. Testing the plausible repairs against our genus-2
// quotient decides it: only `-3t -> -3s` gives an isomorphic curve; `-3s^2`, `-3`, and dropping
// the term all fail. So the published `-3t` is a typo for `-3s`.
// With that reading the whole V_4 diagram matches (tests/GuoYangEquations.m):
//     GY y-quotient   y^2 = AB    genus 2  ==  our [1,93]     IsIsomorphic
//     GY product      y^2 = ABC   genus 3  ==  our [1,3]      IsIsomorphic
//     GY conic        x^2 = C     genus 0  ==  our [1,31]     same class in Q*/Q*^2
// ⚠ SCOPE: that is the QUOTIENT diagram, not the genus-5 full curve. A direct IsIsomorphic on the
// W={1} CRV pair is the 10h+ regime (cf. 26_3) and is not run. Three of the four cover keys are
// pinned exactly and the fourth is their fibre product, so this is strong but not a full-curve
// proof. ModelChecks validates the file structurally.
//
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,3]] := [* <3, P![ -7/729, 46/729, -53/243, 134/243, -758/729, 332/243, -1007/729, 68/81, -16/81 ], P![]> *];
models[[Integers()|1,31]] := [* <0, P![ -63, 36, -144 ], P![]> *];
models[[Integers()|1,93]] := [* <2, P![ 1/9, -2/3, 17/9, -34/9, 50/9, -4, 1 ], P![]> *];
models[[Integers()|1]] := [* <5, "CRV", [ Strings() | "y^2 - s^6 + 4*s^5*z - 50/9*s^4*z^2 + 34/9*s^3*z^3 - 17/9*s^2*z^4 + 2/3*s*z^5 - 1/9*z^6", "x^2 + 144*s^2 - 36*s*z + 63*z^2" ]> *];
