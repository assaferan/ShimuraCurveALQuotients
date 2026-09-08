// ⚠⚠ KNOWN DEFECT (found 2026-09-08): THE [1,29] ENTRY BELOW IS WRONG. It is NOT the quotient of
// X_0^87(1) by w_29. Evidence: our [1] full curve IS isomorphic to Guo-Yang's published equation,
// and our [1,3] and [1,87] both match the corresponding quotients derived from their curve and
// involutions -- but [1,29] disagrees with theirs in POINT COUNT at 8 of 10 small primes
// (ours/theirs 10/4 at p=7, 8/20 at 11, 12/10 at 13, 16/14 at 17, 26/22 at 19, 21/18 at 23,
// 44/30 at 37), which REFUTES isomorphism rather than merely failing to establish it. Both are
// genus 3, so genus does not catch it.
// ⇒ This was invisible until tests/GuoYangQuotientOracle.m, because Guo-Yang publish only the FULL
// curve for 87_1 and tests/GuoYangEquations.m could therefore compare only that. The defect is
// pinned there as expected-to-mismatch; fix the entry and the test will tell you to delist it.
// ⚠ NOT diagnosed. 87_1 already had an open question about multiple bases per cover.

// Subhyperelliptic cover models for X_0(87,1)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,29]] := [* <3, P![ -5, 14, 23, -81, -36, 93, 70, 18, 3 ], P![ 1, 0, 1, 1 ]> *];
models[[Integers()|1,3]] := [* <2, P![ -129140163/3444736, 1190959281/1722368, -15635525661/3444736, 10581521751/861184, -34231709133/3444736, -8451506223/1722368, -10460353203/3444736 ], P![]> *];
models[[Integers()|1,87]] := [* <0, P![ 0, -27 ], P![]> *];
models[[Integers()|1]] := [* <5, P![ -129140163/3444736, 0, -44109603/1722368, 0, -21447909/3444736, 0, -537597/861184, 0, -64413/3444736, 0, 589/1722368, 0, -27/3444736 ], P![]> *];
