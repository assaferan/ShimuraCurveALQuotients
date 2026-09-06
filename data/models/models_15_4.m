// Subhyperelliptic cover models for X_0(15,4)*
//
// ⚠⚠ THIS FILE IS NOT A PRODUCT OF THIS PIPELINE. A FOURTH PROVENANCE CATEGORY.
// X_0^15(4) is OUTSIDE the Guo-Yang method by the authors' own statement -- see their published
// Remark 39 (Compositio Math. 153 (2017) 1-40; ABSENT from arXiv v1):
//   "there is a curve, namely, X = X_0^15(4), whose equation is not obtained using our method.
//    This is because the normalizer of the Eichler order in this case is larger than the
//    Atkin-Lehner group. For this special curve, we use the result of Tu [Tu14]."
// Since N^+_B(O) strictly contains W_{15,4}, the star quotient our pipeline forms is the WRONG
// OBJECT. Do NOT try to make genmodels.m produce this base; it also trips
// `assert IsSquarefree(N)` and, behind NONSQFREE=1, the Hall-divisor star check.
//
// HOW IT WAS OBTAINED, and WHAT IS ACTUALLY VERIFIED HERE (2026-09-06).
// INPUT, transcribed: Tu (Pacific J. Math. 269 (2014), Lemma 13), as QUOTED by Guo-Yang -- a
// Hauptmodul t4 on X/<w_3,w_5> taking values +-1/sqrt(-3), +-sqrt(-15)/5, (+-1+-sqrt(-15))/8 at
// CM discriminants -12, -15, -60. Tu's paper itself was not needed: GY quote all of it we use.
//
// DERIVED AND CHECKED IN CODE (not taken on faith):
//   * the polynomials are FORCED -- they are exactly the minimal polynomials of those CM values
//     ((4x^2-x+1)(4x^2+x+1) from -60, (5x^2+3) from -15, (3x^2+1) from -12), up to the rational
//     scalars 80 and 3;
//   * the ramification assignment (X/w_3 over -15 and -60; X/w_15 over -12) then gives the two
//     equations below;
//   * our own Shimura genus formula independently reproduces GY's quotient structure:
//     W={1} g=5, W={1,3} g=2, W={1,15} g=0, W={1,3,5,15} g=0;
//   * the CRV pair below has genus 5, computed.
//
// ⚠ WHAT IS INDEPENDENTLY CONFIRMED vs WHAT RESTS ON GUO-YANG:
//   * a = -1 (the y-side constant) is CONFIRMED BY US and the check DISCRIMINATES: Eichler-Selberg
//     point counts give 16 checks / 0 failures at a = -1, while a = 1, 2, -2, 5, -5 each FAIL 2-3.
//     That check is independent of the Borcherds/Schofer path.
//   * b = -1 (the conic constant) is NOT confirmed by us. Point counts are VACUOUS here -- every
//     b passes 16/16, because the conic is genus 0 and the counts do not constrain its twist.
//     `b < 0` is forced (the conic must have no real points: X_0^15(4) is not hyperelliptic over
//     R, per Ogg and GY), but b = -1, -2, -3, -5, -6 are ALL admissible on that criterion.
//     b = -1 rests on Guo-Yang's Schofer/CM argument, which we have NOT reproduced.
//     ⇒ If you need b, close it by computing the CM values of t_2 on X_0^15(2)/<w_3,w_5> with our
//       Schofer code and checking GY's relation t_2 = (5t_4^2+2t_4+1)/(7t_4^2-2t_4+3).
//
// VALIDATION: VerifyModelSet(15, 4) -> 16 checks, 0 failures.
//
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,3]] := [* <2, -(4*x^2 - x + 1)*(4*x^2 + x + 1)*(5*x^2 + 3), P![]> *];
models[[Integers()|1,15]] := [* <0, -(3*x^2 + 1), P![]> *];
models[[Integers()|1]] := [* <5, "CRV", [ Strings() |
  "y^2 + (4*x^2 - x*s + s^2)*(4*x^2 + x*s + s^2)*(5*x^2 + 3*s^2)",
  "z^2 + 3*x^2 + s^2" ]> *];
