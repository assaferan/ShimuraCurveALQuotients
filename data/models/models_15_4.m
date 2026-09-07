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
//   * b = -1 (the conic constant) is ALSO CONFIRMED BY US -- 2026-09-06, see below. It was open
//     for one commit (59de486) and is now closed.
//
// HOW b WAS CLOSED, and why the obvious check could NOT do it.
// `VerifyModelSet` is VACUOUS for b: every value passes 16/16, because the conic is genus 0 and
// point counts do not constrain its twist. Two weaker criteria narrowed it:
//   (a) Guo-Yang's own field-of-definition criterion, applied in the RING CLASS FIELD
//       H = Q(sqrt(-3), sqrt(5)) (NOT in K = Q(sqrt(-15)) -- doing it in K gives an answer that
//       CONTRADICTS Guo-Yang, because the CM points of discriminant -15/-60 are defined over H,
//       h(-15) = 2; that was a wrong-object error caught by the contradiction). In H the criterion
//       accepts the square class {-1, 3, -5, 15} and REJECTS 7 of 12 tested values.
//   (b) no real points (X_0^15(4) is not hyperelliptic over R) forces b < 0, killing 3 and 15.
// That left exactly {-1, -5} -- genuinely different conics, (-3,-1) ramified at {3, oo} versus
// (-15,-5) ramified at {5, oo}.
// THE DECIDER: the FULL genus-5 curve depends on b, and `ModelChecks` never tests it because it
// SKIPS CRV entries. Counting points on the fibre product over F_p and comparing with
// `ComputePointsViaTrace` (the Eichler-Selberg trace formula on the W-fixed part of the D-new
// space -- independent of Borcherds/Schofer entirely) settles it:
//     b = -1 matches at ALL 12 primes 7..47;
//     b = -2,-3,-5,-6,-7,-10,-15,-30,1,2,3,5,15 EACH FAIL at 4-6 of them.
// So b = -1 is uniquely determined by our own machinery. tests/CRV_15_4.m runs this check.
//
// VALIDATION: VerifyModelSet(15, 4) -> 16 checks, 0 failures; plus the full-curve point counts
// above, which are what pin b.
//
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,3]] := [* <2, -(4*x^2 - x + 1)*(4*x^2 + x + 1)*(5*x^2 + 3), P![]> *];
models[[Integers()|1,15]] := [* <0, -(3*x^2 + 1), P![]> *];
models[[Integers()|1]] := [* <5, "CRV", [ Strings() |
  "y^2 + (4*x^2 - x*s + s^2)*(4*x^2 + x*s + s^2)*(5*x^2 + 3*s^2)",
  "z^2 + 3*x^2 + s^2" ]> *];
