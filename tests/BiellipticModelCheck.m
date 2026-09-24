// tests/BiellipticModelCheck.m
//
// CI checks for BiellipticModelCheck.m: deciding between candidate genus-2 bielliptic models
// y^2 = A x^6 + B x^4 + C x^2 + D of an Atkin-Lehner quotient X_0(D,N)/W by the fields of
// definition of the fixed points of the residual Atkin-Lehner involutions (Shimura reciprocity).
//
// Every expected value is labelled with where it comes from:
//   HAND        worked out by hand from the definitions; the argument is in the comment
//   FILE        read off the vendored data/bielliptic_candidates.m (a fixed copy of
//               fsaia/GenusAtMost2, fetched 2026-09-17), checked against the raw text
//   MAGMA       a Magma built-in (not this repo's code)
//   CROSS-CHECK two different methods are compared; neither is an independent truth
//   CODE-DERIVED a regression value produced by this code, with no independent source; it pins
//               behaviour, it does not validate it
//
//   [1] ALGroupFromGenerators closes a generating set under Atkin-Lehner multiplication,
//   [2] MatchFixedPointOrbits: Galois orbits of observed fixed points must be absorbed by the
//       predicted (CM discriminant, count, possible fields) rows -- multiset, not "some row",
//   [3] ReadBiellipticCandidates parses the vendored candidate file,
//   [4] X_0(34,3)/w_102: candidate 1 contradicts (independent ramification argument),
//       candidate 2 is consistent (code-derived),
//   [5] CheckBiellipticEntry returns "determined" there and "not attempted" for a non-squarefree
//       level,
//   [6] ModelInvolutionCheck: which bielliptic involution of the model is Atkin-Lehner, when the
//       reduced automorphism group is larger than C_2; and the downgrade it triggers.

QQ := Rationals();
P<x> := PolynomialRing(QQ);

// [1]  HAND: w_a w_b = w_{ab/gcd(a,b)^2}.  {2,3} in 210 closes to {1,2,3,6};
//      {2,3,133} in 798 = 2*3*7*19 closes to all products of 2, 3, 133.
assert ALGroupFromGenerators({3}, 102) eq {1, 3};
assert ALGroupFromGenerators({102}, 102) eq {1, 102};
assert ALGroupFromGenerators({2, 3}, 210) eq {1, 2, 3, 6};
assert ALGroupFromGenerators({2, 3, 133}, 798) eq {1, 2, 3, 6, 133, 266, 399, 798};
assert ALGroupFromGenerators({Integers() | }, 102) eq {1};

// [2]  HAND: synthetic inputs; the expected answer follows from the definition (each orbit of
//      degree k takes k points from one row whose field list contains its field; every row must
//      be filled exactly).  Rows are <disc, count on the quotient, [* possible fields *]>.
Qm3 := QuadraticField(-3);
K4 := NumberField(x^4 + 1/4*x^2 - 1/4);
exp_ell := [* <-51, 2, [* Qm3 *]> *];
assert MatchFixedPointOrbits([* Qm3 *], exp_ell);                     // one quadratic orbit
assert not MatchFixedPointOrbits([* QQ, QQ *], exp_ell);              // two rational points: no
assert not MatchFixedPointOrbits([* QuadraticField(-222) *], exp_ell);
exp_rat := [* <-3, 2, [* QQ *]> *];
assert MatchFixedPointOrbits([* QQ, QQ *], exp_rat);
assert not MatchFixedPointOrbits([* Qm3 *], exp_rat);                 // Q(sqrt-3) is not Q
exp_hyp := [* <-24, 2, [* Qm3 *]>, <-68, 4, [* K4 *]> *];
assert MatchFixedPointOrbits([* Qm3, K4 *], exp_hyp);
assert not MatchFixedPointOrbits([* Qm3, Qm3, QuadraticField(17) *], exp_hyp);  // -24 has only 2
assert not MatchFixedPointOrbits([* QuadraticField(3), K4 *], exp_hyp);
// a rational orbit can also sit in a row whose field list contains a degree-1 field
assert MatchFixedPointOrbits([* QQ, QQ *], [* <-3, 2, [* NumberField(x - 1 : DoLinearExtension) *]> *]);

// [3]  FILE: the vendored file has 231 entries.  This is a count of that fixed file, not a
//      mathematical fact; it is cross-checked here against the raw text (one "[* D, N," opener
//      per entry), which does not go through ReadBiellipticCandidates.
entries := ReadBiellipticCandidates("data/bielliptic_candidates.m");
raw := Read("data/bielliptic_candidates.m");
raw := raw[Position(raw, "\ngenus_2_bielliptics_eqn_not_determined")..#raw];
assert #entries eq 231;
assert #[i : i in [1..#raw-1] | raw[i] eq "[" and raw[i+1] eq "*"] eq 231;   // one "[*" per entry
// FILE: the first entry is [* 6, 17, { 3 }, [ two sextics ] *]; entry (34,3,{102}) has
//       36*x^6 + 117*x^4 + 18*x^2 - 27 as its second candidate (read off the file text).
e := entries[1];
assert e[1] eq 6 and e[2] eq 17 and e[3] eq {3} and #e[4] eq 2 and Degree(e[4][1]) eq 6;
e102 := [* e : e in entries | e[1] eq 34 and e[2] eq 3 and e[3] eq {102} *][1];
assert e102[4][1] eq -1152/37*x^6 + 2709/37*x^4 + 2529/37*x^2 - 864/37;
assert e102[4][2] eq 36*x^6 + 117*x^4 + 18*x^2 - 27;

// [4]  X_0(34,3)/w_102, residual AL group W_full/W of order 4 (HAND: W_full has 2^3 elements,
//      W = {1, 102}), so three nontrivial cosets.
X := ALQuotientFromGenerators(34, 3, {102});
assert X`g eq 2;                                   // FILE: the entry is a genus-2 quotient
expected := ExpectedALFixedPointData(X);
assert #expected eq 3;                             // HAND: index of W in W_full is 4
// CROSS-CHECK: genera from GenusShimuraCurveQuotient; a V_4 containing iota acting on a genus-2
// curve has quotients of genus 0, 1, 1 (Riemann-Hurwitz), so the two must agree.
assert Sort([r`QuotientGenus : r in expected]) eq [0, 1, 1];
//
// HAND: candidate 1 is NOT the model.  Its geometric automorphism group is V_4 (checked below via
// MAGMA), so the AL involutions are sigma : x -> -x and sigma*iota, whose fixed points are
// (0, +-sqrt(D)) and the two points at infinity with y/x^3 = +-sqrt(A).  Here A = -1152/37 is
// -74 times a square and D = -864/37 is -222 times a square, so these fixed points would be
// defined over Q(sqrt(-74)) and Q(sqrt(-222)), both ramified at 37.  But a fixed point of a
// residual AL involution on C is the image of a fixed point of some w_m, m | 102, on
// X_0(34,3); those are CM points by orders containing Z[sqrt(-m)] or Z[(1+sqrt(-m))/2] (and
// Z[i] for m = 2), of discriminant dividing 4m | 408 = 2^3*3*17; they are defined over ring
// class fields, which are unramified outside the discriminant, and so is the field of their
// image on C.  37 does not divide 408.  (There are no cusps: 34 is a quaternion discriminant.)
f1 := e102[4][1];
assert IsSquare(Coefficient(f1, 6) / -74) and IsSquare(Coefficient(f1, 0) / -222);        // MAGMA
assert &and[Discriminant(MaximalOrder(QuadraticField(d))) mod 37 eq 0 : d in [-74, -222]]; // MAGMA
assert #GeometricAutomorphismGroup(HyperellipticCurve(f1)) eq 4;                          // MAGMA
// CROSS-CHECK the premise against the code's own expected data: every predicted CM
// discriminant divides 408
assert &and[&and[408 mod AbsoluteValue(row[1]) eq 0 : row in r`Rows] : r in expected];
ok1, rep1 := CheckBiellipticCandidate(X, f1, expected);
assert not ok1 and "CONTRADICTION" in rep1;
assert "Q(sqrt(-222))" in rep1 and "Q(sqrt(-74))" in rep1;   // it fails for the reason above
// CODE-DERIVED: candidate 2 is consistent.  No independent source -- "consistent" only means the
// code found no contradiction; this pins the behaviour.
ok2, rep2 := CheckBiellipticCandidate(X, e102[4][2], expected);
assert ok2 and "CONTRADICTION" notin rep2;

// [5]  "determined" = candidate 1 excluded (HAND, [4]) + candidate 2 not excluded (CODE-DERIVED).
v := CheckBiellipticEntry(e102);
assert v`Status eq "determined" and v`Consistent eq [2] and v`Untrusted eq [];
// HAND: N = 9 is not squarefree, outside the stated scope
v9 := CheckBiellipticEntry([* 34, 9, {2, 153}, [x^6 + 1] *]);
assert v9`Status eq "not attempted";

// [6]  ModelInvolutionCheck.  The AL involutions are defined over Q, so they sit in Aut_Q(C);
//      the model's sigma is guaranteed to be one of them only if every Klein four-subgroup of
//      Aut_Q containing iota is conjugate to <iota, sigma>.
//
// (a) the generic case: reduced group C_2.  CODE-DERIVED (the "V_4" verdict is this repo's
//     ModelInvolutionCheck); the geometric group of order 4 it should agree with is MAGMA, asserted in [4].
ok, desc := ModelInvolutionCheck(e102[4][2]);
assert ok and "V_4" in desc;
//
// (b) X_0(14,15)/<w_7,w_30>, candidate 2 -- the ONLY attempted (squarefree-N) candidate in the
//     file with extra automorphisms (CODE-DERIVED survey, see the BiellipticModelCheck.m header).
//     HAND: f = 9x^6 + 90x^4 + 225x^2 + 252 satisfies (x+1)^6 f((x-3)/(x+1)) = 64 f(x), so
//     x -> (x-3)/(x+1), y -> 8y/(x+1)^3 is an automorphism over Q; its matrix [[1,-3],[1,1]] has
//     cube -8*I, so it has order 3 in PGL_2 and the reduced group is at least S_3, not C_2.
//     With #Aut_Q = 12 (MAGMA), Aut_Q = C_2 x S_3, and the three involutions of S_3 are
//     conjugate, so every Klein subgroup containing iota is conjugate to <iota, sigma>: TRUSTED.
e1415 := [* e : e in entries | e[1] eq 14 and e[2] eq 15 and e[3] eq {7, 30} *][1];
f12 := e1415[4][2];
assert f12 eq 9*x^6 + 90*x^4 + 225*x^2 + 252;                                             // FILE
assert P!((x + 1)^6 * Evaluate(f12, (x - 3)/(x + 1))) eq 64*f12;                          // HAND
assert Matrix(QQ, 2, 2, [1, -3, 1, 1])^3 eq ScalarMatrix(QQ, 2, -8);                      // HAND
assert #AutomorphismGroup(HyperellipticCurve(f12)) eq 12;                                 // MAGMA
ok, desc := ModelInvolutionCheck(f12);
assert ok and "<12, 4>" in desc and "1 conjugacy class" in desc;
// CODE-DERIVED: the entry stays determined, with no untrusted candidate
v := CheckBiellipticEntry(e1415);
assert v`Status eq "determined" and v`Consistent eq [2] and v`Untrusted eq [];
//
// (c) an UNTRUSTED model: X_0(15,4)/<w_3>, candidate 1 (non-squarefree N, so never attempted by
//     the pipeline; used here only as an even sextic with the right automorphisms).
//     HAND: f = (576x^6 + 639x^4 + 639x^2 + 576)/11 is palindromic, so tau : (x,y) -> (1/x, y/x^3)
//     is an automorphism.  sigma tau (x,y) = (-1/x, y/x^3) and tau sigma (x,y) = (-1/x, -y/x^3),
//     so tau sigma tau^-1 = iota sigma: the conjugacy class of sigma is {sigma, iota sigma} and
//     <iota, sigma> is normal in <sigma, tau> = D_8, while <iota, tau> is a different Klein
//     subgroup.  The geometric group has order 8 (MAGMA), so Aut_Q = D_8 exactly and the two
//     Klein subgroups containing iota are not conjugate: the AL V_4 could be <iota, tau>.
f8 := entries[148][4][1];
assert entries[148][1] eq 15 and entries[148][2] eq 4;                                    // FILE
assert f8 eq (576*x^6 + 639*x^4 + 639*x^2 + 576)/11;                                      // FILE
assert P!(x^6 * Evaluate(f8, 1/x)) eq f8;                                                 // HAND
assert #GeometricAutomorphismGroup(HyperellipticCurve(f8)) eq 8;                          // MAGMA
ok, desc := ModelInvolutionCheck(f8);
assert not ok and "2 conjugacy class" in desc;
//
// (d) NEGATIVE CONTROL for the downgrade: X_0(14,33)/<w_2,w_3,w_77> has a residual AL group of
//     order 2 whose coset has a genus-1 quotient (no hyperelliptic AL coset, so no Weierstrass
//     check), and candidate 1 is excluded ONLY by the fixed-point fields of sigma.  Forcing the
//     model to count as untrusted must therefore keep candidate 1 alive (HAND, from the rule),
//     turning "determined [2]" (CODE-DERIVED) into "ambiguous [1, 2]".
e1433 := [* e : e in entries | e[1] eq 14 and e[2] eq 33 and e[3] eq {2, 3, 77} *][1];
v := CheckBiellipticEntry(e1433);
assert v`Status eq "determined" and v`Consistent eq [2] and v`Untrusted eq [];
vu := CheckBiellipticEntry(e1433 : AssumeUntrusted := true);
assert vu`Status eq "ambiguous" and vu`Consistent eq [1, 2] and vu`Untrusted eq [1, 2];
assert "NOT EXCLUDED (untrusted model)" in vu`Report;
// ... whereas an involution-independent contradiction still excludes an untrusted model: in
// (34,3)/w_102 candidate 1 also fails the Weierstrass-point check (iota is unique).
vu := CheckBiellipticEntry(e102 : AssumeUntrusted := true);
assert vu`Status eq "determined" and vu`Consistent eq [2];
