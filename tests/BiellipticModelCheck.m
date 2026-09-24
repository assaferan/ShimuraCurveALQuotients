// tests/BiellipticModelCheck.m
//
// CI checks for BiellipticModelCheck.m: deciding between candidate genus-2 bielliptic models
// y^2 = A x^6 + B x^4 + C x^2 + D of an Atkin-Lehner quotient X_0(D,N)/W by the fields of
// definition of the fixed points of the residual Atkin-Lehner involutions (Shimura reciprocity).
//
//   [1] ALGroupFromGenerators closes a generating set under Atkin-Lehner multiplication,
//   [2] MatchFixedPointOrbits: Galois orbits of observed fixed points must be absorbed by the
//       predicted (CM discriminant, count, possible fields) rows -- multiset, not "some row",
//   [3] ReadBiellipticCandidates parses the vendored candidate file,
//   [4] X_0(34,3)/w_102: the second candidate is consistent, the first contradicts (the case
//       worked out by hand in verify_al_fixed_fields.m; both have the trace-formula a_p),
//   [5] CheckBiellipticEntry returns a "determined" verdict there and "not attempted" for a
//       non-squarefree level.

QQ := Rationals();
P<x> := PolynomialRing(QQ);

// [1]
assert ALGroupFromGenerators({3}, 102) eq {1, 3};
assert ALGroupFromGenerators({102}, 102) eq {1, 102};
assert ALGroupFromGenerators({2, 3}, 210) eq {1, 2, 3, 6};
assert ALGroupFromGenerators({2, 3, 133}, 798) eq {1, 2, 3, 6, 133, 266, 399, 798};
assert ALGroupFromGenerators({Integers() | }, 102) eq {1};

// [2]  expected rows are <disc, count on the quotient, [* possible fields *]>
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

// [3]
entries := ReadBiellipticCandidates("data/bielliptic_candidates.m");
assert #entries eq 231;
e := entries[1];
assert e[1] eq 6 and e[2] eq 17 and e[3] eq {3} and #e[4] eq 2 and Degree(e[4][1]) eq 6;
e102 := [* e : e in entries | e[1] eq 34 and e[2] eq 3 and e[3] eq {102} *][1];
assert e102[4][2] eq 36*x^6 + 117*x^4 + 18*x^2 - 27;

// [4]
X := ALQuotientFromGenerators(34, 3, {102});
assert X`g eq 2;
expected := ExpectedALFixedPointData(X);
assert #expected eq 3;
assert Sort([r`QuotientGenus : r in expected]) eq [0, 1, 1];
ok1, rep1 := CheckBiellipticCandidate(X, e102[4][1], expected);
ok2, rep2 := CheckBiellipticCandidate(X, e102[4][2], expected);
assert not ok1 and ok2;
assert "CONTRADICTION" in rep1 and "CONTRADICTION" notin rep2;

// [5]
v := CheckBiellipticEntry(e102);
assert v`Status eq "determined" and v`Consistent eq [2];
v9 := CheckBiellipticEntry([* 34, 9, {2, 153}, [x^6 + 1] *]);
assert v9`Status eq "not attempted";
