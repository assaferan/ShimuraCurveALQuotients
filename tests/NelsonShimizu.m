// tests/NelsonShimizu.m
//
// Checks for NelsonShimizu.m: the W-fixed differentials of an Atkin-Lehner quotient evaluated by
// the Shimizu lift (Nelson, Math. Comp. 84 (2015), Theorem 3.1) with the oldform extension,
// Noether's count of canonical quadrics from the values, and the genus-3 model through a marked
// Atkin-Lehner involution.  Everything the package computes is numerical; what is asserted here
// is that the numerical answer agrees with a stored exact fact.
//
// Every expected value is labelled with where it comes from:
//   FILE         a stored model in data/models/ (read there, not retyped)
//   HAND         follows from the definitions; the argument is in the comment
//   MAGMA        a Magma built-in (not this repo's code)
//   EXACT MODEL  the equation of X_0^1326(1)/<w_17, w_78> obtained by an exact construction
//                (the double cover of the elliptic curve 1326b2 branched at the CM points of
//                discriminants -24 and -408 prescribed by the heights of Heegner points), which
//                agrees with the trace formula at all primes up to 1000
//   CODE-DERIVED a value produced by this code with no independent source; a regression value
//
// Cutoffs are kept small so that the whole file runs in a few minutes; the verdicts need only a
// clean gap, not many digits.  Each verdict fixes its Eichler order through a serialised Z-basis
// (CODE-DERIVED, one order among the conjugates Magma returns), because the count at a small
// cutoff depends on which order, hence which points of the curve, the sample points name; at a
// larger cutoff any order serves, at the cost of a longer test.

QQ := Rationals();
P<x> := PolynomialRing(QQ);

// [1]  HAND: an Eichler order serialised and rebuilt is an Eichler order of the right level,
//      and the rebuilt one contains the basis it was built from.  The point of the round trip is
//      that QuaternionOrder(A, N) is not the same order from one process to the next.
obas := NelsonEichlerOrderBasis(46, 3);
O := NelsonEichlerOrder(46, 3, obas);
assert Discriminant(O) eq 138;
assert #obas eq 4 and forall{ v : v in obas | #v eq 4 };

// [2]  FILE: X_0^35(1) is hyperelliptic -- data/models/models_35_1.m stores y^2 = f(x) of genus 3
//      at the key [1].  Noether's count at genus 3 must find ONE quadric through the canonical
//      image.  The count is taken at a small cutoff with ten points; the gap it reports is
//      CODE-DERIVED and only its order of magnitude is asserted.
m35 := eval (Read("data/models/models_35_1.m") cat "\nreturn models;");
assert m35[[Integers()|1]][1][1] eq 3 and Degree(m35[[Integers()|1]][1][2]) eq 8;
ob35 := [ [ 1, 0, 0, 0 ], [ 1/2, 0, 1/2, 1/2 ], [ 1/2, 0, 1/2, -1/2 ], [ 0, 1/4, 1/2, 1/4 ] ];
verdict, R := NelsonHyperellipticVerdict(35, 1, {1} : Precision := 30, Cutoff := 28, Points := 10, OrderBasis := ob35);
assert verdict eq "HYPERELLIPTIC";
assert R`relations eq 1 and R`gap gt 10^5;
assert R`basis`oldcount eq 0;                  // HAND: at N = 1 every differential is a newform

// [3]  FILE: X_0^46(3)/<w_3> is NOT hyperelliptic -- data/models/models_46_3.m stores it at the
//      key [1,3] as a fibre product of two genus-1 covers (a CRV entry), and the exact Noether
//      count on the function field of that stored curve gives no quadric: its canonical model is
//      a smooth plane quartic.  One of its three differentials is an oldform (level 46), so this is also the
//      control for the oldform extension: without it the package must refuse (HAND: two newforms
//      of level 138 against genus 3), with it the count must find NO quadric, and the measured
//      Atkin-Lehner eigenvalue of every column under w_3 must be the predicted one.  Twelve
//      points at cutoff 30 against the order below: with fewer points, a smaller cutoff, or some
//      other orders, the spectrum is a ramp and the count is inconclusive rather than wrong.
m46 := eval (Read("data/models/models_46_3.m") cat "\nreturn models;");
assert m46[[Integers()|1,3]][1][1] eq 3 and m46[[Integers()|1,3]][1][2] eq "CRV";
refused := false;
try
    _ := NelsonBasis(46, 3, {1,3} : Precision := 30, Cutoff := 30, Oldforms := false);
catch e
    refused := true;
end try;
assert refused;
ob46 := [ [ 1, 0, 0, 0 ], [ 1/2, 1, 1/6, 4/3 ], [ 0, 1/2, 0, 3/2 ], [ 0, 3/2, 0, 3/2 ] ];
verdict, R := NelsonHyperellipticVerdict(46, 3, {1,3} : Precision := 30, Cutoff := 30, Points := 12, CheckW := 1, OrderBasis := ob46);
assert verdict eq "NOT hyperelliptic";
assert R`relations eq 0 and R`singular[1] gt 10^(-3);
assert R`basis`oldcount eq 1 and R`basis`newcount eq 2;
assert R`wcheck lt 10^(-6);

// [4]  FILE: X_0^26(3)/<w_39> is hyperelliptic -- data/models/models_26_3.m stores y^2 = f(x) of
//      genus 3 at the key [1,39].  Two of its three differentials are oldforms (level 26), with
//      characters of opposite sign under w_3, so the relative scale and sign of the pullback
//      g|w_3 enter the count.  One quadric, and every column's w_39-eigenvalue measured as +1.
m26 := eval (Read("data/models/models_26_3.m") cat "\nreturn models;");
assert m26[[Integers()|1,39]][1][1] eq 3 and Degree(m26[[Integers()|1,39]][1][2]) eq 8;
ob26 := [ [ 1, 0, 0, 0 ], [ 1/2, 1, 1/10, 11/10 ], [ 0, 0, 0, 1 ], [ 0, 3/2, 0, 3/2 ] ];
verdict, R := NelsonHyperellipticVerdict(26, 3, {1,39} : Precision := 30, Cutoff := 28, Points := 10, CheckW := 1, OrderBasis := ob26);
assert verdict eq "HYPERELLIPTIC";
assert R`relations eq 1 and R`gap gt 10^5;
assert R`basis`oldcount eq 2 and R`basis`newcount eq 1;
assert R`wcheck lt 10^(-6);

// [5]  EXACT MODEL: the model step on X_0^1326(1)/<w_17, w_78>, from the stored values
//      data/nelson/values_1326_1.m (cutoff 46, 80 digits, 14 circles of 16 points, far too
//      expensive to recompute in a test; produced by
//        B := NelsonBasis(1326, 1, {1,17,78,1326} : Precision := 80, Cutoff := 46, OrderBasis := obas);
//        V := NelsonCircleValues(B, NelsonSamplePoints(B, 14), 16, RealField(80)!3/100);
//      with obas the Z-basis recorded in the file, and base point 0.1103 + 1.0411 i).
//      The involution u = w_2 has a genus-1 quotient; the invariants recognised, the twist
//      constants solved, and the resulting curve must be isomorphic over Q to the exact model
//      y^2 = 272x^8 - 72x^6 - 55x^4 - 18x^2 + 17.  The invariants, kappa and mu themselves are
//      CODE-DERIVED regression values.
V := ReadNelsonValues("data/nelson/values_1326_1.m");
assert V`D eq 1326 and V`N eq 1 and V`W eq [1,17,78,1326] and V`ncirc eq 14 and #V`pts eq 238;
C, diag := NelsonHyperellipticModel(V : Involution := 2);
assert diag`u eq 2;
assert diag`Arat eq -935/324 and diag`Brat eq 1296/3025 and diag`Crat eq -935/324;   // CODE-DERIVED
assert diag`kappa eq 17 and diag`mu eq -34;                                          // CODE-DERIVED
assert Genus(C) eq 3;
Cexact := HyperellipticCurve(272*x^8 - 72*x^6 - 55*x^4 - 18*x^2 + 17);
assert IsIsomorphic(C, Cexact);                                                      // EXACT MODEL
// the diagnostics must show what the construction relies on: a clean conic, u-equivariance of
// the conic, an even octic (the quartic in x^2 predicts the centres it was not fitted on)
assert diag`conicsv[2]/diag`conicsv[1] gt 10^10;
assert diag`equivariance lt 10^(-10);
assert diag`worstprediction lt 10^(-4);

// [6]  The four checks beyond point counts on that model.  (a) HAND + CM theory: the coset w_2 W
//      = {2, 34, 39, 663}; its elements have 16 fixed points on X_0^1326(1) in all, 4 on the
//      quotient as Riemann-Hurwitz demands for a double cover of a genus-1 curve by a genus-3
//      curve; the model puts both pairs over Q(sqrt 17) (kappa = 17 and a leading coefficient
//      17 times a square), and the fields of definition of the CM points must agree -- the
//      discriminant -136 carries all four points, defined over Q(sqrt 17), so it contributes
//      that field twice.  (b) MAGMA: the Jacobian of the quartic has conductor 1326 and the
//      a_p of the W-fixed newform with quaternionic w_2 = +1.  (c) the W-fixed newforms split
//      1 + 2 under w_2.  (d) the zeta numerators over F_p, p <= 60, agree with the trace formula.
ok, R := NelsonModelChecks(1326, 1, {1,17,78,1326}, 2, [diag`Arat, diag`Brat, diag`Crat], [diag`kappa, diag`mu]);
assert ok;
assert R`branchmodel eq {* 17, 17 *} and R`branchcm eq {* 17, 17 *};
assert R`splitting eq [1, 2];
assert Conductor(R`jacobian) eq 1326;
assert R`zetafailed eq 0 and R`zetachecked ge 20;

// [7]  HAND: the other two involutions of that quotient (w_3, w_6) also have genus-1 quotients,
//      but at these fourteen centres the coordinate they fix is badly spread (the moduli |x| at
//      the centres span less than a factor of five) and the invariants are not recognised; the
//      package must refuse rather than return a curve, and NelsonCoordinateSpread must show the
//      spread the refusal comes from.  The spread values are CODE-DERIVED.
spread := NelsonCoordinateSpread(V : Involution := 3);
assert #spread eq 14 and Max(spread)/Min(spread) lt 8;
refused := false;
try
    _ := NelsonHyperellipticModel(V : Involution := 3);
catch e
    refused := true;
end try;
assert refused;

// [8]  Round trip of the value file: written and read back, the record is the same.
tmp := Tempname("nelson_values_");
WriteNelsonValues(V, tmp);
V2 := ReadNelsonValues(tmp);
assert V2`pts eq V`pts and V2`vals eq V`vals and V2`obas eq V`obas and V2`alsig eq V`alsig;
System("rm -f " cat tmp);
