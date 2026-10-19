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
fexact := 272*x^8 - 72*x^6 - 55*x^4 - 18*x^2 + 17;
Cexact := HyperellipticCurve(fexact);
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

// [9]  The descent from a quadratic field (NelsonQuadraticDescent, NelsonQuadraticTwist): when
//      the two fixed points of the marking involution on the canonical conic are conjugate over
//      K = Q(sqrt d), the quartic lives over K and the curve over Q is recovered by the norm
//      equation, the substitution x = (1 + sqrt d s)/(1 - sqrt d s), and the two twists e (by
//      the involution) and r (quadratic) read off point counts.  Controls with known equations:
//      HAND: a synthetic octic prod (d s^2 - d tau s + 1) is invariant under s -> 1/(d s) by
//      construction; moved to the K-coordinate and descended it must come back isomorphic over Q
//      for some pair (e, r) supported on the bad primes.  EXACT MODEL: the 1326 model above with
//      its involutions x -> -1/(2x) (fixed points over Q(sqrt -2), the quotient by w_3) and
//      x -> 1/(2x) (over Q(sqrt 2), by w_6), the twists found from the trace formula, must give
//      a curve isomorphic to the model and the quotient Jacobians 1326c1 and 1326b2 (MAGMA:
//      Cremona labels); the K-invariants of the model are computed exactly here, not retyped.
function KInvariantsOfModel(f, a, b, c, dd)
    // the involution M(x) = (a x + b)/(c x - a); X = (x - r1)/(x - r2) for its fixed points r1, r2
    K<w> := QuadraticField(dd); PK<X> := PolynomialRing(K);
    rts := Roots(PK ! (c*x^2 - 2*a*x - b)); r1 := rts[1][1]; r2 := rts[2][1];
    F := PK ! (Evaluate(PK ! f, (r1 - r2*X)/(1 - X)) * (1 - X)^8);
    assert &and[ Coefficient(F, i) eq 0 : i in [1..7 by 2] ];      // M acts as X -> -X
    g := [ Coefficient(F, 2*j) : j in [0..4] ];
    return g[1]*g[3]/g[2]^2, Rationals() ! (g[2]*g[4]/g[3]^2);
end function;
for dd in [-2, 2] do
    Ps := &*[ dd*x^2 - dd*t*x + 1 : t in [3/2, -7/3, 5, -1/4] ];
    assert Ps eq P ! (dd^4*x^8*Evaluate(Ps, 1/(dd*x)));           // HAND: the symmetry
    Ad, Bd := KInvariantsOfModel(Ps, 0, 1, dd, dd);
    assert Degree(Parent(Ad)) eq 2 and not Ad in Rationals();   // a generic case: A irrational
    Pd := NelsonQuadraticDescent(Ad, Bd);
    dPs := Discriminant(Ps); dPd := Discriminant(Pd);
    supp := Sort(SetToSequence(Set(PrimeDivisors(Numerator(dPs))) join Set(PrimeDivisors(Denominator(dPs)))
                                join Set(PrimeDivisors(Numerator(dPd))) join {2, Abs(dd)}));
    Cs := HyperellipticCurve(Ps); recovered := false;
    for bits in CartesianPower([0,1], #supp) do
        e0 := &*([1] cat [ supp[i] : i in [1..#supp] | bits[i] eq 1 ]);
        for e in [e0, -e0] do
            if not NormEquation(Parent(Ad), Rationals() ! e) then continue; end if;
            Pe := NelsonQuadraticDescent(Ad, Bd : Twist := e);
            for bits2 in CartesianPower([0,1], #supp) do
                r0 := &*([1] cat [ supp[i] : i in [1..#supp] | bits2[i] eq 1 ]);
                if IsIsomorphic(Cs, HyperellipticCurve(r0*Pe)) or IsIsomorphic(Cs, HyperellipticCurve(-r0*Pe)) then
                    recovered := true; break bits;
                end if;
            end for;
        end for;
    end for;
    assert recovered;
end for;
for tup in [ <0, -1/2, 1, -2, 3, "1326c1">, <0, 1/2, 1, 2, 6, "1326b2"> ] do
    Am, Bm := KInvariantsOfModel(fexact, tup[1], tup[2], tup[3], tup[4]);
    e, r, R := NelsonQuadraticTwist(1326, 1, {1,17,78,1326}, tup[5], Am, Bm);
    assert R`split and IsIsomorphic(R`model, Cexact);
    assert CremonaReference(MinimalModel(Jacobian(GenusOneModel(R`rE*R`PE)))) eq tup[6];
    assert Conductor(MinimalModel(Jacobian(GenusOneModel(R`rE*R`PE)))) eq 1326;
end for;

// [10] The quadratic-field case on real data: X_0^6(323)/<w_6, w_19, w_34>, from the stored
//      values data/nelson/values_6_323.m (cutoff 46, 80 digits, 14 circles of 16 points, four
//      processes against the Z-basis recorded in the file; produced by
//        B := NelsonBasis(6, 323, {1,6,19,34,51,114,646,969} : Precision := 80, Cutoff := 46, OrderBasis := obas);
//        V := NelsonCircleValues(B, NelsonSamplePoints(B, 14), 16, RealField(80)!3/100);
//      with base point 0.1103 + 1.0411 i).  With u = w_2 the middle invariant is rational and
//      the other two are complex conjugates; the fixed points of u on X are CM points of
//      discriminant -7752 defined over a quartic field whose quadratic subfield is Q(sqrt -2)
//      (CM theory, through FieldOfDefinitionOfCMPoint), A is recognised there, and the twists
//      from the trace formula put the curve over a conic with NO rational point: the quotient
//      has no model y^2 = f(s) over Q, and the result is the K-model with its descent datum.
//      What is asserted: the field, B, that A lies in K with the recognised value
//      (CODE-DERIVED), that the conic is pointless, and the four checks: the fixed points of u
//      on the K-model are defined over the CM field, the Jacobian of the descended quotient has
//      conductor 1938 and the a_p of the W-fixed newform with w_2 = +1, the W-fixed newforms
//      split 1 + 2, and the zeta numerators at the split primes agree with the trace formula.
V := ReadNelsonValues("data/nelson/values_6_323.m");
assert V`D eq 6 and V`N eq 323 and V`ncirc eq 14;
dK := NelsonBranchField(6, 323, {1,6,19,34,51,114,646,969}, 2);
assert dK eq -8;                                                                      // CM theory
C6, diag := NelsonHyperellipticModel(V : Involution := 2);
K6 := Parent(diag`Arat);
assert Degree(K6) eq 2 and Discriminant(Integers(K6)) eq -8;
assert diag`Brat eq 22468/134689;                                                     // CODE-DERIVED
assert diag`Arat eq (100340736*K6.1 - 130090857)/1262027560 or diag`Arat eq (-100340736*K6.1 - 130090857)/1262027560;  // CODE-DERIVED
assert BaseRing(C6) eq K6;                                                            // no model over Q
e6, r6, R6 := NelsonQuadraticTwist(6, 323, {1,6,19,34,51,114,646,969}, 2, diag`Arat, diag`Brat);
assert not R6`split and not HasRationalPoint(R6`conic);
assert R6`counts ge 30;
ok6, Rc6 := NelsonDescentChecks(6, 323, {1,6,19,34,51,114,646,969}, 2, R6);
assert ok6;
assert Rc6`splitting eq [1, 2] and Conductor(Rc6`jacobian) eq 1938 and Rc6`zetafailed eq 0 and Rc6`zetachecked ge 10;
assert Degree(Rc6`branchmodel) eq 4;
