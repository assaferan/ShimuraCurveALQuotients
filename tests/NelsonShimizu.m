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
//   EXACT MODEL  an equation read from data/models/, obtained by an exact construction: the
//                double cover of an elliptic quotient branched at the CM points prescribed by the
//                heights of Heegner points (for X_0^1326(1)/<w_17, w_78>, over 1326b2 with branch
//                points of discriminants -24 and -408, agreeing with the trace formula at all
//                primes up to 1000).  How each stored model was built is recorded in
//                data/models/PROVENANCE.md
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
//      The non-hyperelliptic reading is the one a precision floor can fake, so the code refuses it
//      without a second precision: the count is repeated at 45 digits at the same points, and the
//      smallest singular value must stay put (it is 8.4e-3 here, nowhere near a floor).  At
//      X_0^4130(1) that same check overturned a reading: 1.9e-3 at three cutoffs at 50 digits,
//      2.0e-10 at 80 (campaign record, second batch).
verdict, R := NelsonHyperellipticVerdict(46, 3, {1,3} : Precision := 30, Cutoff := 30, Points := 12,
                                         CheckW := 1, OrderBasis := ob46, ConfirmPrecision := 45);
assert verdict eq "NOT hyperelliptic";
assert R`confirmed and R`verdictatconfirm eq "NOT hyperelliptic";

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
m1326 := eval (Read("data/models/models_1326_1.m") cat "\nreturn models;");
key1326 := [Integers()|1,17,78,1326];
e1326 := m1326[key1326][1];
assert e1326[1] eq 3 and e1326[3] eq 0;
fexact := e1326[2];
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
assert V`W eq [1,6,19,34,51,114,646,969] and V`g eq 3 and V`cut eq 46 and V`prec eq 80;
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

// [11] HAND: the counting rule at genus 4, where the suite has no curve.  The canonical image of
//      a hyperelliptic curve of genus 4 is a rational normal curve in P^3, which lies on exactly
//      three quadrics -- the 2x2 minors of the 2x3 matrix of its coordinates -- so values
//      (1, t, t^2, t^3) at twelve points must give the count 3.  The three vanishing singular
//      values there differ among themselves by large factors (two are exactly zero at this
//      precision and the first nonzero one is 0.019), which is the configuration a rule reading
//      off the largest ratio anywhere would split into one relation.
//      The second configuration is a sample of the quadric surface u1 u4 = u2 u3 that nearly
//      satisfies a second quadric: its rank is 2, which no canonical image of genus 4 can have,
//      and the rule must refuse the count rather than report the nearest admissible one.
CC := ComplexField(30);
rnc := [ [ CC | sc*1, sc*t, sc*t^2, sc*t^3 ] where sc := CC ! (1 + t/7)
         : t in [ CC | 1/2, -2/3, 3, -5/4, 7/5, -8/3, 11/4, -13/5, 17/6, -19/7, 23/8, 29/9 ] ];
cnt4, sv4, gap4 := NoetherQuadricCount(rnc, 4);
assert cnt4 eq 3 and gap4 gt 10^4;
assert #[ i : i in [1..#sv4] | sv4[i] lt 10^(-4) ] eq 3;
quad := [ [ CC | a, b, c, b*c/a ] where a := CC ! (1 + i/3) where b := CC ! (2 - j/5) where c := CC ! (1 + i*j/11)
          : i in [1..4], j in [1..3] ];
cnt2, sv2, gap2 := NoetherQuadricCount(quad, 4);
assert #[ i : i in [1..#sv2] | sv2[i] lt 10^(-4) ] eq 2;
assert cnt2 eq -1;

// [12] The precision rule: a request whose evaluations would sum more than about a million terms
//      at fewer than 50 digits is refused, which is the X_0^1326(1) setting whose residual
//      stalled (about 3.3 million terms per evaluation at cutoff 46).  HAND: the term count is
//      10 (cut/4pi)^4 sigma_2(DN)/DN, the law measured in the write-up.
refused := false;
try
    _ := NelsonBasis(1326, 1, {1,17,78,1326} : Precision := 30, Cutoff := 46);
catch e
    refused := true;
end try;
assert refused;

// [13] EXACT MODEL: the fit checked against a model known exactly at a base where the certainty
//      criterion REFUSES the invariants, which is what tells a height limit from an error.  The
//      stored model of X_0^6(185)/<w_2, w_15, w_111> is an even octic, so the residual
//      Atkin-Lehner involution is already x -> -x there and its three marked invariants are
//      exact rationals.  The banked values (data/nelson/values_6_185.m, cutoff 46, 80 digits, 22
//      centres; produced by
//        B := NelsonBasis(6, 185, {1,2,15,30,111,185,222,370} : Precision := 80, Cutoff := 46, OrderBasis := obas);
//        V := NelsonCircleValues(B, ctr, 16, RealField(80)!3/100);
//      with obas the Z-basis in the file and ctr fourteen points from NelsonSamplePoints together
//      with eight where |x| is large) must reproduce them, and they do to about 2e-12 relative.
//      Which of the two outer invariants is which is not determined: exchanging the two roots of
//      the binary form is x -> 1/x, which reverses the quartic, so the fit fixes the pair and the
//      middle invariant, and either ordering is correct.  The relative differences are
//      CODE-DERIVED; the exact invariants are not.
//      The refusal is then shown to be correct rather than conservative: the criterion accepts a
//      fraction p/q only when q^2 eps < 1/2 with eps the noise it reads off the imaginary part,
//      and here that number is about 15 for the invariant of denominator 15792200 and 0.8 for the
//      one of denominator 3721, both refused, against 0.003 for the one of denominator 137288,
//      which it accepted.  The margin at 3721 is thin, which is the measure of how much more
//      accuracy that base would need.
V185 := ReadNelsonValues("data/nelson/values_6_185.m");
assert V185`D eq 6 and V185`N eq 185 and V185`g eq 3 and V185`cut eq 46 and V185`prec eq 80;
assert V185`W eq [1,2,15,30,111,185,222,370] and V185`ncirc eq 22;
m185 := eval (Read("data/models/models_6_185.m") cat "\nreturn models;");
e185 := m185[[Integers()|1,2,15,30,111,185,222,370]][1];
assert e185[1] eq 3 and e185[3] eq 0;
f185 := e185[2];
assert &and[ Coefficient(f185, i) eq 0 : i in [1,3,5,7] ];      // even: u acts as x -> -x
cf := [ Coefficient(f185, 2*j) : j in [0..4] ];
exA := cf[1]*cf[3]/cf[2]^2; exB := cf[2]*cf[4]/cf[3]^2; exC := cf[3]*cf[5]/cf[4]^2;
assert Denominator(exA) eq 15792200 and Denominator(exB) eq 3721 and Denominator(exC) eq 137288;
A185, B185, C185, um185, d185 := NelsonMarkedInvariants(V185);
assert um185 eq 3;
CCp := Parent(A185);
rel := func< z, q | Abs(z - CCp ! q)/Abs(CCp ! q) >;
assert rel(B185, exB) lt 10^(-10);
direct := rel(A185, exA) lt 10^(-10) and rel(C185, exC) lt 10^(-10);
reversed := rel(A185, exC) lt 10^(-10) and rel(C185, exA) lt 10^(-10);
assert direct or reversed;
// the certainty criterion's decision, reproduced from the noise and the exact denominators
qeps := func< z, q | Denominator(q)^2 * 20 * Abs(Imaginary(z)) >;
tall := direct select A185 else C185;        // the invariant of denominator 15792200
short := direct select C185 else A185;       // the invariant of denominator 137288
assert qeps(tall, exA) gt 1;                 // refused, and not marginally
assert qeps(B185, exB) gt 1/2;               // refused
assert qeps(short, exC) lt 1/100;            // accepted

// [14] EXACT MODEL against CODE-DERIVED: at the three bases where this pipeline and the exact
//      construction have both produced a model, they must agree over Q.  The equations below are
//      what this code printed (regression values, with the commands in data/nelson/ and in the
//      write-up); the curves they are compared with are read from data/models/, where they were
//      put by the Heegner-height construction, which shares no step with the evaluation here.
//      X_0^1326(1) is also the model recovered in [5] from the banked values.
pairs := [* <[Integers()|1,11,39,429], "data/models/models_39_11.m",
             HyperellipticCurve(-3*x^7-22*x^6-55*x^5-69*x^4-55*x^3-22*x^2-3*x, x^4+x^3+x^2+x+1)>,
           <[Integers()|1,35,55,77], "data/models/models_55_7.m",
             HyperellipticCurve(x^8+5*x^7+13*x^6+21*x^5+24*x^4+19*x^3+x^2-6*x-2, x^4+x+1)>,
           <key1326, "data/models/models_1326_1.m", C> *];
for pr in pairs do
    ms := eval (Read(pr[2]) cat "\nreturn models;");
    es := ms[pr[1]][1];
    assert es[1] eq 3;
    assert IsIsomorphic(pr[3], HyperellipticCurve(es[2] + es[3]^2/4));
end for;

// [15] EXACT MODEL: the boundary of what the numerical route can reach, at
//      X_0^6(385)/<w_5, w_6, w_11, w_14>.  The fit is as clean as anywhere -- conic residual
//      1.6e-15, predictions to 3e-9 -- and the three computed invariants agree with the stored
//      exact model to relative 1.4e-13, 1.3e-12 and 2.1e-11.  But their heights put two of them
//      beyond any recognition: the criterion accepts p/q only when q^2 eps < 1/2, and with
//      denominators 1 805 524 232 and 22 477 081 that number is about 35 000 and 870 here,
//      against 0.37 for the third, of denominator 168 200.  Since eps falls roughly like the
//      accuracy of the values, closing a factor of 35 000 would need four or five more digits
//      everywhere, which no feasible cutoff at D*N = 2310 provides.  So this base is not a case
//      of needing a little more precision, as X_0^6(185) is with its margin of 0.8 -- it is where
//      the route stops, and the model has to come from somewhere else.
//      The banked values (cutoff 46, 90 digits, 14 centres) were evaluated in four processes
//      before the Galois conjugates of the quadratic newform orbit were ordered canonically, and
//      the two conjugate columns of one chunk were exchanged afterwards; the conic residual
//      asserted below is what says that repair is right, since each chunk alone has a conic and
//      the four share one only after it.
V385 := ReadNelsonValues("data/nelson/values_6_385.m");
assert V385`D eq 6 and V385`N eq 385 and V385`g eq 3 and V385`cut eq 46 and V385`prec eq 90;
assert V385`W eq [1,5,6,11,14,21,30,55,66,70,105,154,231,330,770,1155];
m385 := eval (Read("data/models/models_6_385.m") cat "\nreturn models;");
e385 := m385[[Integers()|1,5,6,11,14,21,30,55,66,70,105,154,231,330,770,1155]][1];
assert e385[1] eq 3 and e385[3] eq 0;
f385 := e385[2];
assert &and[ Coefficient(f385, i) eq 0 : i in [1,3,5,7] ];
cg := [ Coefficient(f385, 2*j) : j in [0..4] ];
gA := cg[1]*cg[3]/cg[2]^2; gB := cg[2]*cg[4]/cg[3]^2; gC := cg[3]*cg[5]/cg[4]^2;
assert Denominator(gA) eq 1805524232 and Denominator(gB) eq 22477081 and Denominator(gC) eq 168200;
A385, B385, C385, um385, d385 := NelsonMarkedInvariants(V385);
assert um385 eq 2;
assert d385`conicres lt 10^(-13);                  // the conjugate columns agree across chunks
CC385 := Parent(A385);
rel385 := func< z, q | Abs(z - CC385 ! q)/Abs(CC385 ! q) >;
assert rel385(B385, gB) lt 10^(-9);
assert (rel385(A385, gA) lt 10^(-9) and rel385(C385, gC) lt 10^(-9))
    or (rel385(A385, gC) lt 10^(-9) and rel385(C385, gA) lt 10^(-9));
qeps385 := func< z, q | Denominator(q)^2 * 20 * Abs(Imaginary(z)) >;
assert qeps385(A385, gA) gt 100 and qeps385(B385, gB) gt 100;   // out of reach, not marginal
assert qeps385(C385, gC) lt 1/2;                                // the one the criterion accepts
