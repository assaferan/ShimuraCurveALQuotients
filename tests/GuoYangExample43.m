// EXTERNAL SOURCE: Guo-Yang, "Equations of hyperelliptic Shimura curves", Compositio Math. 153
// (2017), Example 43, pp. 31-32 of the journal version.  Everything compared against below is
// printed there: the four divisors, the table of x and y at ten CM points, the relation between x
// and y, the Weierstrass model y^2 = x^6 - 18x^4 + 113x^2 - 32, and the two elliptic curves
// E_302A1 and E_302C1 whose product the Jacobian is isogenous to.
//
// X = X_0^302(1)/W_{302,1} has genus two, so neither the hauptmodul search nor the genus-one recipe
// of Example 42 reaches it.  The paper builds the model from four Borcherds forms with the divisors
//     div psi_1 = P_-43 + P_-72 - P_-19 - P_-88,    div psi_2 = P_-20 + P_-36 - P_-19 - P_-88,
//     div psi_3 = P_-8 + 2 P_-40 - P_-4 - 2 P_-88,  div psi_4 = P_-11 + P_-19 + P_-43 - P_-4 - 2 P_-88,
// takes x with div x = div psi_1 and y with div y = div psi_3, and fits the relation between them.
// Two features of this divisor shape are new here: the poles are two distinct points rather than
// one double point, and psi_3, psi_4 have a point of multiplicity two.
//
// WHICH FORMS ARE USED HERE.  The zeros of psi_1 and psi_2 include the discriminants -72 and -36,
// which are not fundamental (both have conductor 3), and this repository's CM pool is the
// fundamental discriminants of bounded class number, so those two points are not available.  In
// their place the test uses the forms of two other fibres of the same function, read off the
// paper's own table: x = 1 at -11 and -40, and x = 5/3 at -148 and -232, over the same poles.  The
// coordinate u they define is therefore the paper's x composed with the affine map that sends 1 to
// 0 and 5/3 to 1, namely u = (3x - 3)/2, and the curve is the same.  Guo-Yang's psi_3 and psi_4 are
// used as printed.
//
// ASSUMPTION 40 IS A GENUINE HYPOTHESIS HERE, as the paper states.  A Borcherds form with integral
// principal part has a divisor class of finite order, so when the Jacobian has no rational torsion
// the form is a function and nothing need be assumed; but the Jacobian of this star is isogenous to
// E_302A1 x E_302C1 and E_302A1 has a rational point of order five (the paper's P = (-32, 256)),
// which part 2 checks.  So the model below is conditional on the assumption, exactly as the
// paper's is.
//
// WHAT CERTIFIES THE MODEL.  The relation is not determined by the values at the CM points alone --
// there are ten of them where eleven conditions would be needed -- and the divisor of y supplies
// the rest, so a wrong divisor gives a wrong relation rather than no relation: part 6 exhibits one,
// a different curve of the same genus through the same ten points.  What settles the matter is the
// trace formula, which counts the points of the star over F_p with no model at all: the curve the
// paper prints agrees with it at all sixteen good primes below 60, and the curve from the wrong
// divisor agrees at none of the thirteen where it has good reduction.

procedure test_guoyang_example43()
    printf "Guo-Yang Example 43: X_0^302(1)* (genus two) from Borcherds forms of degree two and three...";
    D := 302; N := 1;
    Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
    Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
    assert Xstar`g eq 2;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};

    // (1) The ten rational CM points of the paper's table, and no others up to class number 2.
    ds := [-4, -8, -11, -19, -20, -40, -43, -88, -148, -232];
    rat, _ := RationalandQuadraticCMPoints(star : bd := 2, coprime_to_level := false);
    assert Sort([p[1] : p in rat]) eq Sort(ds);

    Rt<t> := PolynomialRing(Rationals());
    GYcurve := HyperellipticCurve(t^6 - 18*t^4 + 113*t^2 - 32);
    E302A1 := EllipticCurve([1, 1, 1, -230, 1251]);
    E302C1 := EllipticCurve([1, -1, 1, 0, 3]);
    // the number of points of a hyperelliptic curve over F_p, and -1 where it does not reduce to one
    count := function(f, p)
        if exists{c : c in Coefficients(f) | Valuation(Denominator(c), p) gt 0} then return -1; end if;
        fp := PolynomialRing(GF(p))![GF(p) | c : c in Coefficients(f)];
        if Degree(fp) lt 2*Genus(GYcurve) + 1 or not IsSquarefree(fp) then return -1; end if;
        return #Points(HyperellipticCurve(fp));
    end function;
    good := [p : p in PrimesInInterval(3, 60) | (D*N) mod p ne 0];
    traces := [ComputePointsViaTrace(Xstar, p, 1) : p in good];

    // (2) The paper's table alone, with no Borcherds form computed: its values of x and y and the
    // divisor of its y give back its relation, and completing the square gives its curve.  This
    // checks the fit and the model against the paper independently of anything found below.
    paper_x := [* -1, 3/2, 1, Infinity(), 2, 1, 0, Infinity(), 5/3, 5/3 *];
    paper_y := [* Infinity(), 0, 1, 1, -1, 0, 1, Infinity(), -1/9, -1/2 *];
    divy := [<-8, 1>, <-40, 2>, <-4, -1>, <-88, -2>];
    Fpaper := HyperellipticRelation(ds, paper_x, paper_y, divy);
    PP<xx, yy> := Parent(Fpaper);
    // the fit returns the relation with a positive leading coefficient, so up to the overall sign
    printed := (xx+1)*yy^2 + (-2*xx^3 + 11*xx^2 - 13*xx + 2)*yy + (2*xx-3)*(xx-1)^2;
    assert Fpaper eq printed or Fpaper eq -printed;
    Hpaper := HyperellipticModelFromRelation(Fpaper);
    assert Genus(Hpaper) eq 2;
    assert IsIsomorphic(Hpaper, GYcurve);
    // the Jacobian is the product of the paper's two elliptic curves, and one of them has a
    // rational point of order five, so Assumption 40 is not automatic (see the header)
    assert CremonaReference(E302A1) eq "302a1" and CremonaReference(E302C1) eq "302c1";
    assert #TorsionSubgroup(E302A1) eq 5;
    Jp := Jacobian(Hpaper);
    order_mod := func<E, p | #EllipticCurve([GF(p) | a : a in aInvariants(E)])>;
    assert forall{p : p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 37] |
                  #Points(BaseChange(Jp, GF(p))) eq order_mod(E302A1, p) * order_mod(E302C1, p)};
    // the paper's equation against the trace formula, which needs no model: every good prime below
    // 60 agrees, and the only primes of bad reduction are the two dividing the level
    assert Set(BadPrimes(GYcurve)) eq {2, 151};
    nagree := 0;
    for j->p in good do
        c := count(t^6 - 18*t^4 + 113*t^2 - 32, p);
        if c lt 0 then continue; end if;
        assert c eq traces[j];
        nagree +:= 1;
    end for;
    assert nagree ge 15;

    // (3) Our four forms.  The first two replace the paper's psi_1 and psi_2 (see the header); the
    // last two are the paper's psi_3 and psi_4.  The principal parts must be integral, so that each
    // Borcherds form is a function and not a root of one.
    divisors := [[<-11, 1>, <-40, 1>, <-19, -1>, <-88, -1>],
                 [<-148, 1>, <-232, 1>, <-19, -1>, <-88, -1>],
                 [<-8, 1>, <-40, 2>, <-4, -1>, <-88, -2>],
                 [<-11, 1>, <-19, 1>, <-43, 1>, <-4, -1>, <-88, -2>]];
    fs := BorcherdsForms(star, curves : Divisors := divisors, IntegralSolution := true);
    assert Keys(fs) eq {1, 2, 3, 4};
    for i in [1..4] do
        assert Set(DivisorOfBorcherdsForm(fs[i], star)) eq {<c[1], c[2]> : c in divisors[i]};
        g := qExpansionAtoo(fs[i], 1);
        assert forall{m : m in [Valuation(g)..0] | IsIntegral(Coefficient(g, m))};
    end for;

    // (4) Schofer values at the ten points, the two coordinates, and the paper's table entry by
    // entry.  Our u is the paper's x composed with s -> (3s - 3)/2 (see the header); our y is
    // theirs, since both are normalised by the same divisor and the same value 1 at the zeros of
    // psi_4.  Every disagreement is printed before the assertion fires.
    vals := BorcherdsFormValues([fs[i] : i in [1..4]], ds, star);
    us := CoordinateFromBorcherdsForms(vals[1], vals[2], ds, [-11, -40], [-148, -232]);
    ys := CoordinateFromBorcherdsForms(vals[3], vals[4], ds, [-8, -40], [-11, -19, -43]);
    same := func<a, b | (Type(a) eq Infty and Type(b) eq Infty) or
                        (Type(a) ne Infty and Type(b) ne Infty and a eq b)>;
    moebius := func<s | Type(s) eq Infty select Infinity() else (3*s - 3)/2>;
    bad := [* *];
    for j->d in ds do
        if not same(us[j], moebius(paper_x[j])) then Append(~bad, <"x", d, us[j], moebius(paper_x[j])>); end if;
        if not same(ys[j], paper_y[j]) then Append(~bad, <"y", d, ys[j], paper_y[j]>); end if;
    end for;
    if #bad gt 0 then
        printf "\n  entries disagreeing with the paper's table (coordinate, d, ours, paper): %o\n", bad;
    end if;
    assert #bad eq 0;

    // (5) The relation between u and y, and the curve: the paper's relation in the coordinate u,
    // and a model isomorphic to the one the paper prints.
    F := HyperellipticRelation(ds, us, ys, divy);
    P<u, y> := Parent(F);
    subst := Evaluate(Fpaper, [(2*u + 3)/3, y]);
    subst *:= LCM([Denominator(c) : c in Coefficients(subst)]);
    subst /:= GCD([Integers() | c : c in Coefficients(subst)]);
    assert F eq subst or F eq -subst;
    H, num, den := HyperellipticModelFromRelation(F);
    assert Genus(H) eq 2;
    assert IsIsomorphic(H, GYcurve);
    // the substitution really is a point of the model: the hyperelliptic coordinate at the CM
    // points where both of ours are finite
    for j->d in ds do
        if Type(us[j]) eq Infty or Type(ys[j]) eq Infty then continue; end if;
        Y := Evaluate(num, [us[j], ys[j]]) / Evaluate(den, [us[j], ys[j]]);
        assert Y^2 eq Evaluate(HyperellipticPolynomials(H), us[j]);
    end for;
    // and the model against the trace formula
    fH := HyperellipticPolynomials(IntegralModel(H));
    nagree := 0;
    for j->p in good do
        c := count(fH, p);
        if c lt 0 then continue; end if;
        assert c eq traces[j];
        nagree +:= 1;
    end for;
    assert nagree ge 10;

    // (6) NEGATIVE CONTROL.  The ten points do not determine the relation by themselves, so the
    // conditions the divisor of y contributes have to be doing the work.  Exchanging the
    // multiplicities of its two zeros leaves y vanishing and polar at the same points, so no value
    // can detect it; the fit then returns a relation that passes through all ten points and whose
    // curve has genus two as well, but is not this one, and the trace formula refuses it at every
    // good prime.  A control that produced the same curve, or no curve, would not show that.
    // (two relations from different calls are built in different polynomial rings, which Magma
    // will not compare, so the control is stated on the curves they define)
    Fbad := HyperellipticRelation(ds, us, ys, [<-8, 2>, <-40, 1>, <-4, -1>, <-88, -2>]);
    Hbad := HyperellipticModelFromRelation(Fbad);
    assert Genus(Hbad) eq 2;
    assert not IsIsomorphic(Hbad, GYcurve);
    fbad := HyperellipticPolynomials(IntegralModel(Hbad));
    nref := 0;
    for j->p in good do
        c := count(fbad, p);
        if c lt 0 then continue; end if;
        assert c ne traces[j];
        nref +:= 1;
    end for;
    assert nref ge 10;
    printf " [refused at %o primes]Done!\n", nref;
end procedure;

test_guoyang_example43();
