// X_0^142(1)/w_71 is NOT hyperelliptic -- decided through the genus-one star, by three routes
// that agree, where the pipeline's filters had left the verdict open (data/*: the entry for
// W = {1, 71} carries no IsHyp) and point counts cannot decide it (every count up to p = 23 is
// within the hyperelliptic bound 2p + 2).
//
// The setting.  C = X_0^142(1)/w_71 has genus 3 and E = X_0^142(1)/W has genus 1 (E142A1,
// tests/GuoYangExample42.m).  C -> E is a double cover branched at the four points of E fixed by
// w_2 or w_142: P_-4, P_-8 (Ogg: w_2 fixes the CM points of discriminants -4, -8) and the two
// conjugate points of discriminant -568 = -4*142 (w_142).  If C were hyperelliptic, its
// hyperelliptic involution would commute with w_71 and descend to tau: P -> Q0 - P on E
// preserving the branch locus B, and tau would fix no branch point (a lift of tau to C would
// otherwise have at most seven fixed points instead of eight), so B would split into two pairs
// with equal sums Q0.  Conversely, of the four double covers branched at B, whose classes are
// Q0 + E[2], the one of class Q0 is hyperelliptic and the others are Ciani quartics; here no
// pairing has equal sums, so none of the four is hyperelliptic and the class is not needed.
//
// Guo-Yang's Assumption 40 (trivial character of the Borcherds forms) is not needed either:
// E(Q) is Z with no torsion, the divisor of a Borcherds product is rational and its class has
// finite order, so it is principal; the positions of P_-4, P_-8 and of R + R' follow from the
// Schofer values alone.
//
// Sources of the numbers.  The positions of the ten rational CM points as multiples of
// Q = (1, 1) are Guo-Yang, Compositio 153 (2017), Example 42, the second table of the example
// (-Q, -2Q, O, Q, 2Q, 3Q, 4Q, -3Q, -4Q, -6Q for -3, -4, -8, -19, -20, -24, -40, -43, -148,
// -232).  R + R' = 16 Q, the entry 4 - 16 of the obstruction functional at the cycle, the pole
// order 159 at which it is computed, and the twist c = -15006 below were computed by this
// repository and are regression values, not verifications.  That the field of definition of the
// -568 points is Q(sqrt -71) is checked both from the model and from the CM-field code.
//
// What is checked.
//   1. The model and the positions of the CM points (the paper's table).
//   2. A Borcherds form with divisor Z(-568) + P_-148 - 3 P_-40 (found by the prescribed-divisor
//      search) is pinned inside L(3 P_-40 - P_-148) by its Schofer values at eight CM points; its
//      remaining zeros are the two points R, R' of discriminant -568, defined over Q(sqrt -71),
//      and R + R' = 16 Q.  This is what the obstruction functional predicts (CMDivisorObstruction
//      read in point multiplicities assigns the cycle 2 + 2 - 16 = 4 - 16).
//   3. No pairing of B has equal sums: -2Q + O != 16Q, and a mixed pairing would need R - R' to
//      be 2-torsion.  So the criterion says: not hyperelliptic.
//   4. Independently: the double cover w^2 = phi, phi the generator of
//      L(2 P_-24 + 2 P_-40 - P_-4 - P_-8 - (R + R')), has genus 3 and its three holomorphic
//      differentials satisfy no quadratic relation (canonical image a plane quartic), while the
//      cover branched at the equal-sum locus {P_-4, P_-8, P_-19, P_-43} with the matching class
//      does satisfy one (hyperelliptic) and the criterion says so.
//   5. The cover is X/w_71: with c = -15006 = -2*3*41*61 (3, 41, 61 come from the scaling of
//      the generator phi), the curve w^2 = c phi has the trace formula's number of points of
//      X_0^142(1)/w_71 over F_p at every prime 5 <= p <= 250 where the reduction of the model is
//      a curve of genus 3, which is at least forty primes.

procedure test_bielliptic_142()
    printf "X_0^142(1)/w_71 through the genus-one star: branch pair, criterion, canonical quadric, point counts...";
    D := 142; N := 1;
    Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
    Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};
    X71 := [c : c in curves | c`W eq {1, 71}][1];
    assert X71`g eq 3;
    ds := [-3, -4, -8, -19, -20, -24, -40, -43, -148, -232];
    divisors := [<[-4, -8], -3>, <[-19, -43], -3>, <[-8, -40], -20>, <[-19, -24], -20>, <[-568, -148], -40>];
    fs := BorcherdsForms(star, curves : Divisors := divisors, PoolClassNumber := 4, IntegralSolution := true);
    vals := BorcherdsFormValues([fs[i] : i in [1..5]], ds, star);
    xs := CoordinateFromBorcherdsForms(vals[1], vals[2], ds, [-4, -8], [-19, -43]);
    ys := CoordinateFromBorcherdsForms(vals[3], vals[4], ds, [-8, -40], [-19, -24]);
    F := BidegreeTwoRelation(xs, ys);
    E, toE := EllipticCurveFromRelation(F, xs, ys : BasePoint := Index(ds, -8));
    assert CremonaReference(E) eq "142a1";
    Cp := Domain(toE);
    pts := AssociativeArray();
    for j->d in ds do
        if Type(xs[j]) ne Infty and Type(ys[j]) ne Infty then pts[d] := toE(Cp![xs[j], ys[j], 1]); end if;
    end for;
    pts[-3] := DivisionPoints(pts[-4] + pts[-8], 2)[1];
    pts[-20] := DivisionPoints(pts[-8] + pts[-40], 2)[1];
    Q := E![1, 1, 1];
    mult := func<P | [m : m in [-20..20] | m*Q eq P][1]>;
    // (1) Guo-Yang's positions, up to the sign of the generator
    paper := [-1, -2, 0, 1, 2, 3, 4, -3, -4, -6];
    ours := [mult(pts[d]) : d in ds];
    assert ours eq paper or ours eq [-m : m in paper];
    sgn := (ours eq paper) select 1 else -1;

    // (2) the branch pair of discriminant -568
    FE := FunctionField(E);
    pl := AssociativeArray();
    for d in ds do pl[d] := Place(pts[d]); end for;
    V, mV := RiemannRochSpace(3*pl[-40] - pl[-148]);
    assert Dimension(V) eq 2;
    f1 := mV(V.1); f2 := mV(V.2);
    v5 := vals[5];
    ev := func<f, d | Evaluate(f, pl[d])>;
    good := [d : d in ds | Type(v5[Index(ds, d)]) ne Infty and v5[Index(ds, d)] ne 0 and d notin [-40, -148]];
    assert #good ge 4;
    d1 := good[1]; d2 := good[2];
    r := v5[Index(ds, d1)] / v5[Index(ds, d2)];
    psi5s := [];
    for s in [1, -1] do
        den := ev(f2, d1) - s*r*ev(f2, d2);
        if den eq 0 then continue; end if;
        b := (s*r*ev(f1, d2) - ev(f1, d1))/den;
        g := f1 + b*f2;
        kappa := v5[Index(ds, d1)] / AbsoluteValue(ev(g, d1));
        if forall{d : d in good | kappa*AbsoluteValue(ev(g, d)) eq v5[Index(ds, d)]} then Append(~psi5s, g); end if;
    end for;
    assert #psi5s eq 1;
    supp, mults := Support(Divisor(psi5s[1]));
    deg2 := [supp[i] : i in [1..#supp] | mults[i] eq 1 and Degree(supp[i]) eq 2];
    assert #deg2 eq 1;
    Rpl := deg2[1];
    K := ResidueClassField(Rpl);
    assert IsIsomorphic(K, QuadraticField(-71));
    assert IsIsomorphic(FieldOfDefinitionOfCMPoint(star, -568), QuadraticField(-71));
    EK := BaseChange(E, K);
    R := EK!Eltseq(RepresentativePoint(Rpl));
    sigma := Automorphisms(K)[2];
    Rc := EK![sigma(c) : c in Eltseq(R)];
    S := R + Rc;
    assert forall{c : c in Eltseq(S) | c in Rationals()};
    S := E![Rationals()!c : c in Eltseq(S)];
    assert S eq 16*sgn*Q;     // regression value (computed here)
    // the functional, read on the cycle: its entry at Z(-568) in point multiplicities is 4 - 16
    // (regression values: the entry and the pole order 159)
    Kf, _ := CMDivisorObstruction(star, ds cat [-568] : PoleOrder := 159);
    assert Dimension(Kf) eq 1;
    v := Eltseq(Basis(Kf)[1]);
    scale := v[Index(ds, -19)] / (2 - 1);
    assert v[#v]/scale eq 4 - 16;

    // (3) the criterion
    B := [EK!pts[-4], EK!pts[-8], R, Rc];
    Ls := DivisionPoints(E![Rationals()!c : c in Eltseq(&+B)], 2);
    assert #Ls eq 1;       // E(Q) torsion-free: one rational square-root class
    hyp, equal := BiellipticBranchCriterion(B, Ls[1]);
    assert (not hyp) and #equal eq 0;

    // (4) the explicit cover and the canonical quadric, with its control
    phiV, mphi := RiemannRochSpace(2*pl[-24] + 2*pl[-40] - pl[-4] - pl[-8] - Rpl);
    assert Dimension(phiV) eq 1;
    phi := mphi(phiV.1);
    AFF, toAFF := AlgorithmicFunctionField(FE);
    PT<T> := PolynomialRing(AFF);
    FC := FunctionField(T^2 - toAFF(phi));
    assert Genus(FC) eq 3;
    assert not IsHyperellipticByCanonicalQuadric(FC);
    ctlV, mctl := RiemannRochSpace(pl[-4] + pl[-8] - pl[-19] - pl[-43]);
    assert Dimension(ctlV) eq 1;
    FCtl := FunctionField(T^2 - toAFF(mctl(ctlV.1)));
    assert Genus(FCtl) eq 3;
    assert IsHyperellipticByCanonicalQuadric(FCtl);
    assert BiellipticBranchCriterion([pts[-4], pts[-8], pts[-19], pts[-43]], pts[-4] + pts[-8]);

    // (5) the cover is X/w_71: point counts with the twist c = -15006 (regression value)
    c := -15006;
    pf := ProjectiveFunction(phi);
    Nphi := Numerator(pf); Dphi := Denominator(pf);
    cN := Coefficients(Nphi); cD := Coefficients(Dphi);
    den := LCM([Denominator(x) : x in cN cat cD]);
    checked := 0;
    for p in [q : q in PrimesUpTo(250) | q ge 5] do
        Ep := ChangeRing(E, GF(p));
        FEp := FunctionField(Ep);
        AFFp, toAFFp := AlgorithmicFunctionField(FEp);
        R3 := PolynomialRing(GF(p), 3);
        xp := FEp.1; yp := FEp.2;
        Np := &+[R3 | (GF(p)!(den*MonomialCoefficient(Nphi, m)))*Monomial(R3, Exponents(m)) : m in Monomials(Nphi)];
        Dp := &+[R3 | (GF(p)!(den*MonomialCoefficient(Dphi, m)))*Monomial(R3, Exponents(m)) : m in Monomials(Dphi)];
        if Evaluate(Dp, [xp, yp, 1]) eq 0 or Evaluate(Np, [xp, yp, 1]) eq 0 then continue; end if;
        phip := Evaluate(Np, [xp, yp, 1]) / Evaluate(Dp, [xp, yp, 1]);
        PTp<Tp> := PolynomialRing(AFFp);
        try
            Fc := FunctionField(Tp^2 - (GF(p)!c)*toAFFp(phip));
            gp := Genus(Fc);
        catch e
            continue;
        end try;
        if gp ne 3 then continue; end if;
        assert #Places(Fc, 1) eq ComputePointsViaTrace(X71, p, 1);
        checked +:= 1;
    end for;
    printf " [%o primes]", checked;
    assert checked ge 40;
    printf "Done!\n";
end procedure;

test_bielliptic_142();
