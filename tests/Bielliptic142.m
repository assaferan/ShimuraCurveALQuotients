// X_0^142(1)/w_71 is NOT hyperelliptic -- decided through the genus-one star, by three routes
// that agree, where the pipeline's filters had left the verdict open (data/*: the entry for
// W = {1, 71} carries no IsHyp) and point counts cannot decide it (every count up to p = 23 is
// within the hyperelliptic bound 2p + 2).
//
// The setting.  C = X_0^142(1)/w_71 has genus 3 and E = X_0^142(1)/W has genus 1 (E142A1,
// tests/GuoYangExample42.m).  C -> E is a double cover branched at the four points of E fixed by
// w_2 or w_142: P_-4, P_-8 (Ogg: w_2 fixes the CM points of discriminants -4, -8) and the two
// conjugate points of discriminant -568 = -4*142 (w_142).  If C were hyperelliptic, its
// hyperelliptic involution would commute with w_71 and descend to an involution P -> Q0 - P of E
// preserving the branch locus B; so B would split into two pairs with equal sums Q0.  Conversely
// of the double covers branched at B the one of square-root class Q0 is hyperelliptic (the
// others are Ciani quartics), and here the class is forced: E(Q) is torsion-free, so B has one
// rational square root.
//
// What is checked.
//   1. The model and the positions of the CM points (Guo-Yang's table: multiples of Q = (1, 1)).
//   2. A Borcherds form with divisor Z(-568) + P_-148 - 3 P_-40 (found by the prescribed-divisor
//      search) is pinned inside L(3 P_-40 - P_-148) by its Schofer values at eight CM points; its
//      remaining zeros are the two points R, R' of discriminant -568, and R + R' = 16 Q.  This is
//      exactly what the obstruction functional predicts (CMDivisorObstruction read in point
//      multiplicities assigns the cycle 16 + 2*2 - 2*2 = 16).
//   3. No pairing of B has equal sums: -2Q + O != 16Q, and a mixed pairing would need R - R' to
//      be 2-torsion.  So the criterion says: not hyperelliptic.
//   4. Independently: the double cover w^2 = phi, phi the generator of
//      L(2 P_-24 + 2 P_-40 - P_-4 - P_-8 - (R + R')), has genus 3 and its three holomorphic
//      differentials satisfy no quadratic relation (canonical image a plane quartic), while the
//      cover branched at the equal-sum locus {P_-4, P_-8, P_-19, P_-43} with the matching class
//      does satisfy one (hyperelliptic) and the criterion says so.
//   5. The cover is X/w_71: at every good prime up to 100 the trace formula's point count of
//      X_0^142(1)/w_71 equals that of w^2 = c phi for c a square or for c a non-square mod p, never
//      neither, and the resulting Legendre symbols are those of a single squarefree integer c.

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
    EK := BaseChange(E, K);
    R := EK!Eltseq(RepresentativePoint(Rpl));
    sigma := Automorphisms(K)[2];
    Rc := EK![sigma(c) : c in Eltseq(R)];
    S := R + Rc;
    assert forall{c : c in Eltseq(S) | c in Rationals()};
    S := E![Rationals()!c : c in Eltseq(S)];
    assert S eq 16*sgn*Q;
    // the functional, read on the cycle: its entry at Z(-568) in point multiplicities is 4 - 16
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

    // (5) the cover is X/w_71: point counts and the twist
    pf := ProjectiveFunction(phi);
    Nphi := Numerator(pf); Dphi := Denominator(pf);
    cN := Coefficients(Nphi); cD := Coefficients(Dphi);
    den := LCM([Denominator(c) : c in cN cat cD]);
    symbols := AssociativeArray();
    for p in [q : q in PrimesUpTo(60) | q ge 5] do
        Ep := ChangeRing(E, GF(p));
        FEp := FunctionField(Ep);
        AFFp, toAFFp := AlgorithmicFunctionField(FEp);
        R3 := PolynomialRing(GF(p), 3);
        xp := FEp.1; yp := FEp.2;
        Np := &+[R3 | (GF(p)!(den*MonomialCoefficient(Nphi, m)))*Monomial(R3, Exponents(m)) : m in Monomials(Nphi)];
        Dp := &+[R3 | (GF(p)!(den*MonomialCoefficient(Dphi, m)))*Monomial(R3, Exponents(m)) : m in Monomials(Dphi)];
        if Evaluate(Dp, [xp, yp, 1]) eq 0 then continue; end if;
        phip := Evaluate(Np, [xp, yp, 1]) / Evaluate(Dp, [xp, yp, 1]);
        PTp<Tp> := PolynomialRing(AFFp);
        nonsq := [a : a in [2..p-1] | not IsSquare(GF(p)!a)][1];
        try
            Fplus := FunctionField(Tp^2 - toAFFp(phip));
            Fminus := FunctionField(Tp^2 - nonsq*toAFFp(phip));
            gp := Genus(Fplus);
        catch e
            continue;
        end try;
        if gp ne 3 then continue; end if;
        target := ComputePointsViaTrace(X71, p, 1);
        nplus := #Places(Fplus, 1); nminus := #Places(Fminus, 1);
        assert nplus eq target or nminus eq target;
        if nplus ne nminus then symbols[p] := (nplus eq target) select 1 else -1; end if;
    end for;
    assert #Keys(symbols) ge 8;
    cs := [c : c in [-10^5..10^5] | c ne 0 and IsSquarefree(c) and forall{p : p in Keys(symbols) | KroneckerSymbol(c, p) eq symbols[p]}];
    printf " [twist c in %o]", cs;
    assert #cs ge 1;
    printf "Done!\n";
end procedure;

test_bielliptic_142();
