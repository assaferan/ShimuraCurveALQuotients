// Regression test for the cusp-0 (SchoferFormula0) isometry-invariance fix, anchored to ground truth.
//
// The norm of a Borcherds form at a CM point is an isometry invariant of the CM point, so it must NOT
// depend on which representative of the (Magma-randomized) quaternion order ShimuraCurveLattice returns.
// Before the fix, SchoferFormula0 bucketed discriminant-group elements by an integer-mod-M residue of a
// NON-canonical scaled dual representative (normalized by the level M instead of the invariant scaling
// denom), so it dropped valid gammas whenever denom != M -- making the cusp-0 value depend on the order
// representative. This test uses X0^15(2) (=15_2), the CM-starved even-level base where Magma's
// randomized MaximalOrder genuinely yields distinct representatives, and checks against Guo-Yang
// "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), Table 45 ("CM-values of X^15_0(2)").
//
// Table 45's star hauptmodule s has 0 at disc -40, oo at -12, 1 at -52. Our fs[-1] shares the 0 and oo
// but is normalized to 1 at -120 (our s~ zero), where Table 45 gives s(-120)=2; hence our_s = s_GY / 2.
// We compare |our_s(d)| = raw_s(d) / raw_s(-120) to |s_GY(d) / 2|.

procedure test_Schofer_15_2_table45()
    printf "Testing Schofer cusp-0 fix on X0^15(2) vs Guo-Yang Table 45...";
    D := 15; N := 2;
    Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
    Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};
    Ldata := ShimuraCurveLattice(D, N);
    fs := BorcherdsForms(star, curves : Prec := 100);
    s := fs[-1];

    lam := func<Ld, d | ElementsOfNorm(Ld`Q, [-d], Ld`O, Ld`basis_L)[-d]>;
    raw := func<Ld, d | RationalNumber(AbsoluteValuesAtRationalCMPoint([s], d, star, Ld : Lambda := lam(Ld, d))[1])>;
    normval := func<Ld, d | AbsoluteValue(raw(Ld, d)) / AbsoluteValue(raw(Ld, -120))>;

    // (1) CORRECTNESS vs Table 45, s-column / 2. Includes the odd-fundamental (2-split) discs
    // -7,-15,-60: their level-prime (2|N) contribution -- the outer m=0 term kappa_0(0)'s N-part,
    // sum_{p|N/(N,d_fund)} log p (Yifan Yang arXiv:1503.07971 Lemma 20) -- is now restored, fixing the
    // former +4 log 2 (factor 16) discrepancy at 2-unramified CM discs.
    table45_half := [ <-7,1/4>, <-15,5/4>, <-60,-1/12>,
                      <-52,1>, <-88,4>, <-120,2>, <-132,-1>, <-148,1/25>, <-168,2/3>,
                      <-228,-1/9>, <-232,144/121>, <-280,10>, <-312,2/25>, <-340,9/17>,
                      <-372,-31/9>, <-408,68/25>, <-520,-8/121>, <-708,-841/121>, <-760,450/529> ];
    for t in table45_half do
        d, sGY := Explode(t);
        assert normval(Ldata, d) eq AbsoluteValue(sGY / 2);
    end for;

    // (2) ISOMETRY-INVARIANCE across a genuinely different order representative. This is the property
    // the fix restores; it FAILS on the pre-fix code. We check the odd-fundamental discs -7,-15,-60,
    // which is exactly where the pre-fix cusp-0 value was representative-dependent.
    found := false; Ld2 := Ldata;
    for tries in [1..80] do
        _ := [Random(1, 10^6) : i in [1..(tries mod 11) + 1]];
        cand := ShimuraCurveLattice(D, N);
        if Eltseq(cand`Q) ne Eltseq(Ldata`Q) then found := true; Ld2 := cand; break; end if;
    end for;
    assert found;
    for d in [-7, -15, -60, -52, -148] do
        assert normval(Ldata, d) eq normval(Ld2, d);
    end for;

    printf "Done!\n";
end procedure;

test_Schofer_15_2_table45();
