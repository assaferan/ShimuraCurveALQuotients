// SpecialFiberIsomorphism may rule out X_0(D,Np)/W only through a source X_0(D,N)/W' that is not
// subhyperelliptic over the algebraic closure of F_p, for the same p.  Not hyperelliptic over Q is
// not enough, and the source's own TestInWhichProved says nothing about p.
//
// Regression: the 2026-09-26 rerun ruled out these three genus-3 targets, each of which is
// hyperelliptic over Q (explicit models y^2 = f(x)), via a source that is non-hyperelliptic over Q
// but hyperelliptic mod p (the sources as recorded there, CurveIDs 230, 272, 270):
//   X_0(114)/<w2,w38>  <- X_0(57)/<w19>  at p = 2  (source proved by SpecialFiber at p = 19)
//   X_0(130)/<w2,w26>  <- X_0(65)/<w13>  at p = 2  (source proved by SpecialFiber at p = 5)
//   X_0(195)/<w5,w195> <- X_0(65)/<w5>   at p = 3  (source proved by SpecialFiber at p = 5)
// Here the sources are marked IsSubhyp = false with those labels, as in the rerun; the targets must
// stay undecided.  As a consistency check, a target recorded hyperelliptic must not raise the
// "proven hyperelliptic" error either, i.e. no certificate may be found for its source at p.
//
// Positive control, the paper's example: X_0(10,63)/<w10,w63> is isomorphic over F_7 to
// X_0(10,9)/<w10>, which has more than 2q+2 points over F_7 (the rerun's source 4027 was ruled by
// "Trace, with p^v = 7^1"), so it must still be ruled out, with that certificate in its label.
// Also a twisted-trace certificate at p, and a check that an undecided source is not used.
//
// The AL fixed-point certificate is valid at odd p only, and [FH] Prop 6 (ComplicatedAL) at none:
//  * X_0(89) (w89 has 12 fixed points, g = 7) is certified at p = 3 by the AL fixed-point test, so
//    X_0(267)/<w3> is ruled out through it; at p = 2 no certificate is found (the AL one is refused),
//    so X_0(178)/<w2> is not.
//  * X_0(10,43)/<w2,w5> (4405) is non-hyperelliptic over Q by ComplicatedAL (N2 = 430), which is its
//    only would-be certificate at p = 3.  The rerun used it to rule out X_0(10,129)/W (5040) at p = 3;
//    that must no longer happen.

mk := function(D, N, W, g, id)
    X := CreateShimuraQuot(D, N, W);
    X`g := g;
    X`CurveID := id;
    return X;
end function;

procedure test_regression()
    pairs := [<<1, 57, {1, 19}>, <1, 114, {1, 2, 19, 38}>, 2, "SpecialFiber p=19 case=2 component=X_0(3)/{ 1 }">,
              <<1, 65, {1, 13}>, <1, 130, {1, 2, 13, 26}>, 2, "SpecialFiber p=5 case=1 component=X_0(13)/{ 1, 13 }">,
              <<1, 65, {1, 5}>, <1, 195, {1, 5, 39, 195}>, 3, "SpecialFiber p=5 case=2 component=X_0(13)/{ 1 }">];
    for hyp in [false, true] do
        curves := [];
        for pr in pairs do
            src := mk(pr[1][1], pr[1][2], pr[1][3], 3, #curves + 1);
            assert GenusShimuraCurveQuotient(src`D, src`N, src`W) eq 3;
            src`IsSubhyp := false; src`IsHyp := false; src`TestInWhichProved := pr[4];
            Append(~curves, src);
            tgt := mk(pr[2][1], pr[2][2], pr[2][3], 3, #curves + 1);
            assert GenusShimuraCurveQuotient(tgt`D, tgt`N, tgt`W) eq 3;
            if hyp then
                tgt`IsSubhyp := true; tgt`IsHyp := true; tgt`TestInWhichProved := "explicit model";
            end if;
            Append(~curves, tgt);
            // the source has no certificate at p
            assert not NonHyperellipticAtPrimeCertificate(src, pr[3]);
        end for;
        SpecialFiberIsomorphism(~curves);   // with hyp: must not raise the consistency error
        for j in [2, 4, 6] do
            if hyp then
                assert curves[j]`IsHyp and curves[j]`TestInWhichProved eq "explicit model";
            else
                assert not assigned curves[j]`IsSubhyp;
            end if;
        end for;
    end for;
end procedure;

procedure test_positive_control()
    src := mk(10, 9, {1, 10}, 3, 1);
    assert GenusShimuraCurveQuotient(src`D, src`N, src`W) eq 3;
    src`IsSubhyp := false; src`IsHyp := false; src`TestInWhichProved := "Trace, with p^v = 7^1";
    tgt := mk(10, 63, {1, 10, 63, 630}, 7, 2);
    assert GenusShimuraCurveQuotient(tgt`D, tgt`N, tgt`W) eq 7;
    ok, cert := NonHyperellipticAtPrimeCertificate(src, 7);
    assert ok and cert eq "more than 2q+2 points over F_q, q = 7^1";
    curves := [src, tgt];
    SpecialFiberIsomorphism(~curves);
    assert assigned curves[2]`IsSubhyp and not curves[2]`IsSubhyp and not curves[2]`IsHyp;
    assert curves[2]`TestInWhichProved eq
        "SpecialFiberIsomorphism, isomorphic over F_7 to curve 1 (source: more than 2q+2 points over F_q, q = 7^1)";
    // A twisted certificate at p: in the rerun X_0(33,14)/<w6,w154> (8645) was ruled at p = 7 via
    // X_0(33,2)/<w6> (8478), itself ruled by "TwistedTrace, h = w22 with p^v = 7^1".
    s3 := mk(33, 2, {1, 6}, 3, 1);
    assert GenusShimuraCurveQuotient(s3`D, s3`N, s3`W) eq 3;
    s3`IsSubhyp := false; s3`IsHyp := false; s3`TestInWhichProved := "TwistedTrace, h = w22 with p^v = 7^1";
    t3 := mk(33, 14, {1, 6, 154, 231}, 3, 2);
    assert GenusShimuraCurveQuotient(t3`D, t3`N, t3`W) eq 3;
    curves := [s3, t3];
    SpecialFiberIsomorphism(~curves);
    assert curves[2]`TestInWhichProved eq
        "SpecialFiberIsomorphism, isomorphic over F_7 to curve 1 (source: TwistedTrace, h = w22 with p^v = 7^1)";
    // an unmarked source is never used, even though it would certify
    src2 := mk(10, 9, {1, 10}, 3, 1);
    assert GenusShimuraCurveQuotient(src2`D, src2`N, src2`W) eq 3;
    tgt2 := mk(10, 63, {1, 10, 63, 630}, 7, 2);
    curves := [src2, tgt2];
    SpecialFiberIsomorphism(~curves);
    assert not assigned curves[2]`IsSubhyp;
end procedure;

procedure test_al_fixed_points()
    src := mk(1, 89, {1}, 7, 1);
    assert GenusShimuraCurveQuotient(src`D, src`N, src`W) eq 7;
    src`IsSubhyp := false; src`IsHyp := false; src`TestInWhichProved := "ALFixedPointsOnQuotient, W_89 has 12 fixed points";
    ok, d, fix := TestALFixedPointsOnQuotient(src);
    assert (not ok) and d eq 89 and fix eq 12;
    ok, cert := NonHyperellipticAtPrimeCertificate(src, 3);
    assert ok and cert eq "ALFixedPointsOnQuotient, W_89 has 12 fixed points, p odd";
    assert not NonHyperellipticAtPrimeCertificate(src, 2);
    t2 := mk(1, 178, {1, 2}, GenusShimuraCurveQuotient(1, 178, {1, 2}), 2);
    t3 := mk(1, 267, {1, 3}, GenusShimuraCurveQuotient(1, 267, {1, 3}), 3);
    curves := [src, t2, t3];
    SpecialFiberIsomorphism(~curves);
    assert not assigned curves[2]`IsSubhyp;
    assert curves[3]`TestInWhichProved eq
        "SpecialFiberIsomorphism, isomorphic over F_3 to curve 1 (source: ALFixedPointsOnQuotient, W_89 has 12 fixed points, p odd)";
end procedure;

procedure test_complicated_not_used()
    src := mk(10, 43, {1, 2, 5, 10}, 3, 1);
    assert GenusShimuraCurveQuotient(src`D, src`N, src`W) eq 3;
    // over Q it is ComplicatedAL, and nothing simpler
    assert IsDefined(TestComplicatedALFixedPointsOnQuotient(10, 43), src`W);
    assert TestALFixedPointsOnQuotient(src);
    src`IsSubhyp := false; src`IsHyp := false;
    src`TestInWhichProved := "SpecialFiberD10 p=43 case=1 component=X_0(10,1)/{ 1, 2, 5, 10 }";
    assert not NonHyperellipticAtPrimeCertificate(src, 3);
    W := {1, 2, 5, 10, 129, 258, 645, 1290};
    tgt := mk(10, 129, W, GenusShimuraCurveQuotient(10, 129, W), 2);
    curves := [src, tgt];
    SpecialFiberIsomorphism(~curves);
    assert not assigned curves[2]`IsSubhyp;
end procedure;

// The certificate cache is keyed by <D, N, W, g, p>: g is an attribute that the certificate reads
// (the AL test compares fixed points with 2g+2, the point count runs v up to 4g^2) and does not
// recompute, so two curves that differ only in g must not share a cached answer.  X_0(89) at p = 3
// gives different answers for g = 5 (w89's 12 fixed points are 2g+2, so the AL test is refused,
// and the point count certifies at 3^2) and for its true g = 7 (the AL test).  Queried
// wrong-right-wrong, so a key without g fails whichever answer is cached first.
procedure test_certificate_cache_key()
    wrong := mk(1, 89, {1}, 5, 1);
    right := mk(1, 89, {1}, 7, 1);
    assert GenusShimuraCurveQuotient(right`D, right`N, right`W) eq 7;
    ok, cert := NonHyperellipticAtPrimeCertificate(wrong, 3);
    assert ok and cert eq "more than 2q+2 points over F_q, q = 3^2";
    ok, cert := NonHyperellipticAtPrimeCertificate(right, 3);
    assert ok and cert eq "ALFixedPointsOnQuotient, W_89 has 12 fixed points, p odd";
    ok, cert := NonHyperellipticAtPrimeCertificate(wrong, 3);
    assert ok and cert eq "more than 2q+2 points over F_q, q = 3^2";
end procedure;

// CheckTwistedAtPrime keeps an involution only if it commutes with T_l for every prime l of the
// twisted filters (twTracePrimes and twWeilPrimes) and p, not for p alone.  No example is known
// where that changes a verdict (for l prime to the level every admissible h commutes with T_l), so
// the list itself is pinned: the third return value is read back from the Hecke operators that
// twCurveData actually built.  X_0(33,2)/<w6>, g = 3, L = 66, p = 7: the trace primes are those
// below 4g^2 = 36 prime to 66, and the Weil-table primes for g = 3 (up to 23) add nothing new.
// Reverting to [p] makes this [7].
procedure test_twisted_prime_list()
    X := mk(33, 2, {1, 6}, 3, 1);
    assert GenusShimuraCurveQuotient(X`D, X`N, X`W) eq 3;
    ok, cert, used := CheckTwistedAtPrime(X, 7);
    assert (not ok) and cert eq "TwistedTrace, h = w22 with p^v = 7^1";
    assert used eq [l : l in PrimesUpTo(35) | 66 mod l ne 0];
    assert used eq [5, 7, 13, 17, 19, 23, 29, 31];
end procedure;

test_regression();
test_positive_control();
test_al_fixed_points();
test_complicated_not_used();
test_certificate_cache_key();
test_twisted_prime_list();
