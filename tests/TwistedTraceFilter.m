// FilterByTwistedTrace: Tr((T_(p^v) - p T_(p^(v-2))) o h) < -(q+1) for an involution h defined over Q
// proves non-hyperellipticity.
//
// Known values, from the twisted-trace sweep (sweeps/twisted_trace/twist.m on branch
// twisted-trace-sweep; the (26,45) and (10,153) hits were reproduced by an independent reviewer
// implementation):
//   11719 X_0^142(1)/<w71>      q = 3, h = w2,     tr = -6
//   5124  (10,153)/<2,9,85>     q = 7, h = V3*w5,  tr = -12
//   8387  (26,45)/<5,9,26>      q = 7, h = V3,     tr = -12
//   896   (1,198)/<9,11>        q = 5, h = V3,     tr = -10   (V3 at p = 2 mod 3: needs 9 in W)
// All four are tested on fresh copies, so the test does not depend on their recorded status: on
// main 5124 and 8387 are recorded non-hyperelliptic by GeneralizedComplicatedFixedPoints, and once
// the GC V3-certificate fix (PR #46) lands they are undecided and become TwistedTrace rule-outs.  Controls: recorded-hyperelliptic curves with AL, V2 and V3 twists, which
// must not be ruled out (the sweep found 0 violations on all 914 hyperelliptic controls).

fresh := function(X)
    Y := CreateShimuraQuot(X`D, X`N, X`W);
    Y`g := X`g;
    Y`CurveID := X`CurveID;
    return Y;
end function;

procedure test_TwistedTrace(curves)
    known := [<11719, "w2", 3, 1, -6>, <896, "V3", 5, 1, -10>,
              <8387, "V3", 7, 1, -12>, <5124, "V3*w5", 7, 1, -12>];
    // the trace itself, on one curve
    ok, h, p, v, tr := CheckTwistedTrace(curves[11719]);
    assert not ok and <h, p, v, tr> eq <"w2", 3, 1, -6>;
    // hyperelliptic controls: (6,23)/<2> (AL twists), (1,63)/<9> and (1,126)/<7,9> (V3 with 9 in W),
    // (1,120)/<5,24> (V2)
    ctl := [1577, 263, 613, 583];
    assert &and[curves[c]`IsSubhyp : c in ctl];
    cs := [fresh(curves[k[1]]) : k in known] cat [fresh(curves[c]) : c in ctl];
    FilterByTwistedTrace(~cs);
    for i->k in known do
        X := cs[i];
        printf "  %o: %o\n", X`CurveID, assigned X`TestInWhichProved select X`TestInWhichProved else "-";
        assert assigned X`IsSubhyp and not X`IsSubhyp and not X`IsHyp;
        assert X`TestInWhichProved eq Sprintf("TwistedTrace, h = %o with p^v = %o^%o", k[2], k[3], k[4]);
    end for;
    assert &and[not assigned cs[i]`IsSubhyp : i in [#known+1..#cs]];
    // decided curves are skipped: 222 = X_0(56) (has V2 twists; decided by ALFixedPointsOnQuotient on
    // main and on integration).  ShimuraQuot is a reference type, so compare against a saved string.
    assert assigned curves[222]`IsSubhyp and not curves[222]`IsSubhyp;
    proof := curves[222]`TestInWhichProved;
    cs := [curves[222]];
    FilterByTwistedTrace(~cs);
    assert (not cs[1]`IsSubhyp) and cs[1]`TestInWhichProved eq proof;
end procedure;

// V3 is used only when 9 in W (sigma(V3) = V3 W9): on X_0(90)/<w9, w2> (420) the V3 twists are
// there, on X_0(90)/<w10> they must not be: there V3 preserves the w10-fixed part and is a
// nontrivial involution on it (checked by hand), but it is not defined over Q, so only the
// 9-in-W rule keeps it out.
procedure test_V3Guard(curves)
    names := func<X | [x[1] : x in TwistedWeilPolynomials(X, 7)]>;
    X := curves[420];
    assert <X`D, X`N> eq <1, 90> and 9 in X`W;
    assert "V3" in names(X);
    assert exists(Y){Y : Y in curves | Y`D eq 1 and Y`N eq 90 and Y`W eq {1, 10}};
    assert 9 notin Y`W;
    printf "  (1,90)/<w10> twists: %o\n", names(Y);
    assert not exists{h : h in names(Y) | "V3" in h};
end procedure;

// Star curves (W the full AL group; FilterByTwistedTraceStar): only V2, V3, V2 V3 can occur.
// X_0^*(396) (1272, undecided) is ruled out by V3 at q = 5, trace -10; the hyperelliptic star
// curves X_0^*(171) (813) and X_0^(9)(55)^* (10149), both with V3, are not.  A star curve with no
// V2/V3 has nothing to test and costs nothing (its level's modular symbols are never built).
procedure test_TwistedTraceStar(curves)
    X := curves[1272];
    assert <X`D, X`N> eq <1, 396> and IsStarCurve(X);
    ctl := [813, 10149];
    assert &and[IsStarCurve(curves[c]) and curves[c]`IsSubhyp : c in ctl];
    cs := [fresh(X)] cat [fresh(curves[c]) : c in ctl];
    FilterByTwistedTrace(~cs);
    printf "  star %o: %o\n", X`CurveID, cs[1]`TestInWhichProved;
    assert (not cs[1]`IsSubhyp) and cs[1]`TestInWhichProved eq "TwistedTrace, h = V3 with p^v = 5^1";
    assert &and[not assigned cs[i]`IsSubhyp : i in [2..#cs]];
    assert exists(Y){Y : Y in curves | IsStarCurve(Y) and Y`g ge 3 and Y`N mod 8 ne 0
                                       and Valuation(Y`N, 3) ne 2};
    Z := fresh(Y);
    assert CurveCostProxy(Z, "FilterByTwistedTraceStar") eq 0;
    assert CurveCostProxy(fresh(X), "FilterByTwistedTraceStar") gt 0;
end procedure;

curves := GetHyperellipticCandidates();
test_TwistedTrace(curves);
test_V3Guard(curves);
test_TwistedTraceStar(curves);
