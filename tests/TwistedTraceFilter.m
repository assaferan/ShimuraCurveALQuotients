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
// 5124 and 8387 are recorded non-hyperelliptic by GeneralizedComplicatedFixedPoints; they are tested
// here on fresh copies.  Controls: recorded-hyperelliptic curves with AL, V2 and V3 twists, which
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
    // decided curves are skipped (ShimuraQuot is a reference type, so compare against a saved string)
    proof := curves[8387]`TestInWhichProved;
    cs := [curves[8387]];
    FilterByTwistedTrace(~cs);
    assert (not cs[1]`IsSubhyp) and cs[1]`TestInWhichProved eq proof;
end procedure;

curves := GetHyperellipticCandidates();
test_TwistedTrace(curves);
