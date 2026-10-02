// CurveCostProxy for the Weil-polynomial stage must put the measured-heavy curves first.  The
// estimate is #W * sum p^(g/2) over the good primes the stage will use, which is the number of
// trace-formula terms the dominant call needs (see Utils.m).
//
// Every assertion below is an ORDERING of two curves whose times were measured, never a cost.
// Measured on lovelace with the class-number tables, over 34 curves, on main at 95b19e6 (with
// #56, #57 and #58), by vvdata/weyl-campaign/weil-retime-2026-10-02/weil_timing.m on the
// m0-theta-campaign branch (logs and the fit alongside it):
//
//     curve  shape                     g  #W   measured
//      1436  X_0(690)/W16              6  16   1998 s      the four slowest of the 34
//      3854  X_0^6(511)/W16            4  16   1905 s
//      9325  X_0^35(66)/W32            5  32   1674 s
//      7060  X_0^15(182)/W32           4  32   1621 s
//      7296  X_0^21(20)/W4             7   4    567 s      the three timed in 2026-09-30's run
//      1071  X_0(240)/W4               6   4    372 s
//     13029  X_0^210(73)/W32           4  32     86 s
//     12610  X_0^210(17)/W8            3   8     18 s      the four fastest of the 34
//       976  X_0(210)/W8               3   8     15 s
//      2176  X_0^6(89)/W2              3   2      4 s
//        94  X_0(30)/W1                3   1      4 s
//
// The data no longer decides anything: the proxy returns 0 for a curve already marked decided, so
// the decision is cleared on the loaded copies before ranking.

curves := eval Read("data/curves_after_UpdateCurves7.dat");
for X in curves do
    if assigned X`IsSubhyp then delete X`IsSubhyp; end if;
end for;

stage := "FilterByWeilPolynomial";
slow := [1436, 3854, 9325, 7060];       // 27 to 33 min
fast := [12610, 976, 2176, 94];         // 4 to 18 s
shapes := AssociativeArray();
shapes[1436] := <1,690,6>;  shapes[3854] := <6,511,4>;  shapes[9325] := <35,66,5>;
shapes[7060] := <15,182,4>; shapes[7296] := <21,20,7>;  shapes[1071] := <1,240,6>;
shapes[13029] := <210,73,4>; shapes[12610] := <210,17,3>; shapes[976] := <1,210,3>;
shapes[2176] := <6,89,3>;   shapes[94] := <1,30,3>;
for id -> sh in shapes do
    X := curves[id];
    assert X`CurveID eq id;
    assert <X`D, X`N, X`g> eq sh;       // the IDs still name the curves the times were measured on
end for;

proxy := AssociativeArray();
for id -> sh in shapes do proxy[id] := CurveCostProxy(curves[id], stage); end for;
printf "  proxy: %o\n", [<id, RealField(4)!proxy[id]> : id in Sort([k : k in Keys(proxy)])];

// (1) every curve measured in the tens of minutes ranks above every curve measured in seconds.
for a in slow do
    for b in fast do
        assert proxy[a] gt proxy[b];
    end for;
end for;

// (2) the three curves of the earlier run, in their measured order (567 s, 372 s, 86 s).
assert proxy[7296] gt proxy[1071];
assert proxy[1071] gt proxy[13029];

// (3) ranks among all curves of genus >= 3, printed for information: no rank has been measured,
// so none is asserted.
ranked := [c : c in curves | c`g ge 3];
pr := [<CurveCostProxy(c, stage), c`CurveID> : c in ranked];
Sort(~pr, func<a, b | b[1] - a[1]>);
ids := [t[2] : t in pr];
printf "  ranks among %o curves of genus >= 3: %o\n", #ranked,
       [<id, Position(ids, id)> : id in slow cat [7296, 1071, 13029] cat fast];

// Skipped curves still cost nothing, and a decided curve is 0.
assert CurveCostProxy(curves[1], stage) eq 0;   // genus 0
Xdec := curves[13029]; Xdec`IsSubhyp := true;
assert CurveCostProxy(Xdec, stage) eq 0;        // decided
