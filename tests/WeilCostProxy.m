// CurveCostProxy for the Weil-polynomial stage must put the measured-heavy curves first.  The
// estimate is #W * sum p^(g/2) over the good primes the stage will use, the shape of the
// dominant trace-formula call (see Utils.m).
//
// Every assertion below is an ORDERING of two curves whose times were measured, never a cost.
// Measured on lovelace with the class-number tables, over 36 curves, on main at 95b19e6 (with
// #56, #57 and #58), by vvdata/weyl-campaign/weil-retime-2026-10-02/weil_timing.m on the
// m0-theta-campaign branch (logs and the fit alongside it).  The first column is the curve's
// index in the 18,379-quotient list of that run; the test finds each curve by D, N and W, since
// a change to the enumeration renumbers the list.
//
//     curve  shape                     g  #W   measured
//      1416  X_0(595)/W8               5   8   6087 s      the four slowest of the 36
//      2325  X_0^6(97)/W2              6   2   4412 s
//      1436  X_0(690)/W16              6  16   1998 s
//      3854  X_0^6(511)/W16            4  16   1905 s
//      7296  X_0^21(20)/W4             7   4    567 s      the three timed in 2026-09-30's run
//      1071  X_0(240)/W4               6   4    372 s
//     13029  X_0^210(73)/W32           4  32     86 s
//     12610  X_0^210(17)/W8            3   8     18 s      the four fastest of the 36
//       976  X_0(210)/W8               3   8     15 s
//      2176  X_0^6(89)/W2              3   2      4 s
//        94  X_0(30)/W1                3   1      4 s
//
// ⚠ The estimate orders 246 of the 271 pairs differing by more than a factor 10, so it does
// mis-order real pairs: 2325 (74 min) sits BELOW 13029 (86 s).  Nothing below asserts a pair the
// estimate gets wrong, and no claim here is a cost.
//
// The data no longer decides anything: the proxy returns 0 for a curve already marked decided, so
// the decision is cleared on the loaded copies before ranking.

curves := eval Read("data/curves_after_UpdateCurves7.dat");
for X in curves do
    if assigned X`IsSubhyp then delete X`IsSubhyp; end if;
end for;

stage := "FilterByWeilPolynomial";
slow := [1416, 2325, 1436, 3854];       // 32 min to 101 min
fast := [12610, 976, 2176, 94];         // 4 to 18 s
timed := AssociativeArray();            // <D, N, W, g> of each timed curve, keyed as in the table
timed[1416]  := <1, 595, {1, 5, 7, 17, 35, 85, 119, 595}, 5>;
timed[2325]  := <6, 97, {1, 291}, 6>;
timed[1436]  := <1, 690, Set(Divisors(690)), 6>;
timed[3854]  := <6, 511, Set(Divisors(3066)), 4>;
timed[7296]  := <21, 20, {1, 12, 15, 20}, 7>;
timed[1071]  := <1, 240, {1, 15, 48, 80}, 6>;
timed[13029] := <210, 73, Set(Divisors(15330)), 4>;
timed[12610] := <210, 17, {1, 15, 42, 51, 70, 85, 238, 3570}, 3>;
timed[976]   := <1, 210, {1, 3, 5, 7, 15, 21, 35, 105}, 3>;
timed[2176]  := <6, 89, {1, 534}, 3>;
timed[94]    := <1, 30, {1}, 3>;
curve := AssociativeArray();
for id -> t in timed do
    found := [X : X in curves | X`D eq t[1] and X`N eq t[2] and X`W eq t[3]];
    assert #found eq 1;
    assert found[1]`g eq t[4];
    curve[id] := found[1];
end for;

proxy := AssociativeArray();
for id -> X in curve do proxy[id] := CurveCostProxy(X, stage); end for;
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
       [<id, Position(ids, curve[id]`CurveID)> : id in slow cat [7296, 1071, 13029] cat fast];

// Skipped curves still cost nothing, and a decided curve is 0.
assert CurveCostProxy([X : X in curves | X`g eq 0][1], stage) eq 0;   // genus 0
Xdec := curve[13029]; Xdec`IsSubhyp := true;
assert CurveCostProxy(Xdec, stage) eq 0;        // decided
