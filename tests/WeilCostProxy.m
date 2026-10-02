// CurveCostProxy for the Weil stage must put the measured-heavy shape first.  Measured on
// lovelace with the class-number tables, to the stage's full prime bound, 2026-09-30, with
// vvdata/weyl-campaign/weil-cost-2026-09-30/weil_timing.m on the m0-theta-campaign branch (it
// replays WeilPolynomial(X, p); logs alongside it, lovelace/): curve 1071 = X_0(240)/W4 (g = 6,
// Qmax = 80) took 14.8 h on main at 95cf87b and 87 min with #58; curve 13029 = X_0^210(73)/W32
// (g = 4, Qmax = 15330) took 40 min on main; curve 7296 = X_0^21(20)/W4 (g = 7) is heavier
// still (p = 13 alone: 1137 s with #56).  The order is the same under every tree.
//
// The assertions are order claims on those three measured curves only, never a cost, and
// they do not depend on which curves the data currently marks as decided: the proxy returns 0
// for a decided curve, so the decision is cleared on the loaded copies before ranking.

curves := eval Read("data/curves_after_UpdateCurves7.dat");
for X in curves do
    if assigned X`IsSubhyp then delete X`IsSubhyp; end if;
end for;
X13029 := curves[13029]; X1071 := curves[1071]; X7296 := curves[7296];
assert X13029`CurveID eq 13029 and X1071`CurveID eq 1071 and X7296`CurveID eq 7296;
assert <X13029`D, X13029`N, X13029`g> eq <210, 73, 4>;
assert <X1071`D, X1071`N, X1071`g> eq <1, 240, 6>;
assert <X7296`D, X7296`N, X7296`g> eq <21, 20, 7>;

stage := "FilterByWeilPolynomial";
p13029 := CurveCostProxy(X13029, stage);
p1071  := CurveCostProxy(X1071, stage);
p7296  := CurveCostProxy(X7296, stage);
printf "  proxy: 13029 (g=4, Qmax=15330) %o | 1071 (g=6, Qmax=80) %o | 7296 (g=7) %o\n", p13029, p1071, p7296;
assert p1071 gt p13029;           // 87 min (14.8 h on main) above 40 min
assert p7296 gt p1071;            // and the genus-7 curve, heavier per prime, above both

// The ranks among every curve of genus >= 3 (the stage skips lower genus), decided or not, are
// printed for information; no rank is asserted, because no rank has been measured.
ranked := [c : c in curves | c`g ge 3];
pr := [<CurveCostProxy(c, stage), c`CurveID, c`g> : c in ranked];
Sort(~pr, func<a, b | b[1] - a[1]>);
ids := [t[2] : t in pr];
r1071 := Position(ids, 1071); r7296 := Position(ids, 7296); r13029 := Position(ids, 13029);
printf "  ranks among %o curves of genus >= 3: 7296 -> %o, 1071 -> %o, 13029 -> %o\n", #ranked, r7296, r1071, r13029;

// Skipped curves still cost nothing, and a decided curve is 0.
assert CurveCostProxy(curves[1], stage) eq 0;   // genus 0
Xdec := curves[13029]; Xdec`IsSubhyp := true;
assert CurveCostProxy(Xdec, stage) eq 0;        // decided
