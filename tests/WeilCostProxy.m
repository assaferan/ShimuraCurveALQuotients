// CurveCostProxy for the Weil stage must put the measured-heavy shape first.  Measured
// 2026-09-30 on main at 95cf87b, before #56 and #58, on a Mac without the class-number tables,
// with vvdata/weyl-campaign/weil-cost-2026-09-30/weil_timing.m (which replays WeilPolynomial(X, p)
// over the stage's own prime bound; logs alongside it): curve 13029 = X_0^210(73)/W32 (g = 4,
// Qmax = 15330) took 559 s for all seven of its primes; curve 1071 = X_0(240)/W4 (g = 6,
// Qmax = 80) took 2959 s through p = 17 alone; curve 7296 = X_0^21(20)/W4 (g = 7) is heavier
// still (p = 13 alone: 1137 s with #56).  On lovelace with the class-number tables, to the
// full prime bound (same directory, lovelace/): 13029 took 40 min and 1071 14.8 h on main,
// 87 min with #58; the order is the same under every tree.  The old proxy (sum 4*Qmax*p^g)
// put 13029 above 1071.
//
// The assertions are order claims only, never a cost, and they do not depend on which curves
// the data currently marks as decided: the proxy returns 0 for a decided curve, so the
// decision is cleared on the loaded copies before ranking.

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
assert p1071 gt p13029;           // the 9-h curve above the 9-min one
assert p7296 gt p1071;            // and genus 7 above genus 6

// Ranking every curve of genus >= 3 (the stage skips lower genus), decided or not, so the
// claim does not move with the data.  The proxy sums n = p^g over the admissible primes, so it
// orders by genus first and by the number of admissible primes second:
//   * the two curves measured in hours (genus 6 and 7) are in the top 1100 of 14231, above
//     every genus <= 5 curve and most genus-6 ones with a larger W;
//   * the curve measured in minutes (genus 4) is at least five times further down the list;
//   * the top 1000 are all genus >= 6 -- the shape the hours come from, and no genus-3 or 4
//     curve was ever measured above a few minutes.
ranked := [c : c in curves | c`g ge 3];
pr := [<CurveCostProxy(c, stage), c`CurveID, c`g> : c in ranked];
Sort(~pr, func<a, b | b[1] - a[1]>);
ids := [t[2] : t in pr];
r1071 := Position(ids, 1071); r7296 := Position(ids, 7296); r13029 := Position(ids, 13029);
printf "  ranks among %o curves of genus >= 3: 7296 -> %o, 1071 -> %o, 13029 -> %o\n", #ranked, r7296, r1071, r13029;
assert r1071 le 1100 and r7296 le 1100;
assert r13029 ge 5*r1071;
assert &and[t[3] ge 6 : t in pr[1..1000]];

// Skipped curves still cost nothing, and a decided curve is 0.
assert CurveCostProxy(curves[1], stage) eq 0;   // genus 0
Xdec := curves[13029]; Xdec`IsSubhyp := true;
assert CurveCostProxy(Xdec, stage) eq 0;        // decided
