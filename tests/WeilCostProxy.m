// CurveCostProxy for the Weil stage must put the measured-heavy shape first.  Measured on a Mac
// without the class-number tables (2026-09-30, main at 95cf87b), replaying WeilPolynomial(X, p)
// over the stage's own prime bound: curve 13029 = X_0^210(73)/W32 (g = 4, Qmax = 15330) took
// 559 s for all seven primes; curve 1071 = X_0(240)/W4 (g = 6, Qmax = 80) took 2959 s through
// p = 17 alone (~9 h projected); curve 7296 = X_0^21(20)/W4 (g = 7) is heavier still.  The old
// proxy (sum 4*Qmax*p^g) weighted 13029 at 1.3e11 above 1071 at 7.2e10 and ranked 1071 138th of
// 886 open curves, so heavy-first dispatch put the heavy curve late.  The cost is set by
// n = p^g, not by Qmax, and the proxy now sums p^g.
//
// The assertions are ORDER claims only (measured), never a cost.  On the old code the first
// assert fails: proxy(13029) > proxy(1071).

curves := eval Read("data/curves_after_UpdateCurves7.dat");
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

// Among all open curves the two measured-heavy ones are in the top 10, and the top 50 are all
// genus >= 5: the shape the stage's hours come from.
open := [c : c in curves | not assigned c`IsSubhyp and c`g ge 3];
pr := [<CurveCostProxy(c, stage), c`CurveID, c`g> : c in open];
Sort(~pr, func<a, b | b[1] - a[1]>);
ids := [t[2] : t in pr];
r1071 := Position(ids, 1071); r7296 := Position(ids, 7296); r13029 := Position(ids, 13029);
printf "  ranks of %o open curves: 7296 -> %o, 1071 -> %o, 13029 -> %o\n", #open, r7296, r1071, r13029;
assert r1071 le 10 and r7296 le 10;
assert r13029 gt 100;
assert &and[t[3] ge 5 : t in pr[1..50]];

// Skipped curves still cost nothing, and a decided curve is 0.
assert CurveCostProxy(curves[1], stage) eq 0;   // genus 0, decided
