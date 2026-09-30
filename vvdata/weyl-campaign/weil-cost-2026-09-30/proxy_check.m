// Does CurveCostProxy rank the Weil stage's heavy curves right?  Measured on this Mac (main):
// 13029 (210,73) g=4 #W=32: 559 s for all 7 primes; 1071 X_0(240)/W4 g=6: 2959 s through
// p=17 alone, ~9 h projected.  Print the proxy for both, and the top 10 open curves by proxy.
AttachSpec("ShimuraQuotients.spec");
curves := eval Read("data/curves_after_UpdateCurves7.dat");
for id in [13029, 1071] do
    X := curves[id];
    printf "id %o: D=%o N=%o g=%o #W=%o Qmax=%o  proxy=%o\n", id, X`D, X`N, X`g, #X`W, Maximum(X`W),
        CurveCostProxy(X, "FilterByWeilPolynomial");
end for;
open := [c : c in curves | not assigned c`IsSubhyp and c`g ge 3];
pr := [<CurveCostProxy(c, "FilterByWeilPolynomial"), c`CurveID, c`D, c`N, c`g, #c`W> : c in open];
Sort(~pr, func<a, b | b[1] - a[1]>);
printf "top 12 open curves by Weil proxy (proxy, id, D, N, g, #W):\n";
for t in pr[1..12] do printf "  %o\n", t; end for;
printf "rank of 1071 among %o open curves: %o\n", #open, Position([t[2] : t in pr], 1071);
exit;
