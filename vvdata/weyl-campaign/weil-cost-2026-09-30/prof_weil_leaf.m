// Where does a single deep Weil-stage call spend its time?  Curve 1071 at p = 13, v = 6
// (119 s on the memo code).  Profiler by total time, top 30.
AttachSpec("ShimuraQuotients.spec");
curves := eval Read("data/curves_after_UpdateCurves7.dat");
X := curves[1071];
printf "X: D=%o N=%o W=%o g=%o\n", X`D, X`N, X`W, X`g;
SetProfile(true);
t0 := Cputime();
tr := TraceDNewALFixed(X`D, X`N, 2, 13^6, X`W);
printf "TraceDNewALFixed(.., 13^6) = %o in %o s\n", tr, Cputime(t0);
SetProfile(false);
G := ProfileGraph();
ProfilePrintByTotalTime(G : Max := 30);
exit;
