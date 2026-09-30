// One deep Weil-stage call, no profiler: curve 1071 (X_0(240)/<w15,w48,w80>, g = 6) at n = 13^6.
AttachSpec("ShimuraQuotients.spec");
W := {Integers() | 1, 15, 48, 80};
t0 := Cputime();
tr := TraceDNewALFixed(1, 240, 2, 13^6, W);
printf "DEEP %o: TraceDNewALFixed(1,240,2,13^6,W) = %o in %o s\n", tag, tr, RealField(5)!Cputime(t0);
exit;
