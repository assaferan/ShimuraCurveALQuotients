// Measure the class-number-table speedup on the Weil-polynomial stage (the heavy
// class-number consumer: it requests Hurwitz class numbers of |d| <= 4*Qmax*p^g, all
// within the 2^28 table range by construction of the Weil prime bound).
//
// Run twice, identical except for the CLASS_GROUPS_DIR env var (read by ClassNumberData.m):
//   without tables:  magma -b measure_classnum.m
//   with tables:     CLASS_GROUPS_DIR=/scratch/assaferan/class-groups-quadratic-imaginary-fields magma -b measure_classnum.m
// The per-process class-number cache persists across curves (as in a real worker), so this
// mirrors the stage's actual behavior.  K (sample size) is taken from env K (default 12).
SetQuitOnError(true);
AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 0);

dir := GetEnv("CLASS_GROUPS_DIR");
if dir ne "" then
    printf "=== MODE: TABLES (%o) ===\n", dir;
else
    printf "=== MODE: live fallback (no tables) ===\n";
end if;

Kenv := GetEnv("K");
K := (Kenv ne "") select StringToInteger(Kenv) else 12;

curves := eval Read("data/par/curves_after_UpdateCurves6.dat");
// Optional genus band (env GMIN/GMAX): target curves where the Weil prime bound is large
// enough that the stage actually requests class numbers (lower genus / smaller Qmax).
gmin := (GetEnv("GMIN") ne "") select StringToInteger(GetEnv("GMIN")) else 3;
gmax := (GetEnv("GMAX") ne "") select StringToInteger(GetEnv("GMAX")) else 1000;
// Undetermined curves only (the ones the Weil stage actually works on).
cand := [X : X in curves | not assigned X`IsSubhyp and X`g ge gmin and X`g le gmax];
// Deterministic order: genus desc, then D*N desc, then a stable key, so both runs use the
// SAME curves in the SAME order.
Sort(~cand, func<a,b |
    (b`g ne a`g) select (b`g - a`g)
    else ((b`D*b`N ne a`D*a`N) select (b`D*b`N - a`D*a`N)
    else (Max(b`W) - Max(a`W)))>);
sample := cand[1..Min(K, #cand)];
printf "sample: %o curves (of %o undetermined g>=3)\n", #sample, #cand;

total := 0.0;
for i in [1..#sample] do
    X := sample[i];
    s := [X];
    t := Realtime();
    FilterByWeilPolynomialGenusScaled(~s);
    dt := Realtime() - t;
    total +:= dt;
    printf "[%2o] D=%o N=%o g=%o Qmax=%o  pruned=%o  %.3os  (cum %.3os)\n",
        i, X`D, X`N, X`g, Max(X`W), assigned s[1]`IsSubhyp and not s[1]`IsSubhyp, dt, total;
end for;

printf "=== TOTAL over %o curves: %.3o s ===\n", #sample, total;
quit;
