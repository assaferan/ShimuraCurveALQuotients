// Per-stage profiling driver for the odd-D branch of BorcherdsForms.  THROWAWAY.
//   BFPROF=1 magma DD:=65 NN:=2 bfprof.m
// Mirrors ppint.m's setup exactly, but wraps BorcherdsForms in try/catch so a runtime
// error still prints the SUMMARY line and exits (magma -b hangs on an uncaught error).
AttachSpec("ShimuraQuotients.spec");
D := 65; N := 2;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
t0 := Realtime();
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
printf "BFPROF setup createquot %os\n", Realtime() - t0;
_pt := Realtime();
curves := GetQuotientsAndGenera([Xstar]);
printf "BFPROF setup getquotients %os ncurves %o\n", Realtime() - _pt, #curves;
_ := exists(star){c : c in curves | IsStarCurve(c)};
_pt := Realtime();
try
    fs := BorcherdsForms(star, curves : Prec := 100);
    printf "BFPROF RESULT ok nforms %o bf %os\n", #Keys(fs), Realtime() - _pt;
catch e
    printf "BFPROF RESULT error bf %os\n", Realtime() - _pt;
    print e`Object;
end try;
printf "BFPROF setup total %os\n", Realtime() - t0;
printf "DONE\n";
quit;
