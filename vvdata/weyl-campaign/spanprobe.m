// PROBE_SPAN driver -- throwaway worktree only.
//
//   PROBE_BUMP=0 magma -b DD:=38 NN:=5 spanprobe.m
//
// Question under test: "Failed to find all Borcherds forms" fires at
// BorcherdsForms.m's `found_v` -- the divisor target is not in the image of coeffs_trunc.
// At 38_5 that was measured as a RANK-1 deficiency (rank 35 of 36 columns) at a single key,
// identical across all 96 (infinity, pair) triples.  If the weakly holomorphic space is
// simply one pole order too small, deepening it (PROBE_BUMP >= 1) supplies the missing
// direction and the base becomes generable.  If the deficit survives every bump, the
// obstruction is genuine.
//
// Each key emits one PROBESPAN line; the run ends in SPANPROBE RESULT ... SUCCESS|FAIL.
AttachSpec("ShimuraQuotients.spec");
D := 38; N := 5;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
printf "SPANPROBE BASE %o %o bumpenv [%o]\n", D, N, GetEnv("PROBE_BUMP");
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
try
    fs := BorcherdsForms(star, curves : Prec := 100);
    printf "SPANPROBE RESULT %o %o SUCCESS nforms %o\n", D, N, #Keys(fs);
catch e
    printf "SPANPROBE RESULT %o %o FAIL %o\n", D, N, e`Object;
end try;
printf "DONE\n";
quit;
