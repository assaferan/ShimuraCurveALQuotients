// Parity survey for the even-correction escape hatch.  THROWAWAY branch only.
//
//   magma -b DD:=38 NN:=5 annprobe.m
//
// The escape hatch (see memory borcherds-obstruction-is-real) turns on ONE number per
// obstructed base: phi(target), where phi generates the annihilator of the divisor map's
// image.  If phi(target) is EVEN, an even divisor correction can send it to 0 without
// changing the cover mod 2.  If it is ODD, that base is out of reach of this method.
//
// BorcherdsForms is patched to print PROBEANN* at the first failing key and then abort
// with "PROBE_STOP", so each base costs one key instead of an exhaustive triple sweep.
AttachSpec("ShimuraQuotients.spec");
D := 38; N := 5;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
printf "ANNPROBE BASE %o_%o\n", D, N;
t0 := Realtime();
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
try
    fs := BorcherdsForms(star, curves : Prec := 100);
    printf "ANNPROBE RESULT %o_%o SURJECTIVE nforms %o time %os\n", D, N, #Keys(fs), Realtime() - t0;
catch e
    obj := Sprint(e`Object);
    if "PROBE_STOP" in obj then
        printf "ANNPROBE RESULT %o_%o OBSTRUCTED time %os\n", D, N, Realtime() - t0;
    else
        printf "ANNPROBE RESULT %o_%o OTHER time %os\n", D, N, Realtime() - t0;
        printf "ANNPROBE_ERR %o\n", obj;
    end if;
end try;
printf "ANNPROBE DONE %o_%o\n", D, N;
quit;
