// The parity census of the Borcherds obstruction: run the Borcherds search of a base with the
// OBSTPAIR instrumentation of branch obstruction-pairing (BorcherdsForms.m, after `found_v`), which
// prints, for every (key, anchor triple) whose target divisor is NOT in the image of the divisor map,
// the pairing of that target with each primitive generator of the annihilator of the image.  The
// search still fails at the end ("Failed to find all Borcherds forms"); only the printed lines matter.
//
//   cd worktrees/obstr && OBSTPAIR=1 magma -b DD:=38 NN:=5 ../campaign/vvdata/weyl-campaign/obstruction/obstpair.m
//
// Control: X_0^38(5), key 11, first triple: pairing -22 (memory borcherds-obstruction-is-real).
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
D := StringToInteger(DD); N := StringToInteger(NN);
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
printf "OBSTPAIR BASE %o %o\n", D, N;
t0 := Realtime();
try
    fs := BorcherdsForms(star, curves : Prec := 100);
    printf "OBSTPAIR RESULT %o %o: forms found for keys %o (%o s)\n", D, N, Sort([k : k in Keys(fs)]), Realtime(t0);
catch e
    printf "OBSTPAIR RESULT %o %o: %o (%o s)\n", D, N, e`Object, Realtime(t0);
end try;
quit;
