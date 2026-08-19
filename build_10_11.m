// ODD-N m=0 ORACLE.
//
// The premise "the true m=0 multiplier is 0 for every odd-N form" has never actually been tested; it
// rests only on "main passes the odd-N tests with no m=0 term at all". Here is why that is weak
// evidence: on X0^10(11) every CM disc currently used is coprime to 11, so kzero_N = log 11 would fire
// at EVERY disc, and a term that is uniform across a row is just an overall row scaling -- which
// y^2 = f(s) absorbs. The omission is therefore structurally INVISIBLE on the current point set.
//
// Admitting discs with 11 | d (where the term does NOT fire) breaks that uniformity and makes the
// multiplier IDENTIFIABLE -- and indeed that is exactly when the build starts failing with
// "QuadraticConstraintsOnEquations: no solution found".
//
// So: build with the non-coprime rational CM points admitted, and dump per-cover (ds, s, y2). The fit
// (done outside) solves, per cover, for the factor 11^delta applied at the FIRING discs only. delta = 0
// for every cover => the multiplier really is 0 for odd N. A consistent nonzero delta => we have just
// MEASURED the odd-N m=0 multiplier.

AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 1);
fh := Open("build_10_11_out.txt", "w");
emit := procedure(s) fprintf fh, "%o\n", s; Flush(fh); end procedure;

D := 10; N := 11;
emit(Sprintf("=== building X0^%o(%o) ===", D, N));
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
emit(Sprintf("curves: %o  (star genus %o)", #curves, Xstar`g));
t0 := Cputime();
try
    eqns, covs := AllEquationsAboveCovers(star, curves : Prec := 100);
    emit(Sprintf("SUCCESS: #equation entries = %o", #Keys(eqns)));
catch e
    emit("CRASH:");
    emit(Sprintf("%o", e`Object));
end try;
emit(Sprintf("elapsed %o s", Cputime(t0)));
emit("DONE");
exit;
