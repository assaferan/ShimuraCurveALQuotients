// Diagnostic: attempt the X0^15(2) model build on current main (post PR #19).
// The last recorded blocker was the find_signs_hauptmodul crash, whose deeper root was that discs
// -7, -15, -60 violate s + s~ = 1. Those are exactly the 2-split discs the newly merged m=0 term
// corrects, so that root may already be gone. This just runs it and reports what happens.

AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 1);
fh := Open("build_out.txt", "w");
emit := procedure(s) fprintf fh, "%o\n", s; Flush(fh); end procedure;

D := 15; N := 2;
emit(Sprintf("=== building X0^%o(%o) ===", D, N));
t0 := Cputime();
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
// use the star curve FROM the returned list: GetQuotientsAndGenera is what assigns `Covers`
_ := exists(star){c : c in curves | IsStarCurve(c)};
emit(Sprintf("curves: %o  (star genus %o)", #curves, Xstar`g));
try
    eqns, covs := AllEquationsAboveCovers(star, curves : Prec := 100);
    emit(Sprintf("SUCCESS: #equation entries = %o", #Keys(eqns)));
    for k in Sort([kk : kk in Keys(eqns)]) do
        emit(Sprintf("  [%o] -> %o", k, eqns[k]));
    end for;
catch e
    emit("CRASH:");
    emit(Sprintf("%o", e`Object));
end try;
emit(Sprintf("elapsed %o s", Cputime(t0)));
emit("DONE");
