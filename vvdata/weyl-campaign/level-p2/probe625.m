AttachSpec("ShimuraQuotients.spec");
hall := func<M | {d : d in Divisors(M) | GCD(d, M div d) eq 1}>;
for b in [[6,25],[6,49]] do
    D := b[1]; N := b[2];
    printf "\n=== %o_%o: Hall(D*N) = %o\n", D, N, Sort(SetToSequence(hall(D*N)));
    t := Realtime();
    try
        Ld := ShimuraCurveLattice(D, N);
        printf "  ShimuraCurveLattice: OK in %o s, |discriminant group| = %o, denom %o\n",
               RealField(4)!(Realtime()-t), #Ld`disc_grp, Ld`denom;
    catch e
        printf "  ShimuraCurveLattice FAILED: %o\n", e`Object;
    end try;
    t := Realtime();
    try
        g := GenusShimuraCurveQuotient(D, N, hall(D*N));
        g1 := GenusShimuraCurveQuotient(D, N, {1});
        printf "  genus of the star curve %o, of the full curve %o (%o s)\n", g, g1, RealField(4)!(Realtime()-t);
    catch e
        printf "  genus FAILED: %o\n", e`Object;
    end try;
end for;
