// Dump the log-N content of Kappa0(m,d) for m divisible by N vs not, at firing and non-firing discs.
AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 0);
D := StringToInteger(DD); N := StringToInteger(NN);
Ld := ShimuraCurveLattice(D, N);
Q := Ld`Q;
discs := [StringToInteger(x) : x in Split(DISCS, ",")];
printf "[K] D=%o N=%o\n", D, N;
printf "[K] %-6o %-6o %-8o %-9o %o\n", "d", "fires", "m", "N|m", "coeff of log N in Kappa0(m,d)";
for d in discs do
    fires := not IsEmpty(PrimeDivisors(N div GCD(N, FundamentalDiscriminant(d))));
    lam := ElementOfNorm(Q, -d, Ld`O, Ld`basis_L);
    for m in [1..12] do
        ok := true;
        try k := Kappa0(m, d, Q, lam); catch e ok := false; end try;
        if not ok then continue; end if;
        // extract the coefficient of log N from the LogSum
        c := 0;
        lc := k`log_coeffs;
        if IsDefined(lc, N) then c := lc[N]; end if;
        printf "[K] %-6o %-6o %-8o %-9o %o\n", d, fires, m, (m mod N eq 0), c;
    end for;
end for;
exit;
