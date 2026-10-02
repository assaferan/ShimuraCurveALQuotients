// Genera of the quotients that matter for the obstructed non-squarefree bases: by the Hall-divisor
// Atkin-Lehner group W (the "star" curve), and by W_D alone (the quotient that is the Galois closure
// of X^* over X_0^D(N/h^2)^*).  Also which intermediate levels have models.
AttachSpec("ShimuraQuotients.spec");
hall := func<M | {d : d in Divisors(M) | GCD(d, M div d) eq 1}>;
printf "%-6o %-6o %-10o %-10o %-14o %o\n", "base", "g(X)", "g(X/W)", "g(X/W_D)", "base star", "models for intermediate levels";
for b in [[15,4],[21,4],[33,4],[15,8],[10,9],[14,9],[22,9]] do
    D := b[1]; N := b[2];
    W := hall(D*N); WD := hall(D);
    g := GenusShimuraCurveQuotient(D, N, {1});
    gW := GenusShimuraCurveQuotient(D, N, W);
    gWD := GenusShimuraCurveQuotient(D, N, WD);
    h := (N mod 9 eq 0) select 3 else 2;
    Np := N div h^2;
    gbs := GenusShimuraCurveQuotient(D, Np, hall(D*Np));
    inter := [Sprintf("%o_%o:%o", D, M, Pipe(Sprintf("test -f data/models/models_%o_%o.m && echo yes || echo no", D, M), "")[1..2])
              : M in [Np*h, Np*h^2] | M ne N] cat [Sprintf("%o_%o:%o", D, Np, Pipe(Sprintf("test -f data/models/models_%o_%o.m && echo yes || echo no", D, Np), "")[1..2])];
    printf "%-6o %-6o %-10o %-10o %-14o %o\n", Sprintf("%o_%o", D, N), g, gW, gWD, Sprintf("%o_%o*: g=%o", D, Np, gbs), inter;
end for;
