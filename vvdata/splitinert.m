// Does the needed log N appear at firing discs where N SPLITS in Q(sqrt d)?
// Schofer Thm 4.1 produces logs only at primes RAMIFIED or INERT in k; a split prime contributes
// none.  The pipeline's correction is mult * kzero_N(d) with kzero_N = log N at EVERY firing disc
// (N unramified), split or inert alike.  If both types occur among the discs the models are
// validated on, then no formula of Thm 4.1's shape can supply the correction -- the gap is in the
// theory, not in the identification of the lattice.
AttachSpec("ShimuraQuotients.spec");
for b in [[15,2],[6,5],[10,3],[21,2]] do
    D := b[1]; N := b[2];
    Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
    Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};
    rat, quad := RationalandQuadraticCMPoints(star : coprime_to_level := false, bd := 4);
    ds := Sort(Setseq({a[1] : a in rat} join {a[1] : a in quad}));
    printf "\n=== X0^%o(%o), N = %o : %o CM discs ===\n", D, N, N, #ds;
    tally := AssociativeArray();
    for k in ["ramified","split (FIRING)","inert (FIRING)"] do tally[k] := []; end for;
    for d in ds do
        df := FundamentalDiscriminant(d);
        if df mod N eq 0 then
            Append(~tally["ramified"], d);
        else
            ks := KroneckerSymbol(df, N);
            if ks eq 1 then Append(~tally["split (FIRING)"], d);
            else Append(~tally["inert (FIRING)"], d); end if;
        end if;
    end for;
    for k in ["ramified","split (FIRING)","inert (FIRING)"] do
        printf "  %-16o %-3o : %o\n", k, #tally[k],
               #tally[k] le 12 select tally[k] else tally[k][1..12];
    end for;
end for;
quit;
