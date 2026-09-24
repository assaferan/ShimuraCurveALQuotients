// THE TEST tests/Kappa0.m never does: Kappa0(m, d) at m DIVISIBLE BY N, at a FIRING disc.
//
// The oracle pins the answer.  The correction the pipeline needs is mult * kzero_N(d), with
// mult = sum_m c(-m) A_m and A_m = -b(m)/4 solved exactly from the oracle; kzero_N(d) = log N at
// every firing disc.  So the per-index missing amount is A_m * log N, INDEPENDENT of d.  Since the
// code's log-N coefficient at firing discs is (reportedly) zero, the prediction is
//        true coefficient of log N in Kappa0(m,d)  =  A_m   for every firing d.
// Here we print what the code actually returns, in full, so the structure of the gap is visible.
AttachSpec("ShimuraQuotients.spec");

bases := [
  <15, 2, [-7, -15, -23], [<2,-1>, <10,1>, <30,0>], [1,3,5,15]>,
  < 6, 5, [-8, -19, -24], [<10,3/2>, <15,1/2>, <30,3/2>], [1,2,3,6]>,
  <10, 3, [-8, -20, -35], [<3,1/2>, <12,1/2>, <30,3/2>], [1,2,5,10]>
];

for b in bases do
    D := b[1]; N := b[2]; discs := b[3]; targets := b[4]; ctrl := b[5];
    Ld := ShimuraCurveLattice(D, N);
    Q := Ld`Q;
    printf "\n================ X0^%o(%o)   N = %o ================\n", D, N, N;
    for d in discs do
        df := FundamentalDiscriminant(d);
        firing := (df mod N) ne 0;
        ok := true; lam := 0;
        try
            lam := ElementOfNorm(Q, -d, Ld`O, Ld`basis_L);
        catch e
            ok := false;
        end try;
        if not ok then printf "  d = %-5o : no lambda of norm %o\n", d, -d; continue; end if;
        printf "  --- d = %-5o (d_fund = %-5o, %o) ---\n", d, df,
               firing select "FIRING" else "non-firing";
        for t in targets cat [<m, 0> : m in ctrl] do
            m := t[1]; want := t[2];
            k := Kappa0(m, d, Q, lam);
            co := k`log_coeffs;
            ks := Sort([p : p in Keys(co)]);
            atN := IsDefined(co, N) select co[N] else 0;
            printf "      m = %-4o %-8o code log_%o coeff = %-8o | target A_m = %-6o | all: %o\n",
                   m, (m mod N eq 0) select "(N | m)" else "", N, atN, want,
                   [<p, co[p]> : p in ks];
        end for;
    end for;
end for;
quit;
