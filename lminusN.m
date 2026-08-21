// Theorem 4.1's logs come from the places where the BINARY lattice L_- fails to represent m.
// For a maximal order L_- is unimodular at a split prime, so a split prime never obstructs -- hence
// no split-prime log.  But at a FIRING disc the pipeline's own L_- has discriminant Delta = r^2 |d|
// with N | r, so L_- is N-SCALED, not unimodular, at N.  That is a configuration the maximal-order
// theory never produces.  Does its local Whittaker at N vanish -- i.e. is there a log N channel?
AttachSpec("ShimuraQuotients.spec");

for b in [<15,2,[-7,-15],[2,10,30,1,3,5]>, <6,5,[-19,-24],[10,15,30,1,2,3]>,
          <10,3,[-8,-20],[3,12,30,1,2,5]>] do
    D := b[1]; N := b[2]; discs := b[3]; ms := b[4];
    Ld := ShimuraCurveLattice(D,N);
    Q := ChangeRing(Ld`Q, Integers());
    printf "\n=========== X0^%o(%o), N = %o ===========\n", D, N, N;
    for d in discs do
        ok := true; lam := 0;
        try lam := ElementOfNorm(Ld`Q, -d, Ld`O, Ld`basis_L); catch e ok := false; end try;
        if not ok then printf "  d = %o: no lambda\n", d; continue; end if;
        Lminus := Kernel(Transpose(Matrix(lam*Q)));
        Bm := BasisMatrix(Lminus);
        Delta := Determinant(Bm*Q*Transpose(Bm));
        df := FundamentalDiscriminant(d);
        printf "  --- d = %-5o (d_fund %o, N is %o) : Delta(L_-) = %o = %o, ord_%o = %o\n",
               d, df, (df mod N eq 0) select "ramified" else
                      (KroneckerSymbol(df,N) eq 1 select "SPLIT" else "inert"),
               Delta, Factorization(Delta), N, Valuation(Delta, N);
        zero := Vector(Rationals(), [0,0,0]);
        for m in ms do
            Sm := Sort(Setseq({p : p in PrimeDivisors(Delta)} join {p : p in PrimeDivisors(m)}));
            vals := [];
            for p in Sm do
                w := LocalWhittakerAtOne(Rationals()!m, p, zero, Lminus, Q);
                Append(~vals, <p, w>);
            end for;
            vanish := [v[1] : v in vals | v[2] eq 0];
            printf "      m = %-4o %-8o W(L_-) = %-34o vanishing at %o\n",
                   m, (m mod N eq 0) select "(N | m)" else "", vals, vanish;
        end for;
    end for;
end for;
quit;
