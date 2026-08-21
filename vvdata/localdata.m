// COMPLETE LOCAL INPUT for redoing KRY section 2's case analysis with a level structure.
// KRY sec.6 flags that their Prop 2.8 ("these are the only cases which contribute") is derived for
// sections given by char functions of the completions of O_k, and that a general section needs
// "more elaborate calculations of Whittaker functions and their derivatives".  This computes those:
// the local Whittaker POLYNOMIALS of the binary lattice L_- at every relevant prime, classified by
// the data KRY's case analysis is indexed on -- the splitting type chi_d(p), the level/ramification
// of p, and the conductor exponent -- so the case list can be re-enumerated honestly.
AttachSpec("ShimuraQuotients.spec");

MMAX := 60; if assigned MM then MMAX := StringToInteger(MM); end if;
R<X> := PolynomialRing(Rationals());

for b in [<10,3,-8>, <15,2,-7>, <6,5,-19>] do
    D := b[1]; N := b[2]; d := b[3];
    Ld := ShimuraCurveLattice(D,N);
    Q := ChangeRing(Ld`Q, Integers());
    lam := ElementOfNorm(Ld`Q, -d, Ld`O, Ld`basis_L);
    Lminus := Kernel(Transpose(Matrix(lam*Q)));
    Bm := BasisMatrix(Lminus);
    Delta := Determinant(Bm*Q*Transpose(Bm));
    df := FundamentalDiscriminant(d);
    zero := Vector(Rationals(), [0,0,0]);
    printf "\n=================== X0^%o(%o), d = %o, disc(L_-) = %o ===================\n",
           D, N, d, Factorization(Delta);
    // group the polynomials by their classifying data
    seen := AssociativeArray();
    for m in [1..MMAX] do
        Sm := Sort(Setseq({p : p in PrimeDivisors(Delta)} join {p : p in PrimeDivisors(m)}));
        for p in Sm do
            W := R ! LocalWhittakerPolynomial(Rationals()!m, p, zero, Lminus, Q);
            chi := (df mod p eq 0) select 0 else KroneckerSymbol(df, p);
            role := (N mod p eq 0) select "LEVEL" else ((D mod p eq 0) select "ram(D)" else "good");
            // conductor exponent as KRY/Yang index it: 4m = d_m c^2 for the field of sqrt(m)
            sq := SquarefreeFactorization(m);
            dm := (sq eq 1) select 1 else Discriminant(Integers(QuadraticField(sq)));
            isq, c := IsSquare(Rationals()!(4*m)/dm);
            kp := isq select Valuation(Integers()!c, p) else -1;
            key := <role, chi, Valuation(m,p), kp>;
            if not IsDefined(seen, key) then seen[key] := {@ @}; end if;
            Include(~seen[key], W);
        end for;
    end for;
    ks := Sort(Setseq(Keys(seen)));
    printf "  %-8o %-6o %-8o %-4o %o\n", "role", "chi_d", "ord_p(m)", "k_p", "W_p(X)  (distinct polys in this class)";
    for k in ks do
        ps := [f : f in seen[k]];
        printf "  %-8o %-6o %-8o %-4o %o%o\n", k[1], k[2], k[3], k[4],
               #ps le 2 select ps else ps[1..2],
               #ps gt 2 select Sprintf(" ...(%o distinct)", #ps) else "";
    end for;
end for;
quit;
