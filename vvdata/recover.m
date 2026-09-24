// THE EXPERIMENT.  If W_N's zero is a zeta artifact rather than a Diff place, then cancelling it
// (W~_N = W_N/(1-X) = 1 when N does not divide mm) removes N from the vanishing count, and terms the
// code currently discards under its "two or more vanishing places contribute 0" rule come back.
// The question is AT WHICH PRIME those recovered terms contribute.  The deficiency we must explain is
// purely log N.  If the recovered contributions sit at other primes, the zeta story cannot be the
// fix -- it would break the primes the pipeline already gets right.
AttachSpec("ShimuraQuotients.spec");

for b in [<10,3,-8,[3,12,30]>, <15,2,-7,[2,10,30]>, <6,5,-19,[10,15,30]>] do
    D := b[1]; N := b[2]; d := b[3]; ms := b[4];
    Ld := ShimuraCurveLattice(D,N);
    Q := ChangeRing(Ld`Q, Integers()); Qrat := ChangeRing(Q, Rationals());
    lam := ElementOfNorm(Ld`Q, -d, Ld`O, Ld`basis_L);
    c_Lplus := Content(lam);
    Lminus := Kernel(Transpose(Matrix(lam*Q)));
    Lplus := RSpaceWithBasis(Matrix(lam div c_Lplus));
    L := RSpaceWithBasis(IdentityMatrix(Integers(),3));
    L_quo, L_quo_map := L / (Lplus + Lminus);
    lam_rat := ChangeRing(lam, Rationals());
    Bm := BasisMatrix(Lminus);
    Delta := Determinant(Bm*Q*Transpose(Bm));
    dk := FundamentalDiscriminant(d);
    printf "\n===== X0^%o(%o) d=%o  (N = %o is %o in k) =====\n", D, N, d, N,
           (dk mod N eq 0) select "ramified" else
           (KroneckerSymbol(dk,N) eq 1 select "SPLIT" else "inert");
    for m in ms do
        recovered := AssociativeArray();   // prime -> how many recovered terms differentiate there
        still0 := 0; unchanged := 0;
        for mu_bar in L_quo do
            mu := mu_bar@@L_quo_map;
            c_mu := ((mu*Q, lam)/(lam*Q,lam));
            mu_minus := ChangeRing(mu, Rationals()) - c_mu*lam_rat;
            sqr_bd := m/(-d);
            for k in [Ceiling((-Sqrt(sqr_bd) - c_mu)*c_Lplus)..Floor((Sqrt(sqr_bd) - c_mu)*c_Lplus)] do
                x := (c_mu + k*c_Lplus^(-1))*lam_rat;
                mm := m - (x*Qrat,x)/2;
                if mm eq 0 then continue; end if;
                Sm := Sort(Setseq({p : p in PrimeDivisors(Delta)}
                                  join {p : p in PrimeDivisors(Numerator(mm))}));
                van := [];
                for p in Sm do
                    if LocalWhittakerAtOne(mm, p, mu_minus, Lminus, Q) eq 0 then Append(~van, p); end if;
                end for;
                vanNoN := [p : p in van | p ne N];       // after cancelling the zeta factor at N
                if #van eq #vanNoN then unchanged +:= 1; continue; end if;   // N did not vanish
                if #vanNoN eq 1 then
                    p := vanNoN[1];
                    if not IsDefined(recovered, p) then recovered[p] := 0; end if;
                    recovered[p] +:= 1;
                elif #vanNoN eq 0 then
                    if not IsDefined(recovered, N) then recovered[N] := 0; end if;
                    recovered[N] +:= 1;                   // N alone: this WOULD give a log N
                else
                    still0 +:= 1;
                end if;
            end for;
        end for;
        ks := Sort([p : p in Keys(recovered)]);
        printf "  m = %-4o recovered terms by prime: %-34o | still 0: %-4o unaffected: %o\n",
               m, [<p, recovered[p]> : p in ks], still0, unchanged;
    end for;
end for;
quit;
