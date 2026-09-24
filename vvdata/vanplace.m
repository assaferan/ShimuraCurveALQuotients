// Mirror Kappa's own loop and report, for every (mu, x1) term, WHICH place vanishes.
// I checked the vanishing pattern only at mu = 0 before; but kappaminus is called with
// mu = gamma_- + mu_-, and for a nonzero coset the local Whittaker at N can differ.  If N ever
// becomes the vanishing place when N | m, the log N channel is open inside the existing structure
// and the code is mishandling it.  If it never does, the structure genuinely cannot produce it.
AttachSpec("ShimuraQuotients.spec");

for b in [<10,3,-8,[3,12,30]>, <15,2,-7,[2,10,30]>, <6,5,-19,[10,15,30]>] do
    D := b[1]; N := b[2]; d := b[3]; ms := b[4];
    Ld := ShimuraCurveLattice(D,N);
    Q := ChangeRing(Ld`Q, Integers()); Qrat := ChangeRing(Q, Rationals());
    lambda_v := ElementOfNorm(Ld`Q, -d, Ld`O, Ld`basis_L);
    c_Lplus := Content(lambda_v);
    Lplus := RSpaceWithBasis(Matrix(lambda_v div c_Lplus));
    Lminus := Kernel(Transpose(Matrix(lambda_v*Q)));
    L := RSpaceWithBasis(IdentityMatrix(Integers(),3));
    L_quo, L_quo_map := L / (Lplus + Lminus);
    lambda_rat := ChangeRing(lambda_v, Rationals());
    Bm := BasisMatrix(Lminus);
    Delta := Determinant(Bm*Q*Transpose(Bm));
    printf "\n===== X0^%o(%o) at d = %o : |L_quo| = %o, Delta = %o =====\n",
           D, N, d, #L_quo, Factorization(Delta);
    for m in ms do
        printf "  m = %-4o (N | m)\n", m;
        for mu_bar in L_quo do
            mu := mu_bar@@L_quo_map;
            c_mu_plus := ((mu*Q, lambda_v)/(lambda_v*Q,lambda_v));
            mu_plus := c_mu_plus*lambda_rat;
            mu_minus := ChangeRing(mu, Rationals()) - mu_plus;
            sqr_bd := m/(-d);
            lb := Ceiling((-Sqrt(sqr_bd) - c_mu_plus)*c_Lplus);
            ub := Floor((Sqrt(sqr_bd) - c_mu_plus)*c_Lplus);
            for k in [lb..ub] do
                x := (c_mu_plus + k*c_Lplus^(-1))*lambda_rat;
                mm := m - (x*Qrat,x)/2;
                if mm eq 0 then
                    printf "      mu=%-14o k=%-3o : INNER m=0\n", Eltseq(mu), k;
                    continue;
                end if;
                Sm := Sort(Setseq({p : p in PrimeDivisors(Delta)}
                                  join {p : p in PrimeDivisors(Numerator(mm))}));
                vals := [];
                for p in Sm do
                    Append(~vals, <p, LocalWhittakerAtOne(mm, p, mu_minus, Lminus, Q)>);
                end for;
                van := [v[1] : v in vals | v[2] eq 0];
                printf "      mu=%-14o k=%-3o m-Q(x)=%-6o W=%-30o vanish=%o%o\n",
                       Eltseq(mu), k, mm, vals, van,
                       (N in van) select "   <== N IS THE VANISHING PLACE" else "";
            end for;
        end for;
    end for;
end for;
quit;
