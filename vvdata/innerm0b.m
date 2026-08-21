// Yang eq (12): kappa_eta(m) = sum_mu sum_{x in eta_+ + mu_+ + L_+} kappa^-_{eta_-+mu_-}(m - <x,x>/2).
// By eq (11) the inner kappa^- at argument 0 is kappa_0(0) when eta_-+mu_- = 0 and 0 otherwise --
// and kappa_0(0) is exactly where Yang's Remark 15(2) puts the N-part  sum_{p|N/(N,d)} log p.
// So the log N enters the m > 0 terms through the INNER m = 0 of the x-sum.  The code's branch for
// that case does something ad hoc and `kappaminuszero`, which computes kappa_0(0), is dead code.
// COUNT the firings across ALL firing discriminants (not one, as I did before) and compare with A_m.
AttachSpec("ShimuraQuotients.spec");

targets := AssociativeArray();
targets[[15,2]] := [<2,-1>,<10,1>,<30,0>];
targets[[6,5]]  := [<10,3/2>,<15,1/2>,<30,3/2>];
targets[[10,3]] := [<3,1/2>,<12,1/2>,<30,3/2>];

for b in [[15,2],[6,5],[10,3]] do
    D := b[1]; N := b[2];
    Ld := ShimuraCurveLattice(D,N);
    Q := ChangeRing(Ld`Q, Integers()); Qrat := ChangeRing(Q, Rationals());
    Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
    Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};
    rat, quad := RationalandQuadraticCMPoints(star : coprime_to_level := false, bd := 4);
    ds := Sort(Setseq({a[1] : a in rat} join {a[1] : a in quad}));
    firing := [d : d in ds | FundamentalDiscriminant(d) mod N ne 0];
    printf "\n===== X0^%o(%o) : firing discs %o =====\n", D, N, firing;
    printf "  %-6o %-10o %o\n", "m", "target A_m", "inner-m=0 count per firing disc";
    for t in targets[[D,N]] do
        m := t[1]; want := t[2];
        counts := [];
        for d in firing do
            ok := true; lam := 0;
            try lam := ElementOfNorm(Ld`Q, -d, Ld`O, Ld`basis_L); catch e ok := false; end try;
            if not ok then Append(~counts, <d, "no lam">); continue; end if;
            c_Lplus := Content(lam);
            Lplus := RSpaceWithBasis(Matrix(lam div c_Lplus));
            Lminus := Kernel(Transpose(Matrix(lam*Q)));
            L := RSpaceWithBasis(IdentityMatrix(Integers(),3));
            L_quo, L_quo_map := L / (Lplus + Lminus);
            lam_rat := ChangeRing(lam, Rationals());
            cnt := 0;
            for mu_bar in L_quo do
                mu := mu_bar@@L_quo_map;
                c_mu := ((mu*Q, lam)/(lam*Q,lam));
                mu_minus := ChangeRing(mu, Rationals()) - c_mu*lam_rat;
                // eq (11): the inner kappa^- at 0 is nonzero only for the trivial coset
                if not IsCoercible(Lminus, [Rationals()!x : x in Eltseq(mu_minus)]) then continue; end if;
                sqr_bd := m/(-d);
                for k in [Ceiling((-Sqrt(sqr_bd)-c_mu)*c_Lplus)..Floor((Sqrt(sqr_bd)-c_mu)*c_Lplus)] do
                    x := (c_mu + k*c_Lplus^(-1))*lam_rat;
                    if (m - (x*Qrat,x)/2) eq 0 then cnt +:= 1; end if;
                end for;
            end for;
            Append(~counts, <d, cnt>);
        end for;
        printf "  %-6o %-10o %o\n", m, want, counts;
    end for;
end for;
quit;
