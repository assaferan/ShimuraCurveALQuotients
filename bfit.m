// SOLVE for the weight-3/2 Eisenstein coefficients b^{eta*}_eta(r) instead of guessing them.
//
// Borcherds' obstruction, with G the Eisenstein series whose constant term is e_eta* + e_-eta*:
//     mult_f = (1/2) c_eta*(0) = -(1/4) sum_{eta, r>0} c_eta(-r) b_eta(r).
// By [GY, Lemma 24] (verified numerically to 1e-103), c_eta(-r) is Coefficient(foo,-r) on eta = 0 for
// integer r, plus Coefficient(f0,-Mr) on every eta in the bucket.  So mult_f is a linear functional
//     mult_f = sum_r A_r * Coefficient(foo,-r)  +  sum_j B_j * Coefficient(f0,-j)
// with A_r = -b_0(r)/4  and  B_j = -(1/4) sum_{eta in bucket j} b_eta(j/M), and the oracle's exact
// multipliers give one equation per form.  This script reports the system, solves it, and prints the
// solution alongside what m0_multiplier's own local assembly produces at the same indices.
AttachSpec("ShimuraQuotients.spec");

bases := [[15,2],[6,5],[10,3]];
oracle := AssociativeArray();
// keys in the order Sort(Keys(BorcherdsForms)) = [-2,-1,9,10,11,12,13,14,15]
oracle[[15,2]] := [2, 4, 0, 0, 4, 2, 4, -2, 2];
oracle[[6,5]]  := [0, 0, 6, 3, 12, 9, 3, 0, 3];
oracle[[10,3]] := [0, 0, 6, 3, 3, 3, 3, 0, 6];

for b in bases do
    D := b[1]; N := b[2];
    Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
    Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};
    fs := BorcherdsForms(star, curves : Prec := 100);
    ks := Sort([k : k in Keys(fs)]);
    M := IsOdd(D*N) select 4*D*N else 2*D*N;
    Ld := ShimuraCurveLattice(D, N);
    printf "\n======================== X0^%o(%o)  M = %o ========================\n", D, N, M;

    foos := [qExpansionAtoo(fs[k], 1) : k in ks];
    f0s  := [qExpansionAt0(fs[k], 1) : k in ks];
    Roo := Sort(Setseq({ r : r in [1..-Valuation(foos[i])], i in [1..#ks] | Coefficient(foos[i],-r) ne 0 }));
    R0  := Sort(Setseq({ j : j in [1..-Valuation(f0s[i])],  i in [1..#ks] | Coefficient(f0s[i],-j) ne 0 }));
    printf "occurring oo-indices r : %o\n", Roo;
    printf "occurring  0-indices j : %o\n", R0;
    printf "of these, N | r : %o   and  N | j : %o\n",
           [r : r in Roo | r mod N eq 0], [j : j in R0 | j mod N eq 0];

    // Unknowns, under the LEVEL RULE (support on N | index).  Without a support restriction the
    // system is hopelessly underdetermined (9 equations); the rule is what makes it finite.
    Uidx := [r : r in Roo | r mod N eq 0];
    Vidx := [j : j in R0  | j mod N eq 0];
    nu := #Uidx; nv := #Vidx;
    printf "unknowns: %o (oo) + %o (0) = %o, equations: %o\n", nu, nv, nu+nv, #ks;

    Amat := Matrix(Rationals(), #ks, nu+nv,
        &cat[ [ Rationals() | Coefficient(foos[i], -r) : r in Uidx ]
          cat [ Rationals() | Coefficient(f0s[i], -j)  : j in Vidx ] : i in [1..#ks] ]);
    rhs := Matrix(Rationals(), #ks, 1, [Rationals() | oracle[b][i] : i in [1..#ks]]);
    ok, sol := IsConsistent(Transpose(Amat), Transpose(rhs));
    if not ok then
        printf "*** INCONSISTENT under the level rule -- the rule fails on this base ***\n";
        continue;
    end if;
    ker := Kernel(Transpose(Amat));
    printf "CONSISTENT.  solution space dimension = %o  (spare conditions = %o)\n",
           Dimension(ker), #ks - (nu+nv-Dimension(ker));
    s := Eltseq(sol);
    printf "  A_r  (r, weight, implied b_0(r) = -4A) :\n";
    for t->r in Uidx do
        printf "      r = %-5o A = %-8o b = %o\n", r, s[t], -4*s[t];
    end for;
    printf "  B_j  (j, weight, implied sum_bucket b = -4B) :\n";
    for t->j in Vidx do
        printf "      j = %-5o (r = %o) B = %-8o b = %o\n", j, (Rationals()!j)/M, s[nu+t], -4*s[nu+t];
    end for;
    if Dimension(ker) gt 0 then
        printf "  kernel basis (directions the data cannot see):\n";
        for v in Basis(ker) do printf "      %o\n", Eltseq(v); end for;
    end if;
end for;
quit;
