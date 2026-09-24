// Does ANY support hypothesis fit X0^21(2)'s exact multipliers?  The level rule (N | r on the
// oo-block, nothing on the 0-block) is inconsistent there; test weaker hypotheses before concluding
// the rule is false -- in particular allow the 0-block index, and allow the full unrestricted system.
AttachSpec("ShimuraQuotients.spec");
D := 21; N := 2;
oracle := [2, 2, 4, 2, 2, 2, -4, -4, 0];
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
ks := Sort([k : k in Keys(fs)]);
M := 2*D*N;
foos := [qExpansionAtoo(fs[k], 1) : k in ks];
f0s  := [qExpansionAt0(fs[k], 1) : k in ks];
Roo := Sort(Setseq({ r : r in [1..-Valuation(foos[i])], i in [1..#ks] | Coefficient(foos[i],-r) ne 0 }));
R0  := Sort(Setseq({ j : j in [1..-Valuation(f0s[i])],  i in [1..#ks] | Coefficient(f0s[i],-j) ne 0 }));
printf "oo-indices %o ; 0-indices %o ; oracle %o\n", Roo, R0, oracle;

solve := procedure(name, Uidx, Vidx)
    nu := #Uidx; nv := #Vidx;
    A := Matrix(Rationals(), #ks, nu+nv,
        &cat[ [ Rationals() | Coefficient(foos[i], -r) : r in Uidx ]
          cat [ Rationals() | Coefficient(f0s[i], -j)  : j in Vidx ] : i in [1..#ks] ]);
    rhs := Matrix(Rationals(), #ks, 1, [Rationals() | oracle[i] : i in [1..#ks]]);
    ok, sol := IsConsistent(Transpose(A), Transpose(rhs));
    if not ok then
        printf "  %-28o unknowns=%o : INCONSISTENT\n", name, nu+nv;
        return;
    end if;
    kd := Dimension(Kernel(Transpose(A)));
    printf "  %-28o unknowns=%o : consistent, dim ker = %o, spare = %o\n",
           name, nu+nv, kd, #ks - (nu+nv-kd);
    s := Eltseq(sol);
    printf "        A_r = %o   B_j = %o\n",
           [<Uidx[t], s[t]> : t in [1..nu]], [<Vidx[t], s[nu+t]> : t in [1..nv]];
end procedure;

solve("level rule, oo only", [r : r in Roo | r mod N eq 0], []);
solve("level rule + all 0-block", [r : r in Roo | r mod N eq 0], R0);
solve("all oo, no 0-block", Roo, []);
solve("everything", Roo, R0);
quit;
