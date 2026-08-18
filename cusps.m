// Test the pipeline's standing assumption that the Borcherds forms have poles ONLY at the cusps
// oo and 0. Uses Ligozat's order formula for eta quotients: for f = prod_{delta|M} eta(delta*tau)^{r_delta},
// the order at the cusp with denominator c | M is
//     ord_c(f) = (M / (24 gcd(c^2,M))) * sum_{delta|M} gcd(c,delta)^2 r_delta / delta.
// Sanity: at c = M this gives (sum delta*r_delta)/24, exactly the convention in EtaQuotient.m:93.
// The code's forms are linear COMBINATIONS of eta quotients, so ord >= min over terms; a min >= 0
// therefore PROVES no pole (a min < 0 only indicates a possible pole, since terms can cancel).
//
// The correlation under test: on X0^10(11) the m=0 multiplier is correct (0) exactly for forms
// -2, -1, 14 (the hauptmoduls) and wrong for the cover forms 9,10,11,12,13,15.

AttachSpec("ShimuraQuotients.spec");
fh := Open("cusps_out.txt", "w");
emit := procedure(s) fprintf fh, "%o\n", s; Flush(fh); end procedure;

ord_at := function(M, ds, r, c)
    return (Rationals()!M/(24*GCD(c^2, M))) * (&+[Rationals() | GCD(c, d)^2 * r[i] / d : i->d in ds]);
end function;

emit("START");
for base in [<15,2>, <10,11>] do
    D, N := Explode(base);
    emit(Sprintf("=========== X0^%o(%o) ===========", D, N));
    Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
    Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};
    fs := BorcherdsForms(star, curves : Prec := 100);
    ks := Sort([k : k in Keys(fs)]);
    for k in ks do
        eta := fs[k];
        R := Parent(eta);
        M := R`M; ds := R`ds;
        exps := [r : r in Exponents(eta)];
        bad := [];
        for c in Divisors(M) do
            if c eq 1 or c eq M then continue; end if;   // cusps 0 and oo are allowed poles
            mn := Minimum([ord_at(M, ds, r, c) : r in exps]);
            if mn lt 0 then Append(~bad, <c, mn>); end if;
        end for;
        o_oo := Minimum([ord_at(M, ds, r, M) : r in exps]);
        o_0  := Minimum([ord_at(M, ds, r, 1) : r in exps]);
        emit(Sprintf("  form[%-3o] M=%-4o #terms=%-3o  min-ord: oo=%-8o 0=%-8o | other cusps with min<0: %o",
                     k, M, #exps, o_oo, o_0, IsEmpty(bad) select "NONE" else Sprint(bad)));
    end for;
end for;
emit("DONE");
