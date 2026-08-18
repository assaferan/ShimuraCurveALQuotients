// Build X0^15(2) with the fitted ground-truth m=0 multipliers and run the INDEPENDENT
// verification (ModelVerification.VerifyModelSet): genus self-consistency, genus vs the
// Shimura-curve genus formula, Weil-polynomial divisibility across nested covers, and the
// trace-formula point count. None of these uses the Borcherds/Schofer CM machinery that
// produced the models, so passing them is a genuine cross-check of the 9 multipliers.

AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 0);
fh := Open("verify_out.txt", "w");
emit := procedure(s) fprintf fh, "%o\n", s; Flush(fh); end procedure;

D := 15; N := 2;
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};

eqns, ws := AllEquationsAboveCovers(star, curves : Prec := 100);
emit(Sprintf("#eqn keys = %o", #Keys(eqns)));

models := AssociativeArray();
nmodels := 0;
for k in Sort([kk : kk in Keys(eqns)]) do
    inner := eqns[k];
    if #Keys(inner) eq 0 then
        emit(Sprintf("  key %o (W = %o): NO equation", k, Sort([w : w in curves[k]`W])));
        continue;
    end if;
    W := Sort([w : w in curves[k]`W]);
    ents := [* *];
    for kk in Keys(inner) do
        C := inner[kk];
        try
            f, h := HyperellipticPolynomials(C);
            Append(~ents, <Genus(C), f, h>);
            emit(Sprintf("  W = %-28o genus %o :  y^2 + (%o) y = %o", W, Genus(C), h, f));
        catch e
            emit(Sprintf("  W = %-28o NON-HYPERELLIPTIC entry: %o", W, C));
        end try;
    end for;
    if #ents gt 0 then
        models[[Integers() | w : w in W]] := ents;
        nmodels +:= #ents;
    end if;
end for;
emit(Sprintf("assembled %o model entries over %o cover keys", nmodels, #Keys(models)));

emit("---- VerifyModelSet ----");
nchk, nfail := VerifyModelSet(models, D, N : Verbose := true);
emit(Sprintf("VerifyModelSet: %o checks, %o FAILURES", nchk, nfail));
emit("DONE");
exit;
