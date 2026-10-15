// For each open genus-3 quotient X and each residual involution w with X/w elliptic: the newform
// of E = X/w through the traces of T_n on the D-new W'-fixed space, its isogeny class in Cremona's
// tables, and the rank and torsion of E(Q).
AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
cs := eval Read("data/curves_after_UpdateCurves8.dat");
open := [X : X in cs | not assigned X`IsSubhyp and X`g eq 3];
DB := CremonaDatabase();
for X in open do
    D := X`D; N := X`N; W := X`W; DN := D*N;
    full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
    seen := {};
    for w in full diff W do
        coset := {AtkinLehnerMul(w, u, DN) : u in W};
        if coset in seen then continue; end if; Include(~seen, coset);
        Wb := W join coset;
        if GenusShimuraCurveQuotient(D, N, Wb) ne 1 then continue; end if;
        ps := [p : p in PrimesUpTo(60) | DN mod p ne 0];
        ap := [TraceDNewALFixed(D, N, 2, p, Wb) : p in ps];
        found := [];
        for E in EllipticCurves(DB, DN) do
            if [TraceOfFrobenius(E, p) : p in ps] eq ap then Append(~found, E); end if;
        end for;
        if #found eq 0 then printf "X_0(%o,%o)/%o w=%o: E of conductor %o not found (a_p %o)\n", D, N, Sort(Setseq(W)), w, DN, ap[1..4]; continue; end if;
        E := found[1];
        printf "X_0(%o,%o)/%o w=%o: E = %o (class of %o curves), rank %o, torsion %o\n", D, N, Sort(Setseq(W)), w, CremonaReference(E), #found, Rank(E), #TorsionSubgroup(E);
    end for;
end for;
quit;
