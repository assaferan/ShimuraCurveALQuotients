// Split the gcd(n,N) > 1 discrepancy: is it already in Popa's FULL-space formula
// (TraceFormulaGamma0HeckeAL vs modular symbols on S_k(N)), or only in the newform recursion?
AttachSpec("ShimuraQuotients.spec");
_ := ClassNumberLU(-4);
import "TraceFormula.m" : TraceFormulaGamma0HeckeAL, TraceFormulaGamma0HeckeALNew, get_trace_hecke_AL;
badF := 0; badN := 0; total := 0; badN_onlyF_ok := 0;
for N in [12, 18, 20, 36, 45, 72, 108] do
    for Q in [Q : Q in Divisors(N) | GCD(Q, N div Q) eq 1] do
        for n in [2, 3, 4, 5, 9] do
            if GCD(n, N) eq 1 then continue; end if;
            for k in [2, 4] do
                total +:= 1;
                f := TraceFormulaGamma0HeckeAL(N, k, n, Q);
                m := get_trace_hecke_AL(N, k, n, Q);
                fn := TraceFormulaGamma0HeckeALNew(N, k, n, Q);
                mn := get_trace_hecke_AL(N, k, n, Q : New);
                if f ne m then badF +:= 1; end if;
                if fn ne mn then badN +:= 1; if f eq m then badN_onlyF_ok +:= 1; end if; end if;
                if f ne m or fn ne mn then
                    printf "N=%o Q=%o n=%o k=%o  FULL formula %o modsym %o %o | NEW formula %o modsym %o %o\n",
                        N, Q, n, k, f, m, f eq m select "ok" else "BAD", fn, mn, fn eq mn select "ok" else "BAD";
                end if;
            end for;
        end for;
    end for;
end for;
printf "gcd(n,N)>1 tuples: %o; full-space Popa wrong: %o; newform wrong: %o (of which full-space right: %o)\n",
    total, badF, badN, badN_onlyF_ok;
exit;
