AttachSpec("ShimuraQuotients.spec");
_ := ClassNumberLU(-4);
import "TraceFormula.m" : TraceFormulaGamma0HeckeALNew, get_trace_hecke_AL;
bad := 0; total := 0;
for N in [12, 18, 20, 36, 45, 72, 108] do
    for Q in [Q : Q in Divisors(N) | GCD(Q, N div Q) eq 1] do
        for n in [1, 2, 3, 4, 5, 9] do
            for k in [2, 4] do
                f := TraceFormulaGamma0HeckeALNew(N, k, n, Q);
                m := get_trace_hecke_AL(N, k, n, Q : New);
                total +:= 1;
                if f ne m then
                    bad +:= 1;
                    printf "MISMATCH N=%o Q=%o n=%o k=%o gcd(n,N)=%o : formula %o, modsym %o\n", N, Q, n, k, GCD(n,N), f, m;
                end if;
            end for;
        end for;
    end for;
end for;
printf "%o mismatches of %o\n", bad, total;
exit;
