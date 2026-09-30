// Exact count, on the OLD Sfast (main), of mismatches vs brute force over the loop that
// tests/SfastUnitCondition.m actually runs -- so the test header can state it precisely.
AttachSpec("ShimuraQuotients.spec");
_ := ClassNumberLU(-4);
import "TraceFormula.m" : Sfast, S;
bad := 0; tot := 0;
for N in [4, 8, 9, 12, 16, 18, 20, 24, 27, 36, 45, 50, 72] do
    for u in Divisors(N) do
        for n in [2, 3, 4, 5, 6, 8, 9, 12, 18, 25] do
            for t in [-Floor(SquareRoot(4*n))..Floor(SquareRoot(4*n))] do
                if (t^2 - 4*n) mod u^2 ne 0 then continue; end if;
                tot +:= 1;
                if Sfast(N, u, t, n) ne #S(N, u, t, n) then bad +:= 1; end if;
            end for;
        end for;
    end for;
end for;
printf "OLD Sfast on the test loop: %o of %o counts differ from brute force\n", bad, tot;
exit;
