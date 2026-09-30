// Hypothesis: the Popa-layer discrepancy at gcd(n,N) > 1 is Sfast ignoring the UNIT condition
// on alpha (alpha^2 - t alpha + n = 0 mod Nu with alpha in (Z/NZ)^x).  For p | (n, N) the root
// alpha = 0 mod p is not a unit.  (a) Sfast vs brute-force S; (b) Popa's formula with the
// brute-force count against modular symbols on the same 208-tuple grid.
AttachSpec("ShimuraQuotients.spec");
_ := ClassNumberLU(-4);
import "TraceFormula.m" : Sfast, S, Lemma4_5, P, H, Phil, get_trace_hecke_AL;

// (a)
badS := 0; totS := 0; firstS := [];
for N in [12, 18, 20, 36, 45] do
    for n in [2, 3, 4, 5, 9] do
        if GCD(n, N) eq 1 then continue; end if;
        for t in [-Floor(SquareRoot(4*n))..Floor(SquareRoot(4*n))] do
            totS +:= 1;
            a := Sfast(N, 1, t, n); b := #S(N, 1, t, n);
            if a ne b then badS +:= 1; if #firstS lt 6 then Append(~firstS, <N, t, n, a, b>); end if; end if;
        end for;
    end for;
end for;
printf "(a) Sfast vs brute-force S at p | (n,N): %o of %o differ; first: %o\n", badS, totS, firstS;
// control: they agree when gcd(n, N) = 1
badC := 0; totC := 0;
for N in [12, 18, 20] do for n in [5, 7, 11, 25] do if GCD(n,N) ne 1 then continue; end if;
    for t in [-Floor(SquareRoot(4*n))..Floor(SquareRoot(4*n))] do
        totC +:= 1; if Sfast(N,1,t,n) ne #S(N,1,t,n) then badC +:= 1; end if;
    end for; end for; end for;
printf "    control gcd(n,N) = 1: %o of %o differ\n", badC, totC;

// (b) Popa with brute-force |S_N(t,n)|
function PopaBrute(N, k, n, Q)
    S1 := 0; Q_prime := N div Q; w := k - 2;
    max_abst := Floor(SquareRoot(4*Q*n)) div Q;
    for tQ in [-max_abst..max_abst] do
        t := tQ*Q;
        for u in Divisors(Q) do for u_prime in Divisors(Q_prime) do
            if ((4*n*Q-t^2) mod (u*u_prime)^2 eq 0) then
                C := #S(Q_prime, 1, t, Q*n) * Lemma4_5(Q_prime, u_prime, t^2 - 4*Q*n);
                S1 +:= P(k,t,Q*n)*H((4*Q*n-t^2) div (u*u_prime)^2)*C*MoebiusMu(u) / Q^(w div 2);
            end if;
        end for; end for;
    end for;
    S2 := 0;
    for d in Divisors(n*Q) do
        a := n*Q div d;
        if (a+d) mod Q eq 0 then S2 +:= Minimum(a,d)^(k-1)*Phil(N,Q,a,d) / Q^(w div 2); end if;
    end for;
    ret := -S1/2 - S2/2;
    if k eq 2 then ret +:= &+[n div d : d in Divisors(n) | GCD(d,N) eq 1]; end if;
    return ret;
end function;
bad := 0; total := 0; shown := 0;
for N in [12, 18, 20, 36, 45, 72, 108] do
    for Q in [Q : Q in Divisors(N) | GCD(Q, N div Q) eq 1] do
        for n in [2, 3, 4, 5, 9] do
            if GCD(n, N) eq 1 then continue; end if;
            for k in [2, 4] do
                total +:= 1;
                f := PopaBrute(N, k, n, Q); m := get_trace_hecke_AL(N, k, n, Q);
                if f ne m then
                    bad +:= 1;
                    if shown lt 25 then shown +:= 1; printf "STILL-BAD N=%o Q=%o n=%o k=%o : brute-S Popa %o, modsym %o\n", N, Q, n, k, f, m; end if;
                end if;
            end for;
        end for;
    end for;
end for;
printf "(b) Popa with brute-force S vs modular symbols at gcd(n,N)>1: %o of %o differ (was 80 with Sfast)\n", bad, total;
exit;
