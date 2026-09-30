// Candidate fix for Sfast: for p | N with p | n, the roots of alpha^2 - t alpha + n mod p are
// 0 and t; only t can be a unit, and it is a simple root iff p does not divide t, lifting
// uniquely by Hensel.  So the local factor at such p is 1 if p ∤ t and 0 if p | t, for every
// exponent.  Check the patched count against brute force on a wide grid (u > 1, p = 2, prime
// powers), then Popa's formula with the patched count against modular symbols.
AttachSpec("ShimuraQuotients.spec");
_ := ClassNumberLU(-4);
import "TraceFormula.m" : Sfast, S, Lemma4_5, P, H, Phil, get_trace_hecke_AL;

function SfastFixed(N, u, t, n)
    fac := Factorization(N*u);
    num_sols := 1;
    y := t^2-4*n;
    for f in fac do
        p,e := Explode(f);
        if n mod p eq 0 then                       // p | (n, N): the unit condition bites
            if t mod p eq 0 then return 0; end if;
            continue;                              // exactly one unit root, lifts uniquely
        end if;
        if (y eq 0) then num_sols *:= p^(e div 2); continue; end if;
        e_y := Valuation(y, p);
        y_0 := y div p^e_y;
        if p eq 2 then
            if (e_y le e + 1) and IsOdd(e_y) then return 0; end if;
            if (e le e_y-2) then num_sols *:= 2^(e div 2);
            elif (e_y le e - 1) and IsEven(e_y) then num_sols *:= 2^(e_y div 2 - 1)*(1 + KroneckerSymbol(y_0, 2))*(1 + KroneckerSymbol(-1, y_0));
            elif (e eq e_y) and IsEven(e_y) then num_sols *:= 2^(e_y div 2 - 1)*(1 + KroneckerSymbol(-1, y_0));
            elif (e eq e_y - 1) and IsEven(e_y) then num_sols *:= 2^(e_y div 2 - 1);
            else error "Not Implemented!"; end if;
        else
            if e_y lt e then
                if IsOdd(e_y) then return 0; end if;
                if not IsSquare(Integers(p)!y_0) then return 0; end if;
                num_sols *:= p^(e_y div 2) * 2;
            else
                num_sols *:= p^(e div 2);
            end if;
        end if;
    end for;
    return Integers()!num_sols div u;
end function;

// (a) count vs brute force, including u > 1 and prime powers
bad := 0; tot := 0; first := [];
for N in [4, 8, 9, 12, 16, 18, 20, 24, 27, 36, 45, 50, 72, 108] do
    for u in Divisors(N) do
        for n in [2, 3, 4, 5, 6, 8, 9, 12, 18, 25, 27] do
            for t in [-Floor(SquareRoot(4*n))..Floor(SquareRoot(4*n))] do
                if (t^2 - 4*n) mod u^2 ne 0 then continue; end if;
                tot +:= 1;
                a := SfastFixed(N, u, t, n); b := #S(N, u, t, n);
                if a ne b then bad +:= 1; if #first lt 8 then Append(~first, <N, u, t, n, a, b>); end if; end if;
            end for;
        end for;
    end for;
end for;
printf "(a) SfastFixed vs brute force: %o of %o differ; first %o\n", bad, tot, first;

// (b) Popa with SfastFixed vs modular symbols, wider grid
function PopaFixed(N, k, n, Q)
    S1 := 0; Q_prime := N div Q; w := k - 2;
    max_abst := Floor(SquareRoot(4*Q*n)) div Q;
    for tQ in [-max_abst..max_abst] do
        t := tQ*Q;
        for u in Divisors(Q) do for u_prime in Divisors(Q_prime) do
            if ((4*n*Q-t^2) mod (u*u_prime)^2 eq 0) then
                C := SfastFixed(Q_prime, 1, t, Q*n) * Lemma4_5(Q_prime, u_prime, t^2 - 4*Q*n);
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
function ModSymTrace(N, k, n, Q)
    M := ModularSymbols(N, k, 1); C := CuspidalSubspace(M);
    if Dimension(C) eq 0 then return 0; end if;
    B := Matrix(Basis(VectorSpace(C)));
    return Trace(Solution(B, B*(HeckeOperator(M, n) * AtkinLehner(M, Q))));
end function;
bad := 0; tot := 0; shown := 0;
for N in [8, 12, 16, 18, 20, 24, 27, 36, 45, 50, 72, 108] do
    for Q in [Q : Q in Divisors(N) | GCD(Q, N div Q) eq 1] do
        for n in [2, 3, 4, 5, 6, 8, 9, 12, 25] do
            if GCD(n, N) eq 1 then continue; end if;
            for k in [2, 4, 6] do
                tot +:= 1;
                f := PopaFixed(N, k, n, Q); m := ModSymTrace(N, k, n, Q);
                if f ne m then bad +:= 1; if shown lt 20 then shown +:= 1; printf "BAD N=%o Q=%o n=%o k=%o : fixed Popa %o, modsym %o\n", N, Q, n, k, f, m; end if; end if;
            end for;
        end for;
    end for;
end for;
printf "(b) Popa with SfastFixed vs modular symbols at gcd(n,N)>1: %o of %o differ\n", bad, tot;
exit;
