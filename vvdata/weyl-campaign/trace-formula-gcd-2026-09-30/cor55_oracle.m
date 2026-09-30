// Isolate [Assaf, Cor. 5.5] (and the code's implementation of it) from Popa's formula:
// run the SAME newform recursion as TraceFormulaGamma0HeckeALNew / ...NewSmaller, but with
// the full-space traces taken from modular symbols instead of Popa.  If this matches the
// modular-symbol newform traces at gcd(n, N) > 1, the corollary and its implementation are
// fine and the whole discrepancy is in the Popa layer.
AttachSpec("ShimuraQuotients.spec");
_ := ClassNumberLU(-4);
import "TraceFormula.m" : alpha, IsRelevantNprime, get_trace_hecke_AL;
// Magma's HeckeOperator on a SUBSPACE dies at (N,k,n) = (6,4,2) ("incompatible coefficient
// rings"), so take the operators on the ambient space and restrict, as tests/trace_formula.m does.
function ModSymTrace(N, k, n, Q : New := false)
    M := ModularSymbols(N, k, 1);
    C := CuspidalSubspace(M);
    if New then C := NewSubspace(C); end if;
    if Dimension(C) eq 0 then return 0; end if;
    T := HeckeOperator(M, n) * AtkinLehner(M, Q);
    B := Matrix(Basis(VectorSpace(C)));
    return Trace(Solution(B, B*T));
end function;


forward NewOracle;
function SmallerOracle(N, k, n, Q)
    trace := 0;
    n_Q := GCD(n, Q);
    n_NQ := n div n_Q;
    n_p_NQs := [n_p_NQ : n_p_NQ in Divisors(n_NQ) | IsSquare(n_p_NQ)];
    for n_p_Q in Divisors(n_Q) do
        for d_Q in Divisors(n_p_Q) do
            for n_p_NQ in n_p_NQs do
                n_p := n_p_Q * n_p_NQ;
                if (n_p eq 1) then continue; end if;
                _, d_NQ := IsSquare(n_p_NQ);
                d := d_Q * d_NQ;
                Q_primes := [Q_p : Q_p in Divisors(Q) | IsRelevantNprime(Q_p, Q, d_Q, n_p_Q, n_Q)];
                NQ_primes := [NQ_p : NQ_p in Divisors(N div Q) | (GCD(d_NQ, NQ_p) eq 1) and ((N div (Q*NQ_p)) mod d_NQ eq 0)];
                N_primes := [Q_p * NQ_p : Q_p in Q_primes, NQ_p in NQ_primes];
                weights := [#[x : x in Divisors(GCD(N div Q, N div N_p)) | GCD(x,n) eq 1] : N_p in N_primes];
                traces := [NewOracle(N_p, k, n div n_p, GCD(N_p, Q)) : N_p in N_primes];
                term := &+[Integers() | weights[i]*traces[i] : i in [1..#N_primes]];
                term *:= n_p^(k div 2) div d;
                term *:= MoebiusMu(d);
                trace +:= term;
            end for;
        end for;
    end for;
    return trace;
end function;

function NewOracle(N, k, n, Q)
    trace := 0;
    for N_prime in Divisors(N) do
        a := alpha(Q, n, N div N_prime);
        if a eq 0 then continue; end if;
        Q_prime := GCD(N_prime, Q);
        term := ModSymTrace(N_prime, k, n, Q_prime);        // modular symbols, full space
        term -:= SmallerOracle(N_prime, k, n, Q_prime);
        trace +:= a*term;
    end for;
    return trace;
end function;

bad := 0; total := 0;
for N in [12, 18, 20, 36, 45, 72, 108] do
    for Q in [Q : Q in Divisors(N) | GCD(Q, N div Q) eq 1] do
        for n in [2, 3, 4, 5, 9] do
            if GCD(n, N) eq 1 then continue; end if;
            for k in [2, 4] do
                total +:= 1;
                f := NewOracle(N, k, n, Q);
                m := ModSymTrace(N, k, n, Q : New);
                if f ne m then
                    bad +:= 1;
                    printf "COR55-BAD N=%o Q=%o n=%o k=%o : recursion-on-modsym %o, modsym new %o\n", N, Q, n, k, f, m;
                end if;
            end for;
        end for;
    end for;
end for;
printf "Cor 5.5 recursion on modular-symbol full traces: %o mismatches of %o at gcd(n,N)>1\n", bad, total;
exit;
