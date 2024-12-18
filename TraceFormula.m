// This file implements the trace formulas from [Popa] and [Assaf]
import "Caching.m" : cached_traces, SetCache, GetCache, class_nos;

function Sfast(N, u, t, n)
// Returns |S_N(u,t,n)| as defined in [Popa, 4.1]
// Recal S_N(u,t,n) = { \alpha in (Z/NZ)^{\times} : \alpha^2 - t \alpha + n = 0 (mod Nu)}
// It is called Sfast because it does not produce the set, buut merely counts using straight-forward formulas coming from 
// local-global principle and Hensel lifting.
    fac := Factorization(N*u);
    num_sols := 1;
    y := t^2-4*n;
    for f in fac do
        p,e := Explode(f);
        if (y eq 0) then
	        num_sols *:= p^(e div 2);
            continue;
        end if;
        e_y := Valuation(y, p);
        y_0 := y div p^e_y;
        // following [Assaf]
        if p eq 2 then
            if (e_y le e + 1) and IsOdd(e_y) then
                return 0;
            end if;
        if (e le e_y-2) then
            num_sols *:= 2^(e div 2);
        elif (e_y le e - 1) and IsEven(e_y) then
            num_sols *:= 2^(e_y div 2 - 1)*(1 + KroneckerSymbol(y_0, 2))*(1 + KroneckerSymbol(-1, y_0));
        elif (e eq e_y) and IsEven(e_y) then
            num_sols *:= 2^(e_y div 2 - 1)*(1 + KroneckerSymbol(-1, y_0));
        elif (e eq e_y - 1) and IsEven(e_y) then
            num_sols *:= 2^(e_y div 2 - 1);
        else
            Error("Not Implemented!");
        end if;
        else
            if e_y lt e then
                if IsOdd(e_y) then
                    return 0;
                end if;
                is_sq := IsSquare(Integers(p)!y_0);
                if not is_sq then
                    return 0;
                end if;
                num_sols *:= p^(e_y div 2) * 2;
            else
                num_sols *:= p^(e div 2);  
            end if;
        end if;
    end for;
    return Integers()!num_sols div u;
end function;

function S(N, u, t, n)
// Returns S_N(u,t,n) as in [Popa, 4.1]. 
// The size of the result should be Sfast(N,u,t,n).
    assert N mod u eq 0;
    assert (t^2 - 4 * n) mod u^2 eq 0;
    return [x : x in [0..N-1] | GCD(x,N) eq 1 and (x^2 - t*x + n) mod (N*u) eq 0];
end function;

function phi1(N)
// Returns |P^1(Z/NZ)|, denoted phi_1(N) in [Popa, 2.3]
    primes := [f[1] : f in Factorization(N)];
    if IsEmpty(primes) then return N; end if;
    return Integers()!(N * &*[1 + 1/p : p in primes]);
end function;

function B(N, u, t, n)
// Returns |B_{N,1}(u,t,n)|, as defined in [Popa, 2.3] (with trivial Nebentypus)
// assert Sfast(N,u,t,n) eq #S(N,u,t,n);
    return #S(N,u,t,n) * phi1(N) div phi1(N div u);
  // return Sfast(N,u,t,n) * phi1(N) div phi1(N div u);  
end function;

function C(N, M)
// Returns C_{N,1}(M) as defined in [Popa, 2.3], for a matrix M
    a, b, c, d := Explode(M);
    G_M := GCD(c, d-a, b);
    return B(N, GCD(G_M, N), Trace(M), Determinant(M));
end function;

function C(N, u, t, n)
// Returns C_{N,1}(u,t,n) as defined in [Popa, (2.18)]
    return &+[B(N, u div d, t, n) * MoebiusMu(d) : d in Divisors(u)];
end function;

function Lemma4_5(N, u, D)
// Returns C_N(u,D) as defined in [Popa, Lemma 4.5], with minor correction as mentioned in [Assaf]
    assert N mod u eq 0;
    assert D mod (u^2) eq 0;
    ret := 1;
    fac := Factorization(N);
    for pa in fac do
        p, a := Explode(pa);
        i := Valuation(u, p);
        if (i eq 0) then continue; end if;
        b := Valuation(D, p);
        if IsOdd(p) then
            if (i eq a) then ret *:= p^(Ceiling(a/2)); continue; end if;
            if ((i le b-a) and IsEven(a-i)) then
                ret *:= (p^Ceiling(i/2) - p^(Ceiling(i/2)-1));
            elif ((i eq b-a+1) and IsEven(a-i)) then
                ret *:= - p^(Ceiling(i/2)-1);
            elif ((i eq b-a+1) and IsOdd(a-i)) then
                ret *:= p^Floor(i/2) * LegendreSymbol(D div p^b, p);
            else
                return 0;
            end if;
        else // p = 2
            if (i eq a) then
                if (b ge 2*a+2) or ((b eq 2*a) and (((D div 2^b) mod 4) eq 1)) then 
                    ret *:= p^(Ceiling(a/2));
                elif (b eq 2*a+1) or ((b eq 2*a) and (((D div 2^b) mod 4) eq 3)) then
                    ret *:= -p^(Ceiling(a/2)-1);
                end if;
                continue;
            end if;
            if ((i le b-a-2) and IsEven(a-i)) then
                // print "case 1";
                ret *:= p^(Ceiling(i/2)-1);
            elif ((i eq b-a-1) and IsEven(a-i)) then
                // print "case 2";
                ret *:= - p^(Ceiling(i/2)-1);
            elif ((i eq b-a) and IsEven(a-i)) then
                // print "case 3";
                ret *:= p^(Ceiling(i/2)-1) * KroneckerSymbol(-4,D div p^b);
            elif ((i eq b-a+1) and IsOdd(a-i) and ((D div p^b) mod 4 eq 1) ) then
                // print "case 4";
                ret *:= p^Floor(i/2) * KroneckerSymbol(D div p^b, p);
            else
                // print "returning 0";
                return 0;
            end if;
        end if;
    end for;
    return ret;
end function;

function Cfast(N, u, t, n)
// Returns C_N(u,t,n), computed using [Popa, Lemma 4.5]
    // S := [x : x in [0..N-1] | (GCD(x,N) eq 1) and (((x^2 - t*x + n) mod N) eq 0)];
    // nS1 := #S(N, 1, t, n);
    nS2 := Sfast(N, 1, t, n);
    // assert nS1 eq nS2;
    return nS2 * Lemma4_5(N, u, t^2 - 4*n);
end function;

function Hurwitz(n)
// Returns the Hurwitz class number of n
    if n eq 0 then
	    return -1/12;
    end if;
    t_sum := 0;
    
    for d in Divisors(n) do
        is_sq, f := IsSquare(d);
        if is_sq then
            D := -n div f^2;
            if D mod 4 in [0,1] then
                O := QuadraticOrder(BinaryQuadraticForms(D));
                h := #PicardGroup(O);
                w := #TorsionSubgroup(UnitGroup(O));
                t_sum +:= (h/w);
            end if;
        end if;
    end for;
    return 2*t_sum;
end function;

function H(n)
// Returns the Kronecker-Hurwitz class number of n, as extended by Zagier, see [Popa, 2.2]
    if n lt 0 then
        is_sq, u := IsSquare(-n);
        return (is_sq select -u/2 else 0);
    end if;
    if n eq 0 then
        return -1/12;
    end if;
    if n mod 4 in [1,2] then
        return 0;
    end if;

    ret := 0;
    for d in Divisors(n) do
        if IsSquare(d) and  (n div d) mod 4 in [0,3] then
            b, v := GetCache(-n div d, class_nos);
            if not b then
                v := ClassNumber(-n div d);
                SetCache(-n div d, v, class_nos); 
            end if;
            ret +:= v;
        end if;
    end for;
    // ret := &+[ClassNumber(-n div d) : d in Divisors(n)
    // 		       | IsSquare(d) and (n div d) mod 4 in [0,3] ];
    if IsSquare(n) and IsEven(n) then
        ret -:= 1/2;
    end if;
    if n mod 3 eq 0 and IsSquare(n div 3) then
        ret -:= 2/3;
    end if;
    return ret;
end function;

function Phi(N, a, d)
// Returns Phi_{N,1}(a,d) as defined in [Popa, Theorem 3]
    ret := 0;
    for r in Divisors(N) do
        s := N div r;
        // scalar := r eq s select 1/2 else 1;
        g := GCD(r,s);
        if GCD(N, a-d) mod g eq 0 then
            alpha := CRT([a,d],[r,s]);
            if (GCD(alpha,N) eq 1) then
                ret +:= EulerPhi(g);
            end if;
        end if;
    end for;
    return ret;
end function;

// Gegenbauer polynomials
function P(k, t, m)
// Returns p_k(t,m) - the value of the k-th Gegenbauer polynomial at (t,n), as defined in [Popa] before Theorem 1
    R<x> := PowerSeriesRing(Rationals(), k-1);
    return Coefficient((1 - t*x+m*x^2)^(-1), k-2);
end function;

// This should yield -2*(A1 + A2)
function S1Popa(n, N, k)
// Returns the first summand in [Popa, Theorem 3]
    S1 := 0;
    max_abst := Floor(SquareRoot(4*n));
    for t in [-max_abst..max_abst] do
        for u in Divisors(N) do
            if ((4*n-t^2) mod u^2 eq 0) then
                // print "u = ", u, "t = ", t;
                S1 +:= P(k,t,n)*H((4*n-t^2) div u^2)*C(N,u,t,n);
                // assert H((4*n-t^2) div u^2) eq Hurwitz((4*n-t^2) div u^2);
                // S1 +:= P(k,t,n)*Hurwitz((4*n-t^2) div u^2)*C(N,u,t,n);
                // print "S1 = ", S1;
            end if;
        end for;
    end for;
    return S1;
end function;

// This should yield -2*A3
function S2Popa(n, N, k)
// Returns the second summand in [Popa, Theorem 3]
    S2 := 0;
    for d in Divisors(n) do
        a := n div d;
        S2 +:= Minimum(a,d)^(k-1)*Phi(N,a,d);
    end for;
    return S2;
end function;

function TraceFormulaGamma0(n, N, k)
// Returns the Trace of T_n on S_k(N) using the formula in [Popa, Theorem 3]
    S1 := S1Popa(n,N,k);
    S2 := S2Popa(n,N,k);
    ret := -S1 / 2 - S2 / 2;
    if k eq 2 then
	    ret +:= &+[n div d : d in Divisors(n) | GCD(d,N) eq 1];
    end if;
    return ret;
end function;

function Phil(N, l, a, d)
// Returns Phi_{N,l}(a,d) as defined in [Popa, Theorem 4]
    l_prime := N div l;
    ret := 0;
    for r in Divisors(l_prime) do
        s := l_prime div r;
        g := GCD(r,s);
        if (((a-d) mod g eq 0) and (GCD(a,r) eq 1) and (GCD(d,s) eq 1)) then
            ret +:= EulerPhi(g);
        end if;
    end for;
    return EulerPhi(l) * ret / l;
end function;

function alpha(n)
// Returns the arithmetic function alpha defined in [SZ88, p.117]
    fac := Factorization(n);
    ret := 1;
    for fa in fac do
        if (fa[2] in [1,2]) then ret := -ret; end if;
        if (fa[2] ge 4) then return 0; end if;
    end for;
    return ret;
end function;

// Now this seems to work - tested for N <= 100, even k, 2<=k<=12 and 1 <= n <= 10, n = N
// This formula follows Popa - On the Trace Formula for Hecke Operators on Congruence Subgroups, II
// Theorem 4. 
// (Also appears in Skoruppa-Zagier, but this way of stating the formula was easier to work with).
function TraceFormulaGamma0AL(n, N, k)
// Returns the trace of T_n \circ W_N on S_k(N) using [Popa, Theorem 4]
    if (n eq 0) then return 0; end if; // for compatibility with q-expansions
    S1 := 0;
    //  max_abst := Floor(SquareRoot(4*N*n));
    max_abst := Floor(SquareRoot(4*N*n)) div N;
    // ts := [t : t in [-max_abst..max_abst] | t mod N eq 0];
    // for t in ts do
    for tN in [-max_abst..max_abst] do
        t := tN*N;
        for u in Divisors(N) do
            if ((4*n*N-t^2) mod u^2 eq 0) then
                S1 +:= P(k,t,N*n)*H((4*N*n-t^2) div u^2)*C(1,1,t,N*n)
                        *MoebiusMu(u) / N^(k div 2-1);
            end if;
        end for;
    end for;
    S2 := 0;
    for d in Divisors(n*N) do
        a := n*N div d;
        if (a+d) mod N eq 0 then 
            S2 +:= Minimum(a,d)^(k-1)*EulerPhi(N) / N^(k div 2);
        end if;
    end for;
    ret := -S1 / 2 - S2 / 2;
    if k eq 2 then
        ret +:= &+[n div d : d in Divisors(n) | GCD(d,N) eq 1];
    end if;
    return ret;
end function;

function TraceFormulaGamma0ALTrivialNew(N, k)
// Returns the trace of W_N on S_k(N)^{new} by Mobius inversion
    ms := [d : d in Divisors(N) | N mod d^2 eq 0];
    trace := &+[Integers() | MoebiusMu(m)*TraceFormulaGamma0AL(1, N div m^2, k) : m in ms];
    return trace;
end function;

// At the moment only works for Hecke operators at primes
function TrivialContribution(N, k, p)
// Returns the contribution to the trace of T_p \circ W_N on S_k(N)^{new} when p is prime.
// coming from traces of W_{N_prime} on S_k(N_prime)^{new} as in [Assaf, Cor. 5.14]
    assert IsPrime(p);
    N_primes_2 := [d : d in Divisors(N) | IsSquare(N*p div d) and ((N div d) mod p^3 ne 0) and (d mod p ne 0)];
    // trace_2 := &+[Integers() | get_trace(N_prime, k, 1 : New) : N_prime in N_primes_2];
    trace_2 := &+[Integers() | TraceFormulaGamma0ALTrivialNew(N_prime, k) : N_prime in N_primes_2];
    trace_3 := 0;
    if (N mod p eq 0) then
        N_primes_3 := [d : d in Divisors(N div p) | IsSquare(N div (d*p))];
        trace_3 := &+[Integers() | TraceFormulaGamma0ALTrivialNew(N_prime, k) :  N_prime in N_primes_3];
    end if;
    return p^(k div 2) * trace_3  - p^(k div 2 - 1)*trace_2;
end function;

// At the moment only works for Hecke operators at primes
function TraceFormulaGamma0ALNew(p, N, k)
// Returns the trace of T_p \circ W_N on S_k(N)^{new} when p is prime, using [Assaf, Cor. 5.14]
    if (p eq 1) then return TraceFormulaGamma0ALTrivialNew(N, k ); end if;
    assert IsPrime(p);
    ms := [d : d in Divisors(N) | (N mod d^2 eq 0) and (d mod p ne 0)];
    trace := &+[Integers() | MoebiusMu(m)*(TraceFormulaGamma0AL(p, N div m^2, k) - TrivialContribution(N div m^2, k, p)) : m in ms];
    return trace;
end function;

function A1(n,N,k)
// Returns A1 from [Oestrle, Theorem 3]
    if (not IsSquare(n)) or (GCD(n,N) ne 1) then
        return 0;
    end if;
    return n^(k div 2 - 1)*phi1(N)*(k-1)/12;
end function;

function phi1(N)
// !!! TODO - This function repeats a previous one. Check if we can remove one of them
    return N * &*[ Rationals() | 1 + 1/p : p in PrimeDivisors(N)];
end function;

function mu(N,t,f,n)
// Returns mu(t,f,n) as described in [Oestrle, Theorem 3]
    N_f := GCD(N,f);
    primes := [x[1] : x in Factorization(N) | (N div N_f) mod x[1] ne 0];
    s := #[x : x in [0..N-1] | (GCD(x,N) eq 1) and ((x^2 - t*x+n) mod (GCD(f*N, N^2)) eq 0)];
    prod := IsEmpty(primes) select 1 else &*[ 1 + 1/p : p in  primes];
    assert N_f * prod eq (phi1(N) / phi1(N div GCD(N,f)));
    return N_f * prod * s;
end function;

function A2(n,N,k)
// Returns A2 from [Oestrle, Theorem 3]
    R<x> := PolynomialRing(Rationals());
    max_abst := Floor(SquareRoot(4*n));
    if IsSquare(n) then max_abst -:= 1; end if;
    ret := 0;
    for t in [-max_abst..max_abst] do
        // print "t = ", t;
        F<rho> := NumberField(x^2 - t*x+n);
        rho_bar := t - rho;
        p := Rationals()!((rho^(k-1) - rho_bar^(k-1)) / (rho - rho_bar));
        assert p eq P(k,t,n);
        t_sum := 0;
        for d in Divisors(4*n - t^2) do
            is_sq, f := IsSquare(d);
            if is_sq then
                D := (t^2-4*n) div f^2;
                if D mod 4 in [0,1] then
                    O := QuadraticOrder(BinaryQuadraticForms(D));
                    h := #PicardGroup(O);
                    w := #TorsionSubgroup(UnitGroup(O));
                    t_sum +:= (h/w) * mu(N,t,f,n);
                end if;
            end if;
        end for;
        // print "t_sum = ", t_sum;
        // print "p = ", p;
        ret -:= p*t_sum;
        // print "ret = ", ret;
    end for;
    return ret;
end function;

function A3(n,N,k)
// Returns A3 from [Oestrle, Theorem 3]
    g := 1;
    ds := [d : d in Divisors(n) | d^2 lt n];
    ret := 0;
    for d in ds do
        // print "d = ", d;
        cs := [c : c in Divisors(N) |
            (GCD(N div g, n div d - d) mod GCD(c, N div c) eq 0)
            and (GCD(d mod c, c) eq 1) and (GCD((n div d) mod (N div c), (N div c)) eq 1)];
        ret -:= d^(k-1) * &+[Integers() | EulerPhi(GCD(c, N div c)) : c in cs];
        // print "ret = ", ret;
    end for;
    is_sq, d := IsSquare(n);
    if is_sq then
        cs := [c : c in Divisors(N) |
            (GCD(N div g, n div d - d) mod GCD(c, N div c) eq 0)
            and (GCD(d mod c, c) eq 1) and (GCD((n div d) mod (N div c), (N div c)) eq 1)];
        ret -:= 1/2 * d^(k-1) * &+[Integers() | EulerPhi(GCD(c, N div c)) : c in cs];
    end if;
    return ret;
end function;

function A4(n,N,k)
// Returns A4 from [Oestrle, Theorem 3]
    if k eq 2 then
        return &+[t : t in Divisors(n) | GCD(N, n div t) eq 1];
    end if;
    return 0;
end function;

function TraceCohen(n,N,k)
// Returns the trace of T_n on S_k(N) using Cohen-Oesterle's formula [Oestrle, Theorem 3]
    return A1(n,N,k) + A2(n,N,k) + A3(n,N,k) + A4(n,N,k);
end function;

function get_trace(N, k, n : New := false)
// Returns the trace of T_n \circ W_N on S_k(N). If New is set to true returns the trace on S_k(N)^{new}.
// Computes the trace using modular symbols
    C := CuspidalSubspace(ModularSymbols(N,k,1));
    if New then
	    C := NewSubspace(C);
    end if;
    al := AtkinLehner(C,N);
    T := HeckeOperator(C,n);
    return Trace(T*al);
end function;

procedure testInverseRelationNewSubspacesTrivial(N, k)
// Verifies that computing the trace of W_N on S_k(N)^{new} using Mobius inversion 
// on the next formula yields the correct result
    ms := [d : d in Divisors(N) | N mod d^2 eq 0];
    trace_1 := &+[Integers() | MoebiusMu(m) * get_trace(N div m^2, k, 1) : m in ms];
    assert trace_1 eq get_trace(N, k, 1 : New);
end procedure;

procedure testRelationNewSubspacesTrivial(N, k)
// Verifies that the trace of W_N on S_k(N) is equal to the sum of traces of W_{N'} on S_k(N')^{new}
// where N_prime runs over divisors of N such that N/N' is a square
    N_primes := [d : d in Divisors(N) | IsSquare(N div d)];
    trace_1 := &+[Integers() | get_trace(N_prime, k, 1 : New) : N_prime in N_primes];
    assert trace_1 eq get_trace(N, k, 1);
end procedure;

function TrivialContribution(N, k, p)
// !! TODO - this repeats a previous function, check if we can remove one.
    N_primes_2 := [d : d in Divisors(N) | IsSquare(N*p div d) and ((N div d) mod p^3 ne 0) and (d mod p ne 0)];
    trace_2 := &+[Integers() | get_trace(N_prime, k, 1 : New) : N_prime in N_primes_2];
    trace_3 := 0;
    if (N mod p eq 0) then
        N_primes_3 := [d : d in Divisors(N div p) | IsSquare(N div (d*p))];
        trace_3 := &+[Integers() | get_trace(N_prime, k, 1 : New) :  N_prime in N_primes_3];
    end if;
    return p^(k div 2) * trace_3  - p^(k div 2 - 1)*trace_2;
end function;

procedure testRelationNewSubspaces(N, k, p)
// Verifies that the trace of T_p \circ W_N on S_k(N) satisfies [Assaf, (4.28)]
    N_primes_1 := [d : d in Divisors(N) | IsSquare(N div d) and ((N div d) mod p ne 0)];
    trace_1 := &+[ Integers() | get_trace(N_prime, k, p : New) : N_prime in N_primes_1];
    assert trace_1 + TrivialContribution(N, k, p) eq get_trace(N, k, p);
end procedure;

procedure testInverseRelationNewSubspaces(N, k, p)
// Verifies that the trace of T_p \circ W_N on S_k(N)^{new} satisfies [Assaf, Cor. 5.5]
    ms := [d : d in Divisors(N) | (N mod d^2 eq 0) and (d mod p ne 0)];
    // N_primes := [d : d in Divisors(N) | IsSquare(N div d) and ((N div d) mod p ne 0)];
    trace := &+[Integers() | MoebiusMu(m)*(get_trace(N div m^2, k, p) - TrivialContribution(N div m^2, k, p)) : m in ms];
    assert trace eq get_trace(N, k, p : New);
end procedure;

// Formula from Popa
function TraceFormulaGamma0HeckeAL(N, k, n, Q)
// Returns the trace of T_n \circ W_Q on S_k(N) using [Popa, Theorem 4]
    assert k ge 2;
    if (n eq 0) then return 0; end if; // for compatibility with q-expansions
    S1 := 0;
    Q_prime := N div Q;
    assert GCD(Q, Q_prime) eq 1;
    w := k - 2;
    max_abst := Floor(SquareRoot(4*Q*n)) div Q;
    for tQ in [-max_abst..max_abst] do
	    t := tQ*Q;
	    for u in Divisors(Q) do
	        for u_prime in Divisors(Q_prime) do
		        if ((4*n*Q-t^2) mod (u*u_prime)^2 eq 0) then
		            // print "u =", u, " u_prime = ", u_prime, "t = ", t;
		            S1 +:= P(k,t,Q*n)*H((4*Q*n-t^2) div (u*u_prime)^2 )*Cfast(Q_prime,u_prime,t,Q*n)
			                *MoebiusMu(u) / Q^(w div 2);
		            // print "S1 = ", S1;
		        end if;
	        end for;
	    end for;
    end for;
    S2 := 0;
    for d in Divisors(n*Q) do
	    a := n*Q div d;
	    if (a+d) mod Q eq 0 then
	        // print "a = ", a, "d = ", d;
	        S2 +:= Minimum(a,d)^(k-1)*Phil(N,Q,a,d) / Q^(w div 2);
	        // print "S2 = ", S2;
	    end if;
    end for;
    // print "S2 = ", S2;
    ret := -S1 / 2 - S2 / 2;
    if k eq 2 then
	    ret +:= &+[n div d : d in Divisors(n) | GCD(d,N) eq 1];
    end if;
    return ret;
end function;

function get_trace_hecke_AL(N, k, n, Q : New := false)
// Returns the trace of T_n \circ W_Q on S_k(N). If New is set to true, returns instead the trace on S_k(N)^{new}.
// Computes it using Modular Symbols
    C := CuspidalSubspace(ModularSymbols(N,k,1));
    if New then
	    C := NewSubspace(C);
    end if;
    al := AtkinLehner(C,Q);
    T := HeckeOperator(C,n);
    return Trace(T*al);
end function;

function d_prime(d, N, Q, N_prime)
// Returns d' as defined in [Assaf, Lemma 3.1]
    return GCD(d, N div Q) * GCD(N div (N_prime*d), Q);
end function;

function Q_prime(N, Q, N_prime)
// Returns Q' as defined in [Assaf, Lemma 3.1]
    return GCD(Q, N_prime);
end function;

function dd_prime(d, N, Q, N_prime, n)
    d_p := d_prime(d, N, Q, N_prime);
    if (GCD(d_p, n) * d) mod d_p ne 0 then
	    return 0;
    end if;
    return (GCD(d_p, n) * d) div d_p;
end function;

function n_prime(d, N, Q, N_prime, n)
    d_p := d_prime(d, N, Q, N_prime);
    dd_p := dd_prime(d, N, Q, N_prime, n) * GCD(d_p, n);
    return n div dd_p;
end function;

function get_ds(N, Q, N_prime, n)
// Returns the list of divisors d of n appearing in the sum in [Assaf, Cor. 4.27]
    divs := Divisors(N div N_prime);
    ret := [];
    for d in divs do
        d_p := d_prime(d, N, Q, N_prime);
        dd_p := dd_prime(d, N, Q, N_prime, n);
        if (dd_p eq 0) then
            continue;
        end if;
        if (GCD(dd_p, N_prime) eq 1) and (n mod (dd_p * GCD(d_p, n)) eq 0) then
            Append(~ret, d);
        end if;
    end for;
    return ret;
end function;

procedure testRelationNewSubspacesGeneral(N, k, n, Q)
// Verifies the formula for trace of T_n \circ W_Q on S_k(N) in [Assaf, Cor. 4.27]
    s := 0;
    for N_prime in Divisors(N) do
        ds := get_ds(N, Q, N_prime, n);
        for d in ds do
            n_p := n_prime(d, N, Q, N_prime, n);
            d_p := d_prime(d, N, Q, N_prime);
            dd_p := dd_prime(d, N, Q, N_prime, n);
            Q_p := Q_prime(N, Q, N_prime);
            term := (n div n_p)^(k div 2 - 1);
            term *:= GCD(d_p, n);
            term *:= MoebiusMu(dd_p);
            term *:= get_trace_hecke_AL(N_prime, k, n_p, Q_p : New);
            s +:= term;
        end for;
    end for;
    assert s eq get_trace_hecke_AL(N, k, n, Q);
end procedure;

procedure testBatchRelationNewSubspacesGeneral(Ns, ns, ks)
// Run multiple tests of the previsou verification, where
// Ns is a list of values for N, ns is a list of values for n, and ks is a list of values for k
// For each triple (N,k,n), runs the above verification for all Q || N
    // printf "(N,k,n,Q)=";
    printf "(N,n,Q)=";
    for N in Ns do
        Qs := [Q : Q in Divisors(N) | GCD(Q, N div Q) eq 1];
        for Q in Qs do
            for n in ns do
                printf "(%o,%o,%o),", N, n, Q;
                for k in ks do
                    // printf "(%o,%o,%o,%o),", N, k, n, Q;
                    testRelationNewSubspacesGeneral(N,k,n,Q);
                end for;
            end for;
        end for;
    end for;
end procedure;

function IsRelevantNprime(N_p, N, a, n_p, n)
// Returns true if N_prime appears in the sum in [Assaf, Cor. 5.5], false otherwise.
    if (GCD(N div N_p, n*a) ne n_p) then
	    return false;
    end if;
    if (((N div N_p) mod a) ne 0) then
	    return false;
    end if;
    if (GCD(N_p, a) ne 1) then
	    return false;
    end if;
    return IsSquare(N div (N_p * n_p));
end function;

procedure testRelationNewSubspaces2(N, k, n)
// Verifies that [Assaf, Cor. 5.5] holds
    s := 0;
    for n_p in Divisors(n) do
        for a in Divisors(n_p) do
            N_primes := [N_p : N_p in Divisors(N) | IsRelevantNprime(N_p, N, a, n_p, n)];
            term := &+[Integers() | get_trace_hecke_AL(N_p, k, n div n_p, N_p : New) : N_p in N_primes];
            term *:= n_p^(k div 2) div a;
            term *:= MoebiusMu(a);
            s +:= term;
        end for;
    end for;
    printf "\n";
    assert s eq get_trace_hecke_AL(N, k, n, N);
end procedure;

procedure testBatchRelationNewSubspaces2(Ns, ns, ks)
// Run multiple tests of the previsou verification, where
// Ns is a list of values for N, ns is a list of values for n, and ks is a list of values for k
// For each triple (N,k,n), runs the above verification for all Q || N
    printf "(N,n)=";
    for N in Ns do
        for n in ns do
            printf "(%o,%o),", N, n;
            for k in ks do
                testRelationNewSubspaces2(N,k,n);
            end for;
        end for;
    end for;
    printf "\n";
end procedure;

procedure testRelationNewSubspaces3(N, k, n, Q)
    s := 0;
    n_Q := GCD(n, Q);
    n_NQ := n div n_Q;
    n_p_NQs := [n_p_NQ : n_p_NQ in Divisors(n_NQ) | IsSquare(n_p_NQ)];
    for n_p_Q in Divisors(n_Q) do
        for a_Q in Divisors(n_p_Q) do
            for n_p_NQ in n_p_NQs do
                n_p := n_p_Q * n_p_NQ;
                _, a_NQ := IsSquare(n_p_NQ);
                // print a_NQ;
                a := a_Q * a_NQ;
                Q_primes := [Q_p : Q_p in Divisors(Q) | IsRelevantNprime(Q_p, Q, a_Q, n_p_Q, n_Q)];
                NQ_primes := [NQ_p : NQ_p in Divisors(N div Q) | (GCD(a_NQ, NQ_p) eq 1) and ((N div (Q*NQ_p)) mod a_NQ eq 0)];
                N_primes := [Q_p * NQ_p : Q_p in Q_primes, NQ_p in NQ_primes];
                // print N_primes;
                weights := [#[x : x in Divisors(GCD(N div Q, N div N_p)) | GCD(x,n) eq 1] : N_p in N_primes];
                // print weights;
                weights2 := [#[x : x in Divisors(GCD(N div Q, N div N_p)) | GCD(x,n) eq a_NQ] : N_p in N_primes];
                // print weights;
                // These should be the same
                assert weights eq weights2;
                traces := [get_trace_hecke_AL(N_p, k, n div n_p, GCD(N_p, Q) : New) : N_p in N_primes];
                // print traces;
                term := &+[Integers() | weights[i]*traces[i] : i in [1..#N_primes]];
                term *:= n_p^(k div 2) div a;
                term *:= MoebiusMu(a);
                s +:= term;
            end for;
        end for;
    end for;
    assert s eq get_trace_hecke_AL(N, k, n, Q);
end procedure;

procedure testBatchRelationNewSubspaces3(Ns, ns, ks)
    printf "(N,n,Q)=";
    for N in Ns do
        Qs := [Q : Q in Divisors(N) | GCD(Q, N div Q) eq 1];
        for Q in Qs do
            for n in ns do
                printf "(%o,%o,%o),", N, n, Q;
                for k in ks do
                    testRelationNewSubspaces3(N,k,n,Q);
                end for;
            end for;
        end for;
    end for;
    printf "\n";
end procedure;

function alpha(Q, n, m)
// Returns the arithmetic function alpha_{Q,n}(m) defined in [Assaf, Lemma 5.3]
    if m eq 1 then return 1; end if;
    fac := Factorization(m);
    if #fac gt 1 then
	    return &*[alpha(Q,n,pe[1]^pe[2]) : pe in fac];
    end if;
    p := fac[1][1];
    e := fac[1][2];
    if (e eq 2) and ((Q*n mod p) ne 0) then
	    return 1;
    end if;
    if (e eq 1) and ((Q*n mod p) ne 0) then
	    return -2;
    end if;
    if (e eq 2) and ((Q mod p) eq 0) and ((n mod p) ne 0) then
	    return -1;
    end if;
    if (e eq 1) and (n mod p eq 0) and (Q mod p ne 0) then
	    return -1;
    end if;
    return 0;
end function;

forward TraceFormulaGamma0HeckeALNew;

function TraceFormulaGamma0HeckeALNewSmaller(N, k, n, Q)
// Returns T_{<n,k}(Q,N) as defined in [Assaf, Cor. 4.27]
    trace := 0;
    n_Q := GCD(n, Q);
    n_NQ := n div n_Q;
    n_p_NQs := [n_p_NQ : n_p_NQ in Divisors(n_NQ) | IsSquare(n_p_NQ)];
    for n_p_Q in Divisors(n_Q) do
        for d_Q in Divisors(n_p_Q) do
            for n_p_NQ in n_p_NQs do
                n_p := n_p_Q * n_p_NQ;
                if (n_p eq 1) then
                    continue;
                end if;
                _, d_NQ := IsSquare(n_p_NQ);
                d := d_Q * d_NQ;
                Q_primes := [Q_p : Q_p in Divisors(Q) | IsRelevantNprime(Q_p, Q, d_Q, n_p_Q, n_Q)];
                NQ_primes := [NQ_p : NQ_p in Divisors(N div Q) | (GCD(d_NQ, NQ_p) eq 1) and ((N div (Q*NQ_p)) mod d_NQ eq 0)];
                N_primes := [Q_p * NQ_p : Q_p in Q_primes, NQ_p in NQ_primes];
                weights := [#[x : x in Divisors(GCD(N div Q, N div N_p)) | GCD(x,n) eq 1] : N_p in N_primes];
                traces := [TraceFormulaGamma0HeckeALNew(N_p, k, n div n_p, GCD(N_p, Q) ) : N_p in N_primes];
                term := &+[Integers() | weights[i]*traces[i] : i in [1..#N_primes]];
                term *:= n_p^(k div 2) div d;
                term *:= MoebiusMu(d);
                trace +:= term;
            end for;
        end for;
    end for;
    return trace;
end function;

function TraceFormulaGamma0HeckeALNew(N, k, n, Q)
// Returns the trace of T_n \circ W_Q on S_k(N)^{new} using [Assaf, Cor. 5.5]
    trace := 0;
    for N_prime in Divisors(N) do
	    a := alpha(Q, n, N div N_prime);
	    Q_prime := GCD(N_prime, Q);
	    term := TraceFormulaGamma0HeckeAL(N_prime, k, n, Q_prime);
	    term -:= TraceFormulaGamma0HeckeALNewSmaller(N_prime, k, n, Q_prime );
	    trace +:= a*term;
    end for;
    return trace;
end function;

procedure testBatchTraceFormulaGamma0HeckeALNew(Ns, ns, ks)
// Verifies [Assaf, Cor. 5.5] for all triples (N,n,k) 
    printf "(N,n,Q)=";
    for N in Ns do
        Qs := [Q : Q in Divisors(N) | GCD(Q, N div Q) eq 1];
        for Q in Qs do
            for n in ns do
                printf "(%o,%o,%o),", N, n, Q;
                for k in ks do
                    assert TraceFormulaGamma0HeckeALNew(N,k,n,Q ) eq get_trace_hecke_AL(N,k,n,Q : New);
                end for;
            end for;
        end for;
    end for;
    printf "\n";
end procedure;

// [Assaf] - E. Assaf, a note on the trace formula
// [Oesterle] - J. Oesterle - Sur la Trace des Operateurs de Hecke (Thesis)
// [Popa] - A. Popa, On the trace formula for Hecke operators on congruence subgroups, II
// [SZ88] - N.P. Skoruppa and D. Zagier, Jacobi forms and a certain space of modular forms