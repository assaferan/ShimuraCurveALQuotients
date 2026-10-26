// X_0^77(10)* IS HYPERELLIPTIC, PART TWO: FIVE CM DIVISORS OF DEGREE TWO ARE PAIRWISE TORSION-RELATED.
// About twenty-five minutes (the weight-3/2 cusp forms of level 1540 to precision 8000).
//
// Setting (paper/hyperelliptic-21-58.tex, Section 3, and paper/hyperelliptic-77-10.tex).  With
// X* = X_0^77(10)* and X'* = X_0^77(5)*, both at scalar level M = 1540, every Heegner component of X*
// is (U(4) + 2) applied to a cusp form of the Kohnen plus space of level 1540 that vanishes at every
// index n for which H(n) = sum_s P'_{-n/s^2} is the zero divisor on X'*.  So a degree-zero divisor
// Delta = sum_n c(-n) Z(n) on X* with c orthogonal to V = (U(4) + 2) W+,
//     W+ = { g in S_{3/2}(Gamma_0(1540)) : a_g(n) = 0 unless n = 0, 3 mod 4 and some -n/s^2 carries CM on X'* },
// is torsion in Jac(X*)(Q).  This test computes W+ (dimension 3, the rank) and certifies that the
// differences of the five divisors
//     F2 = P_-16 + P_-100,  F3 = P_-280 + P_-660,  F7 = P_-760 + P_-1540,  F8 = P_-840 + P_-1320,  Q = Z_-4180
// with F2 are torsion, while five differences of points that are not fibres are refused.  With the
// torsion bound 2 of tests/Hyperelliptic77_10.m, two of the five divisors are linearly equivalent
// and the genus-three curve carries a g^1_2.

procedure test_hyperelliptic_77_10_relations()
    printf "X_0^77(10)*: five CM divisors of degree two are pairwise torsion-related (level 1540, plus space)...";
    D := 77; N := 10; Np := 5; M := 1540; B := 8000;
    function HasCM(d, D, N)
        d0 := FundamentalDiscriminant(d); _, f := IsSquare(d div d0);
        if exists{p : p in PrimeDivisors(D) | f mod p eq 0 or KroneckerSymbol(d0, p) eq 1} then return false; end if;
        if exists{p : p in PrimeDivisors(N) | f mod p ne 0 and KroneckerSymbol(d0, p) eq -1} then return false; end if;
        return true;
    end function;
    function Stab(d, D, N)
        e := d eq -3 select 3 else (d eq -4 select 2 else 1);
        if exists{m : m in Divisors(D*N) | m gt 1 and GCD(m, (D*N) div m) eq 1 and (d eq -4*m or (m mod 4 eq 3 and d eq -m))} then e *:= 2; end if;
        return e;
    end function;
    discs := func< n | [ -4*n div (r^2) : r in [1..Isqrt(4*n)] | (4*n) mod (r^2) eq 0 and ((-4*n div (r^2)) mod 4) in {0, 1} ] >;
    pdiscs := func< n | [ -n div (s^2) : s in [1..Isqrt(n)] | n mod (s^2) eq 0 and ((-n div (s^2)) mod 4) in {0, 1} ] >;
    S := CuspidalSubspace(HalfIntegralWeightForms(Gamma0(M), 3/2));
    assert Dimension(S) eq 137;
    Bq := Basis(S, B + 2);
    A0 := Matrix(Rationals(), #Bq, B, [ Coefficient(g, n) : n in [1..B], g in Bq ]);
    Hzero := [ n : n in [1..B] | (n mod 4) in {1, 2} or forall{d : d in pdiscs(n) | not HasCM(d, D, Np)} ];
    K := Kernel(Matrix(Rationals(), #Bq, #Hzero, [ A0[i, n] : n in Hzero, i in [1..#Bq] ]));
    Wp := [ Vector(k) * A0 : k in Basis(K) ];
    assert #Wp eq 3;
    B4 := B div 4;
    V := [ [ g[4*n] + 2*g[n] : n in [1..B4] ] : g in Wp ];
    cmd := AssociativeArray();
    for n in [1..B4] do cmd[n] := [ d : d in discs(n) | HasCM(d, D, N) ]; end for;
    full := [ n : n in [1..B4] | not IsEmpty(cmd[n]) ];
    cmds := Sort(Setseq(&join{ Set(cmd[n]) : n in full }));
    function Certified(target)
        tg := AssociativeArray(); for t in target do tg[t[1]] := t[2]; end for;
        assert forall{t : t in target | t[1] in cmds};
        nrows := #V + #cmds;
        A := ZeroMatrix(Rationals(), nrows, B4); rhs := ZeroMatrix(Rationals(), nrows, 1);
        for i->v in V do for n in [1..B4] do A[i, n] := v[n]; end for; end for;
        for j->d in cmds do
            for n in full do if d in cmd[n] then A[#V + j, n] +:= 1; end if; end for;
            rhs[#V + j, 1] := IsDefined(tg, d) select Stab(d, D, N) * tg[d] else 0;
        end for;
        return IsConsistent(Transpose(A), Transpose(rhs));
    end function;
    F2 := [<-16, 1>, <-100, 1>]; F3 := [<-280, 1>, <-660, 1>]; F7 := [<-760, 1>, <-1540, 1>]; F8 := [<-840, 1>, <-1320, 1>];
    Q := [<-4180, 1>];
    dif := func< X, Y | X cat [ <t[1], -t[2]> : t in Y ] >;
    assert Certified(dif(F2, F3)) and Certified(dif(F2, F7)) and Certified(dif(F2, F8)) and Certified(dif(F2, Q));
    // not vacuous: differences of two points that do not form a fibre are refused
    assert not Certified([<-16, 1>, <-100, -1>]) and not Certified([<-280, 1>, <-660, -1>]) and not Certified([<-16, 1>, <-280, -1>])
       and not Certified([<-760, 1>, <-1540, -1>]) and not Certified([<-840, 1>, <-1320, -1>]);
    printf " [dim W+ = 3; four relations certified, five controls refused] Done!\n";
end procedure;

test_hyperelliptic_77_10_relations();
