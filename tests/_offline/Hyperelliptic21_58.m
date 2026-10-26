// X_0^21(58)* IS HYPERELLIPTIC, PART TWO: FOUR CM DIVISORS OF DEGREE TWO ARE PAIRWISE TORSION-RELATED.
// About half an hour (the weight-3/2 cusp forms of level 2436 to precision 1600).
//
// Setting (paper/hyperelliptic-21-58.tex).  X* = X_0^21(58)*, X'* = X_0^21(29)*, both with Borcherds
// forms coming from weight-1/2 forms of level M = 2436.  For a degree-zero divisor Delta on X*
// supported on CM points, write Delta = sum_n c(-n) Z(n) with Z(n) = sum_r (1/e_d) P_d over the
// discriminants d = -4n/r^2 carrying CM points of X*.  The lemma of the note says: every Heegner
// component of X* is (U(4) + 2) applied to a form of the Kohnen plus space of level 2436 that
// vanishes at every index n for which H(n) = sum_s P'_{-n/s^2} is the zero divisor on X'*.  So if
// c is orthogonal to V = (U(4) + 2) W+, with
//     W+ = { g in S_{3/2}(Gamma_0(2436)) : a_g(n) = 0 unless n = 0, 3 mod 4 and some -n/s^2 carries CM on X'* },
// then Delta is torsion in Jac(X*)(Q).  This test computes W+ (its dimension is 3, the rank of the
// Jacobian), and checks that the differences of the four fibres
//     F2 = P_-16 + P_-120,  F4 = P_-88 + P_-100,  F5 = P_-112 + P_-372,  F6 = P_-168 + P_-228
// are certified torsion, while five differences of points that are not fibres are refused (so the
// test is not vacuous).  Part one (tests/Hyperelliptic21_58.m) checks that the eight points are
// rational, that the Jacobian is the one the lemma needs, and that its torsion has order at most 2:
// four effective divisors of degree two in one torsion coset of a group of order at most 2 means two
// of them are linearly equivalent, and a genus-three curve with two distinct effective divisors of
// degree two in one class carries a g^1_2.

procedure test_hyperelliptic_21_58_relations()
    printf "X_0^21(58)*: the four even fibres are pairwise torsion-related (level 2436, plus space)...";
    D := 21; N := 58; Np := 29; M := 2436; B := 1600;
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
    assert Dimension(S) eq 233;
    Bq := Basis(S, B + 2);
    A0 := Matrix(Rationals(), #Bq, B, [ Coefficient(g, n) : n in [1..B], g in Bq ]);
    // W+: vanishing where H vanishes on X'*
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
    F2 := [<-16, 1>, <-120, 1>]; F4 := [<-88, 1>, <-100, 1>]; F5 := [<-112, 1>, <-372, 1>]; F6 := [<-168, 1>, <-228, 1>];
    dif := func< X, Y | X cat [ <t[1], -t[2]> : t in Y ] >;
    assert Certified(dif(F2, F4)) and Certified(dif(F2, F5)) and Certified(dif(F2, F6)) and Certified(dif(F4, F5));
    // not vacuous: differences of two points that do not form a fibre are refused
    assert not Certified([<-16, 1>, <-120, -1>]) and not Certified([<-88, 1>, <-100, -1>]) and not Certified([<-16, 1>, <-88, -1>])
       and not Certified([<-112, 1>, <-372, -1>]) and not Certified([<-16, 1>, <-168, -1>]);
    printf " [dim W+ = 3; four fibre differences certified, five controls refused] Done!\n";
end procedure;

test_hyperelliptic_21_58_relations();
