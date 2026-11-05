// X_0^39(34)* IS HYPERELLIPTIC, PART TWO: FIVE CM DIVISORS OF DEGREE TWO ARE PAIRWISE TORSION-RELATED.
// About thirty-five minutes (the weight-3/2 cusp forms of level 2652 to precision 4000).
//
// Setting (paper/hyperelliptic-21-58.tex, Section 3, and paper/hyperelliptic-39-34.tex).  With
// X* = X_0^39(34)* and X'* = X_0^39(17)*, both at scalar level M = 2652, let W+ be the plus-space cusp
// forms of level 2652 vanishing at every index n where H(n) = sum_s P'_{-n/s^2} is the zero divisor on
// X'*.  Every component of the plus-type Heegner series H of X'* lies in W+, and every Heegner
// component of X* is (U(4) + 2) applied to one.  Two forms of the test follow, for a degree-zero CM
// divisor Delta on X*:
//   (i)  Delta = sum c(-n) Z(n) in the pole-at-infinity dictionary of X* (even discriminants), and c
//        orthogonal to V = (U(4) + 2) W+;
//   (ii) T Delta = sum c'(n) H(n), T the correspondence to X'* given on CM divisors by the fibre count
//        T C(d) = (1 + chi_d(2)) C'(d) + (2 - chi_{d/4}(2)) C'(d/4) (C = P/e the orbifold divisors), and c'
//        orthogonal to W+; T is injective on the 2-old Jacobian, and in this dictionary odd
//        discriminants can be prescribed on their own.
// Either certifies Delta torsion.  The divisors: F = P_-84 + P_-132, F1 = P_-15 + P_-52,
// F7 = P_-1092 + P_-372, Q1 = Z_-3432, Q2 = Z_-1560; the test certifies F7 - F, Q1 - F, Q2 - F by (i) and
// F1 - F by (ii), and refuses five differences of points that are not fibres.  With the torsion bound 4
// of tests/Hyperelliptic39_34.m, two of the five classes coincide and the curve carries a g^1_2.

procedure test_hyperelliptic_39_34_relations()
    printf "X_0^39(34)*: five CM divisors of degree two are pairwise torsion-related (level 2652, plus space)...";
    D := 39; N := 34; Np := 17; M := 2652; B := 4000;
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
    function Chi2(d)
        d0 := FundamentalDiscriminant(d); _, f := IsSquare(d div d0);
        return f mod 2 eq 0 select 0 else KroneckerSymbol(d0, 2);
    end function;
    discs := func< n | [ -4*n div (r^2) : r in [1..Isqrt(4*n)] | (4*n) mod (r^2) eq 0 and ((-4*n div (r^2)) mod 4) in {0, 1} ] >;
    pdiscs := func< n | [ -n div (s^2) : s in [1..Isqrt(n)] | n mod (s^2) eq 0 and ((-n div (s^2)) mod 4) in {0, 1} ] >;
    S := CuspidalSubspace(HalfIntegralWeightForms(Gamma0(M), 3/2));
    assert Dimension(S) eq 245;
    Bq := Basis(S, B + 2);
    A0 := Matrix(Rationals(), #Bq, B, [ Coefficient(g, n) : n in [1..B], g in Bq ]);
    Hzero := [ n : n in [1..B] | (n mod 4) in {1, 2} or forall{d : d in pdiscs(n) | not HasCM(d, D, Np)} ];
    K := Kernel(Matrix(Rationals(), #Bq, #Hzero, [ A0[i, n] : n in Hzero, i in [1..#Bq] ]));
    Wp := [ Vector(k) * A0 : k in Basis(K) ];
    assert #Wp eq 3;
    // (i) the pole-at-infinity dictionary of X*, against V = (U(4) + 2) W+
    B4 := B div 4;
    V := [ [ g[4*n] + 2*g[n] : n in [1..B4] ] : g in Wp ];
    cmd := AssociativeArray();
    for n in [1..B4] do cmd[n] := [ d : d in discs(n) | HasCM(d, D, N) ]; end for;
    full := [ n : n in [1..B4] | not IsEmpty(cmd[n]) ];
    cmds := Sort(Setseq(&join{ Set(cmd[n]) : n in full }));
    function CertifiedInf(target)
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
    // (ii) the plus-type dictionary of X'*, against W+
    plusidx := [ n : n in [1..B] | n notin Hzero ];
    hd := AssociativeArray();
    for n in plusidx do hd[n] := [ d : d in pdiscs(n) | HasCM(d, D, Np) ]; end for;
    dprime := Sort(Setseq(&join{ Set(hd[n]) : n in plusidx }));
    function CertifiedPlus(target)
        mu := AssociativeArray();
        for t in target do
            d := t[1]; alpha := Stab(d, D, N) * t[2];
            assert HasCM(d, D, N) and HasCM(d, D, Np);
            mu[d] := (IsDefined(mu, d) select mu[d] else 0) + alpha * (1 + Chi2(d));
            if d mod 4 eq 0 and ((d div 4) mod 4) in {0, 1} and HasCM(d div 4, D, Np) then
                mu[d div 4] := (IsDefined(mu, d div 4) select mu[d div 4] else 0) + alpha * (2 - Chi2(d div 4));
            end if;
        end for;
        assert forall{d : d in Keys(mu) | d in dprime};
        nrows := #Wp + #dprime;
        A := ZeroMatrix(Rationals(), nrows, #plusidx); rhs := ZeroMatrix(Rationals(), nrows, 1);
        for i->w in Wp do for j->n in plusidx do A[i, j] := w[n]; end for; end for;
        for k->d in dprime do
            for j->n in plusidx do if d in hd[n] then A[#Wp + k, j] +:= 1; end if; end for;
            rhs[#Wp + k, 1] := IsDefined(mu, d) select mu[d] else 0;
        end for;
        return IsConsistent(Transpose(A), Transpose(rhs));
    end function;
    F := [<-84, 1>, <-132, 1>]; F1 := [<-15, 1>, <-52, 1>]; F7 := [<-1092, 1>, <-372, 1>]; Q1 := [<-3432, 1>]; Q2 := [<-1560, 1>];
    dif := func< X, Y | X cat [ <t[1], -t[2]> : t in Y ] >;
    assert CertifiedInf(dif(F7, F)) and CertifiedInf(dif(Q1, F)) and CertifiedInf(dif(Q2, F));
    assert CertifiedPlus(dif(F1, F));
    // the two dictionaries agree where both apply, and neither is vacuous
    assert CertifiedPlus(dif(F7, F));
    assert not CertifiedInf([<-84, 1>, <-132, -1>]) and not CertifiedInf([<-1092, 1>, <-372, -1>]) and not CertifiedInf([<-84, 1>, <-372, -1>]);
    assert not CertifiedPlus([<-15, 1>, <-52, -1>]) and not CertifiedPlus([<-15, 1>, <-84, -1>]);
    printf " [dim W+ = 3; four relations certified, five controls refused] Done!\n";
end procedure;

test_hyperelliptic_39_34_relations();
