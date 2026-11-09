// X_0^4130(1)* IS HYPERELLIPTIC, PART TWO: A FIBRE DIFFERENCE IS TORSION.  Many hours
// (the cusp forms of weight 3/2 on Gamma_0(8260) to precision 4000).
//
// On Y = X_0^4130(1)* (genus three, W the full Atkin-Lehner group) let Z(n) = sum_r (1/e_d) P_d over the discriminants
// d = -4n/r^2 with CM points (Guo-Yang Lemma 25).  For a linear form phi on Pic(Y)_Q vanishing on the Hodge class,
// sum_n phi(Z(n)) q^n is a cusp form of weight 3/2 on Gamma_0(8260) with the theta^3 multiplier vanishing at every index
// carrying no CM point (paper/hyperelliptic-15-38.tex, Section 2); the space W of such forms carries T(p^2) with the
// characteristic polynomial of the traces of Frobenius of Y at p (checked at p = 3 against the point counts).  A
// degree-zero CM divisor is torsion as soon as some principal part with that divisor pairs to zero with W.  The rational
// torsion of J(Y) is trivial (tests/Hyperelliptic4130_1.m), so torsion is principal, and ONE certified difference of two
// distinct degree-two divisors is a g^1_2.  This test certifies the six differences of the four fibres
// D_1 = P_-35 + P_-532, D_2 = P_-163 + P_-595, D_3 = P_-315 + P_-520, D_4 = P_-340 + P_-420 found by the Shimizu screen
// (any one suffices) and refuses two non-relations.  P_-35 is fixed by w_35 (e = 2) and -315 = 9 * (-35) is not
// fundamental; both are carried by the dictionary (paper/hyperelliptic-4130-1.tex).

procedure test_hyperelliptic_4130_1_torsion()
    printf "X_0^4130(1)*: the differences of the four fibres P_-35 + P_-532, P_-163 + P_-595, P_-315 + P_-520, P_-340 + P_-420 are torsion (level 8260, the space W)...";
    D := 4130; N := 1; M := 8260; B := 4000;
    function HasCM(d)
        d0 := FundamentalDiscriminant(d); _, f := IsSquare(d div d0);
        if exists{p : p in PrimeDivisors(D) | f mod p eq 0 or KroneckerSymbol(d0, p) eq 1} then return false; end if;
        if exists{p : p in PrimeDivisors(N) | f mod p ne 0 and KroneckerSymbol(d0, p) eq -1} then return false; end if;
        return true;
    end function;
    function Stab(d)
        e := d eq -3 select 3 else (d eq -4 select 2 else 1);
        if exists{m : m in Divisors(D*N) | m gt 1 and GCD(m, (D*N) div m) eq 1 and (d eq -4*m or (m mod 4 eq 3 and d eq -m))} then e *:= 2; end if;
        return e;
    end function;
    discs := func< n | [ -4*n div (r^2) : r in [1..Isqrt(4*n)] | (4*n) mod (r^2) eq 0 and ((-4*n div (r^2)) mod 4) in {0, 1} ] >;
    S := CuspidalSubspace(HalfIntegralWeightForms(Gamma0(M), 3/2));
    Bq := Basis(S, B + 2);
    A0 := Matrix(Rationals(), #Bq, B, [ Coefficient(g, n) : n in [1..B], g in Bq ]);
    cmd := AssociativeArray();
    for n in [1..B] do cmd[n] := [ d : d in discs(n) | HasCM(d) ]; end for;
    empty := [ n : n in [1..B] | IsEmpty(cmd[n]) ];
    K := Kernel(Matrix(Rationals(), #Bq, #empty, [ A0[i, n] : n in empty, i in [1..#Bq] ]));
    Wb := [ Vector(k) * A0 : k in Basis(K) ];
    assert #Wb ge 1;                       // every form of W is used, which is the sufficient direction
    // T(p^2) on W: a(p^2 n) + (-n/p) a(n) + p a(n/p^2), read on the indices up to B/p^2; its characteristic polynomial
    // divides prod (x - a_p(f)) over the three newforms of Y, which is x^3 + c1 x^2 + (c2 - 3p) x + (c3 - 2p c1) for
    // L_p(x) = 1 + c1 x + c2 x^2 + c3 x^3 + ...
    W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
    Y := CreateShimuraQuot(D, N, W); Y`g := GenusShimuraCurveQuotient(D, N, W); Y`CurveID := 0;
    p := 3; Bp := B div p^2;
    s := [ p^i + 1 - ComputePointsViaTrace(Y, p, i) : i in [1..3] ];
    e := [Rationals() | 1];
    for k in [1..3] do Append(~e, (&+[ (-1)^(i-1) * e[k-i+1] * s[i] : i in [1..k] ]) / k); end for;
    c1 := -e[2]; c2 := e[3]; c3 := -e[4];
    Tp := func< g | [ g[p^2*n] + KroneckerSymbol(-n, p) * g[n] + (n mod p^2 eq 0 select p * g[n div p^2] else 0) : n in [1..Bp] ] >;
    Wmat := Matrix(Rationals(), [ [ g[n] : n in [1..Bp] ] : g in Wb ]);
    ok, Tm := IsConsistent(Wmat, Matrix(Rationals(), [ Tp(g) : g in Wb ]));
    determined := ok and Rank(Wmat) eq #Wb;      // T(p^2) on W is read off the indices up to B/p^2, which must separate W
    R<x> := PolynomialRing(Rationals());
    if determined then assert IsDivisibleBy(x^3 + c1 * x^2 + (c2 - 3*p) * x + (c3 - 2*p*c1), CharacteristicPolynomial(Tm)); end if;
    cmds := Sort(Setseq(&join{ Set(cmd[n]) : n in [1..B] }));
    function Certified(target)
        tg := AssociativeArray(); for t in target do tg[t[1]] := t[2]; end for;
        assert forall{t : t in target | t[1] in cmds};
        nrows := #Wb + #cmds;
        A := ZeroMatrix(Rationals(), nrows, B); rhs := ZeroMatrix(Rationals(), nrows, 1);
        for i->g in Wb do for n in [1..B] do A[i, n] := g[n]; end for; end for;
        for j->d in cmds do
            for n in [1..B] do if d in cmd[n] then A[#Wb + j, n] +:= 1; end if; end for;
            rhs[#Wb + j, 1] := IsDefined(tg, d) select Stab(d) * tg[d] else 0;
        end for;
        return IsConsistent(Transpose(A), Transpose(rhs));
    end function;
    fibres := [ [<-35, 1>, <-532, 1>], [<-163, 1>, <-595, 1>], [<-315, 1>, <-520, 1>], [<-340, 1>, <-420, 1>] ];
    for i in [1..4], j in [i+1..4] do
        assert Certified(fibres[i] cat [ <t[1], -t[2]> : t in fibres[j] ]);
    end for;
    assert not Certified([<-35, 1>, <-532, -1>]);
    assert not Certified([<-35, 1>, <-163, 1>, <-532, -1>, <-595, -1>]);
    printf " [dim W = %o, T(9) on it %o; six differences certified, two controls refused] Done!\n", #Wb, determined select "divides the curve's" else "not determined by the coefficients available";
end procedure;

test_hyperelliptic_4130_1_torsion();
