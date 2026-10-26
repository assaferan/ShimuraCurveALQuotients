// X_0^15(38)/<w_2, w_15, w_19> IS HYPERELLIPTIC, PART TWO: THE BRANCH PAIRING IS TORSION.  About ten minutes
// (the 113 cusp forms of weight 3/2 on Gamma_0(1140) to precision 4000).
//
// On the genus-one star E = X_0^15(38)/W let Z(n) = sum_r (1/e_d) P_d over the discriminants d = -4n/r^2 with
// CM points (the divisor of the Borcherds lift of q^-n, Guo-Yang Lemma 25).  For a linear form phi on
// Pic(E)_Q vanishing on the Hodge class, Phi_phi = sum_n phi(Z(n)) q^n is a cusp form of weight 3/2 on
// Gamma_0(1140) with the theta^3 multiplier (the e_0-component of Borcherds' generating series), and it
// vanishes at every index n carrying no CM point.  The space W of such forms is three-dimensional here; the
// Hecke operator T(49) acts on it with eigenvalues -4, -2, 4, and Phi_phi lies in the a_7 = -4 eigenspace
// W_f, which is one-dimensional, since a_7(E) = -4 (570e4).  A degree-zero CM divisor Delta is torsion as
// soon as some principal part c with divisor Delta pairs to zero with W_f (paper/hyperelliptic-15-38.tex).
// This test certifies Delta = P_-228 + P_-760 - P_-12 - P_-40 and refuses the other two pairings of the four
// branch points, the differences within the pairs, and the difference P_-532 - P_-40 of two points that
// reduce to the same component at 3 and at 5.

procedure test_hyperelliptic_15_38_torsion()
    printf "X_0^15(38)*: P_-228 + P_-760 - P_-12 - P_-40 is torsion (level 1140, the a_7 = -4 eigenspace)...";
    D := 15; N := 38; M := 1140; B := 4000;
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
    assert Dimension(S) eq 113;
    Bq := Basis(S, B + 2);
    A0 := Matrix(Rationals(), #Bq, B, [ Coefficient(g, n) : n in [1..B], g in Bq ]);
    cmd := AssociativeArray();
    for n in [1..B] do cmd[n] := [ d : d in discs(n) | HasCM(d) ]; end for;
    empty := [ n : n in [1..B] | IsEmpty(cmd[n]) ];
    K := Kernel(Matrix(Rationals(), #Bq, #empty, [ A0[i, n] : n in empty, i in [1..#Bq] ]));
    Wb := [ Vector(k) * A0 : k in Basis(K) ];
    assert #Wb eq 3;
    // T(49) on W: a(49 n) + (-n/7) a(n) + 7 a(n/49), read on the indices up to B/49
    p := 7; Bp := B div p^2;
    Tp := func< g | [ g[p^2*n] + KroneckerSymbol(-n, p) * g[n] + (n mod p^2 eq 0 select p * g[n div p^2] else 0) : n in [1..Bp] ] >;
    Wmat := Matrix(Rationals(), [ [ g[n] : n in [1..Bp] ] : g in Wb ]);
    ok, Tm := IsConsistent(Wmat, Matrix(Rationals(), [ Tp(g) : g in Wb ]));
    assert ok;
    R<x> := PolynomialRing(Rationals());
    assert CharacteristicPolynomial(Tm) eq (x + 4) * (x + 2) * (x - 4);
    Ef := Kernel(Tm + 4 * IdentityMatrix(Rationals(), 3));
    Wf := [ &+[ v[i] * Wb[i] : i in [1..3] ] : v in Basis(Ef) ];
    assert #Wf eq 1;
    cmds := Sort(Setseq(&join{ Set(cmd[n]) : n in [1..B] }));
    function Certified(target)
        tg := AssociativeArray(); for t in target do tg[t[1]] := t[2]; end for;
        assert forall{t : t in target | t[1] in cmds};
        nrows := #Wf + #cmds;
        A := ZeroMatrix(Rationals(), nrows, B); rhs := ZeroMatrix(Rationals(), nrows, 1);
        for i->g in Wf do for n in [1..B] do A[i, n] := g[n]; end for; end for;
        for j->d in cmds do
            for n in [1..B] do if d in cmd[n] then A[#Wf + j, n] +:= 1; end if; end for;
            rhs[#Wf + j, 1] := IsDefined(tg, d) select Stab(d) * tg[d] else 0;
        end for;
        return IsConsistent(Transpose(A), Transpose(rhs));
    end function;
    assert Certified([<-228, 1>, <-760, 1>, <-12, -1>, <-40, -1>]);
    assert not Certified([<-228, 1>, <-12, 1>, <-760, -1>, <-40, -1>]);
    assert not Certified([<-228, 1>, <-40, 1>, <-760, -1>, <-12, -1>]);
    assert not Certified([<-228, 1>, <-760, -1>]) and not Certified([<-12, 1>, <-40, -1>]);
    assert not Certified([<-532, 1>, <-40, -1>]);
    printf " [dim W = 3, dim W_f = 1; one pairing certified, five controls refused] Done!\n";
end procedure;

test_hyperelliptic_15_38_torsion();
