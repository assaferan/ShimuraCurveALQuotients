// The dual test of tests/_offline/Hyperelliptic3390_1.m from a SAVED basis (basisS.m writes basis_<M>_<B>.txt):
//   magma -b BASIS:=basis_6780_1100.txt tools/proof-3390-1/withbasis.m
// Seconds once the basis exists; the offline test recomputes the basis instead.
//
// Y = X_0^3390(1)/H, H = <w_5, w_6, w_113>, has genus four and is a double cover of the genus-two star X* = X_0^3390(1)/W,
// branched at the two rational CM points P_-3 and P_-40 (the fixed points of w_3 and w_10, the only involutions outside
// H with fixed points; tests/Hyperelliptic3390_1.m).  Y is hyperelliptic if and only if P_-3 + P_-40 is a fibre of the
// hyperelliptic map of X* (paper/star-2002-proof.tex, Lemma 3), i.e. iff P_-3 + P_-40 - F is principal for a fibre F.
// The rational torsion of J(X*) is trivial (tests/Hyperelliptic3390_1.m), so principal is the same as torsion, and
// torsion is decided by the dual test of paper/hyperelliptic-15-38.tex, Section 2: on X* let Z(n) = sum_r (1/e_d) P_d
// over the discriminants d = -4n/r^2 with CM points; for a linear form phi on Pic(X*)_Q vanishing on the Hodge class,
// sum_n phi(Z(n)) q^n is a cusp form of weight 3/2 on Gamma_0(6780) with the theta^3 multiplier vanishing at every
// index carrying no CM point.  The space W of such forms has dimension two, the genus, and T(49) acts on it with the
// characteristic polynomial of the star's a_7, so W is the star's; a degree-zero CM divisor is torsion as soon as some
// principal part with that divisor pairs to zero with W.  This test certifies P_-3 + P_-40 - P_-43 - P_-67 and
// P_-3 + P_-40 - P_-123 - P_-168 (two fibres found numerically by the Shimizu screen, which the test does not use), and
// refuses the non-fibres P_-3 - P_-40 and P_-3 + P_-43 - P_-40 - P_-67.

AttachSpec("ShimuraQuotients.spec");
procedure test_hyperelliptic_3390_1_fibre()
    printf "X_0^3390(1)/<w_5, w_6, w_113>: P_-3 + P_-40 is a fibre of the star's g^1_2 (level 6780, the two-dimensional W)...";
    D := 3390; N := 1; M := 6780; B := 1100;
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
    A0 := eval Read(BASIS);               // the cusp forms of weight 3/2 on Gamma_0(M) to precision B, one per row
    assert Ncols(A0) ge B; A0 := Submatrix(A0, 1, 1, Nrows(A0), B); Bq := [1..Nrows(A0)];
    cmd := AssociativeArray();
    for n in [1..B] do cmd[n] := [ d : d in discs(n) | HasCM(d) ]; end for;
    empty := [ n : n in [1..B] | IsEmpty(cmd[n]) ];
    K := Kernel(Matrix(Rationals(), #Bq, #empty, [ A0[i, n] : n in empty, i in [1..#Bq] ]));
    Wb := [ Vector(k) * A0 : k in Basis(K) ];
    assert #Wb eq 2;                       // the genus of the star
    // T(p^2) on W: a(p^2 n) + (-n/p) a(n) + p a(n/p^2), read on the indices up to B/p^2; its characteristic polynomial
    // is prod (x - a_p(f)) over the two newforms of the star, i.e. x^2 + c1 x + (c2 - 2p) for L_p(x) = 1 + c1 x + c2 x^2 + ...
    W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
    Xstar := CreateShimuraQuot(D, N, W); Xstar`g := GenusShimuraCurveQuotient(D, N, W); Xstar`CurveID := 0;
    p := 7; Bp := B div p^2;
    s1 := ComputePointsViaTrace(Xstar, p, 1); s2 := ComputePointsViaTrace(Xstar, p, 2);
    c1 := s1 - p - 1; c2 := (s2 - p^2 - 1 + c1^2) div 2;
    Tp := func< g | [ g[p^2*n] + KroneckerSymbol(-n, p) * g[n] + (n mod p^2 eq 0 select p * g[n div p^2] else 0) : n in [1..Bp] ] >;
    Wmat := Matrix(Rationals(), [ [ g[n] : n in [1..Bp] ] : g in Wb ]);
    ok, Tm := IsConsistent(Wmat, Matrix(Rationals(), [ Tp(g) : g in Wb ]));
    assert ok;
    R<x> := PolynomialRing(Rationals());
    assert CharacteristicPolynomial(Tm) eq x^2 + c1 * x + (c2 - 2*p);
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
    assert Certified([<-3, 1>, <-40, 1>, <-43, -1>, <-67, -1>]);
    assert Certified([<-3, 1>, <-40, 1>, <-123, -1>, <-168, -1>]);
    assert Certified([<-43, 1>, <-67, 1>, <-123, -1>, <-168, -1>]);
    assert not Certified([<-3, 1>, <-40, -1>]);
    assert not Certified([<-3, 1>, <-43, 1>, <-40, -1>, <-67, -1>]);
    printf " [dim W = 2 with the star's T(49); two fibre relations certified, two controls refused] Done!\n";
end procedure;

test_hyperelliptic_3390_1_fibre();
quit;
