// The dual test of tests/_offline/Hyperelliptic6_385.m from a SAVED basis (basisS.m writes basis_<M>_<B>.txt):
//   magma -b BASIS:=basis_4620_4000.txt tools/proof-6-385/withbasis.m
// Minutes once the basis exists; the offline test recomputes the basis instead.
//
// On the genus-one star E = X_0^6(385)/W let Z(n) = sum_r (1/e_d) P_d over the discriminants d = -4n/r^2 with CM
// points (Guo-Yang Lemma 25).  For a linear form phi on Pic(E)_Q vanishing on the Hodge class, Phi_phi =
// sum_n phi(Z(n)) q^n is a cusp form of weight 3/2 on Gamma_0(4620) with the theta^3 multiplier vanishing at every
// index carrying no CM point, and it lies in the a_13(E) = -6 eigenspace W_f of T(13^2) on the space W of such
// forms (paper/hyperelliptic-15-38.tex, Section 2).  A degree-zero CM divisor is torsion as soon as some principal
// part with that divisor pairs to zero with W; the test uses all of W, which is the sufficient direction.  This test certifies Delta = P_-40 + P_-1540 - Z_-9240 and refuses the other
// pairings of the branch places and the differences within them (paper/hyperelliptic-6-385.tex).

AttachSpec("ShimuraQuotients.spec");
procedure test_hyperelliptic_6_385_torsion()
    printf "X_0^6(385)*: P_-40 + P_-1540 - Z_-9240 is torsion (level 4620, the a_13 = -6 eigenspace)...";
    D := 6; N := 385; M := 4620; B := 4000;
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
    // T(p^2) on W: a(p^2 n) + (-n/p) a(n) + p a(n/p^2), read on the indices up to B/p^2
    p := 13; Bp := B div p^2;
    Tp := func< g | [ g[p^2*n] + KroneckerSymbol(-n, p) * g[n] + (n mod p^2 eq 0 select p * g[n div p^2] else 0) : n in [1..Bp] ] >;
    Wmat := Matrix(Rationals(), [ [ g[n] : n in [1..Bp] ] : g in Wb ]);
    ok, Tm := IsConsistent(Wmat, Matrix(Rationals(), [ Tp(g) : g in Wb ]));
    determined := ok and Rank(Wmat) eq #Wb;      // T(p^2) on W is read off the indices up to B/p^2, which must separate W
    // every form of W is used (the sufficient direction: the curve's own form lies in W); when the coefficients determine
    // T(p^2) on W, its -6 eigenspace must be a line, the curve's
    if determined then
        Ef := Kernel(Tm - -6 * IdentityMatrix(Rationals(), #Wb));
        assert Dimension(Ef) eq 1;
    end if;
    Wf := Wb;
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
    assert Certified([<-40, 1>, <-1540, 1>, <-9240, -1>]);
    assert not Certified([<-40, 1>, <-1540, -1>]);
    assert not Certified([<-40, 2>, <-9240, -1>]);
    printf " [dim W = %o, T(p^2) on it %o; the pairing certified, two controls refused] Done!
", #Wb, determined select "determined, with the curve's eigenvalue on a line" else "not determined by the coefficients available";
end procedure;

test_hyperelliptic_6_385_torsion();
quit;
