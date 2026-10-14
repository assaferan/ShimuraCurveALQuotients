// The Atkin-Lehner ideal of norm q = l^k of an Eichler order O of level divisible by l^k.
// O_l = End(L_a) meet End(L_b) for the two ends of a chain of length k in the l-adic tree; a maximal
// order O1 containing O is End(L_0) for a vertex L_0 on the chain, at distances j_a, j_b from the
// ends (j_a + j_b = k).  In coordinates with L_0 = Z_l^2 and the ends cut out by the functionals
// y_1 = 0 mod l^j_a and y_2 = 0 mod l^j_b, the reflection of the chain is
// Pi = [[0, l^j_a], [l^j_b, 0]], of norm l^k, and the ideal is {x in O : x in Pi O_l}.
SetColumns(0);
function al_ideal(O, l, q)
    k := Valuation(q, l);
    O1 := MaximalOrder(O);
    M, phi, phiinv := pMatrixRing(O1, l : Precision := k + 8);
    R := Integers(l^k);
    B := Basis(O);
    ints := func<x | [Integers()!(Rationals()!c) : c in Eltseq(x)]>;
    imgs := [ints(phi(b)) : b in B];
    // stabilised functionals mod l^j, grouped by direction (reduction mod l)
    best := AssociativeArray();                     // direction -> <j, f>
    for j in [1..k] do
        Rj := Integers(l^j);
        for f1, f2 in [0..l^j-1] do
            if GCD(GCD(f1, f2), l) ne 1 then continue; end if;
            f := Matrix(Rj, 1, 2, [f1, f2]);
            if forall{x : x in imgs | (f*Matrix(Rj, 2, 2, x))[1,2]*f[1,1] eq (f*Matrix(Rj, 2, 2, x))[1,1]*f[1,2]} then
                d := [GF(l)!f1, GF(l)!f2]; d := d[1] ne 0 select [1, d[2]/d[1]] else [0, 1];
                if not IsDefined(best, d) or best[d][1] lt j then best[d] := <j, [f1, f2]>; end if;
            end if;
        end for;
    end for;
    dirs := [best[d] : d in Keys(best)];
    assert #dirs in {1, 2} and &+[t[1] : t in dirs] eq k;
    if #dirs eq 1 then
        fb := dirs[1][2]; jb := dirs[1][1]; ja := 0;
        fa := fb[2] mod l ne 0 select [1, 0] else [0, 1];      // any functional independent of fb mod l
    else
        fa := dirs[1][2]; ja := dirs[1][1]; fb := dirs[2][2]; jb := dirs[2][1];
    end if;
    g := Matrix(Rationals(), 2, 2, [fa[1], fa[2], fb[1], fb[2]]);
    assert Valuation(Determinant(g), l) eq 0;
    Pi := g^-1 * Matrix(Rationals(), 2, 2, [0, l^ja, l^jb, 0]) * g;
    assert forall{c : c in Eltseq(Pi) | Valuation(c, l) ge 0};
    // the lattice {n in Z^4 : sum n_i phi(b_i) in Pi O_l}, determined modulo l^k O_l
    G := AbelianGroup([l^k : i in [1..4]]);
    A := hom<G -> G | [G![R!c : c in x] : x in imgs]>;
    Pimgs := [[Integers()!(Rationals()!c) : c in Eltseq(Pi * Matrix(Rationals(), 2, 2, x))] : x in imgs];
    T := sub<G | [G![R!c : c in x] : x in Pimgs]>;
    pre := T @@ A;
    gens := [&+[(Integers()!Eltseq(G!h)[i]) * B[i] : i in [1..4]] : h in Generators(pre)] cat [l^k * b : b in B];
    P := ideal<O | gens>;
    assert Norm(P) eq q and P*P eq q*O and LeftOrder(P) eq O and RightOrder(P) eq O;
    return P;
end function;
AttachSpec("ShimuraQuotients.spec");
Q := QuaternionAlgebra(2); O3 := QuaternionOrder(Q, 3);
printf "level 3: al_ideal = PrimeIdeal: %o\n", al_ideal(O3, 3, 3) eq PrimeIdeal(O3, 3);
function key(I) return MinkowskiGramReduction(GramMatrix(I)/Norm(I) : Canonical := true); end function;
function buckets(reps) b := AssociativeArray(); for i->I in reps do kk := key(I); if IsDefined(b, kk) then Append(~b[kk], i); else b[kk] := [i]; end if; end for; return b; end function;
function find(J, reps, bk)
    ok, cand := IsDefined(bk, key(J));
    if ok then hits := [i : i in cand | IsIsomorphic(J, reps[i])]; if #hits eq 1 then return hits[1]; end if; end if;
    hits := [i : i in [1..#reps] | IsIsomorphic(J, reps[i])]; assert #hits eq 1; return hits[1];
end function;
function build(D, N, p)
    Dp := D div p; Qa := QuaternionAlgebra(Dp); O := QuaternionOrder(Qa, N); OO := Order(O, p);
    L := LeftIdealClasses(O); LL := LeftIdealClasses(OO); bL := buckets(L); bLL := buckets(LL);
    org := [find(lideal< O | Basis(II) >, L, bL) : II in LL];
    alE := []; alV := [];
    for l in PrimeDivisors(Dp*N*p) do
        q := l^Valuation(D*N, l);
        P := Valuation(q, l) eq 1 select PrimeIdeal(OO, l) else al_ideal(OO, l, q);
        Append(~alE, <q, [find(P*II, LL, bLL) : II in LL]>);
        if l ne p then P0 := Valuation(q, l) eq 1 select PrimeIdeal(O, l) else al_ideal(O, l, q); Append(~alV, <q, [find(P0*I, L, bL) : I in L]>); end if;
    end for;
    return <#L, #LL, org, alE, alV, [#UnitGroup(RightOrder(I)) : I in L], [#UnitGroup(RightOrder(II)) : II in LL]>;
end function;
curves := GetHyperellipticCandidates();
for t in [<22, 9, 11>, <22, 9, 2>, <6, 25, 3>, <6, 25, 2>, <15, 16, 3>, <14, 45, 7>, <10, 63, 5>, <6, 35, 3>] do
    D := t[1]; N := t[2]; p := t[3];
    t0 := Realtime();
    data := build(D, N, p);
    nok := 0; bad := [];
    for X in curves do
        if X`D ne D or X`N ne N then continue; end if;
        nv, edges := DualGraphQuotientFromData(data, p, X`W);
        b1 := #edges eq 0 select 0 else #edges - nv + 1;
        if b1 eq X`g then nok +:= 1; else Append(~bad, <Sort(SetToSequence(X`W)), b1, X`g>); end if;
    end for;
    printf "X_0(%o,%o) p=%o: h=%o h'=%o genus %o; quotient identities %o agree, %o disagree %o  (%.0os)\n", D, N, p, data[1], data[2], data[2]-2*data[1]+1, nok, #bad, bad, Realtime()-t0;
end for;
quit;
