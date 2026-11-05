// X_0^6(649)* IS HYPERELLIPTIC, PART ONE: THE TORSION OF ITS JACOBIAN AND THE REDUCTION OF THE THREE FIBRE DIFFERENCES.
// About ten minutes.
//
// Y = X_0^6(649)* = X_0^6(649)/W (W the full Atkin-Lehner group, order 16) has genus three and no further quotient, so
// it is hyperelliptic exactly when it carries two distinct linearly equivalent divisors of degree two.  Three candidates
// come from the Shimizu screen of the invariant differentials at the rational CM points (not used by the proof):
//     D_1 = P_-24 + P_-132,   D_2 = P_-40 + P_-187,   D_3 = P_-52 + P_-660.
// tests/_offline/Hyperelliptic6_649.m certifies, by the dual test of paper/hyperelliptic-15-38.tex, Section 2, that the
// three differences D_i - D_j are TORSION in J(Y)(Q).  This test shows that is enough (paper/hyperelliptic-6-649.tex):
//   (a) J(Y)(Q)_tors has order dividing 6: the orders #J(Y)(F_p) = L_p(1), from the point counts over F_p, F_{p^2},
//       F_{p^3} by the trace formula, have gcd 6 over the good primes below 120;
//   (b) the 3-part of J(Y)(Q)_tors injects into the component group Phi_3 of the Neron model at 3: Y has totally
//       degenerate reduction at 3 (a quotient of the Cerednik-Drinfeld model), the formal group of its Jacobian over Z_3
//       is torsion-free, and the identity component of the special fibre is a split torus, whose F_3-points (F_3^*)^3
//       have no 3-torsion.  Phi_3 is the critical group of the quotient by W of the Cerednik-Drinfeld dual graph
//       (vertices: the left ideal classes of the Eichler order of level 649 in the definite quaternion algebra of
//       discriminant 2, in two copies swapped by w_3; singular points: the classes of the level-1947 order; thickness
//       |O^*/{+-1}| times the order of the stabiliser in W); it is Z/102;
//   (c) each D_i - D_j has class ZERO in Phi_3 (and in Phi_2, as a check): a CM point reduces to the component or
//       singular point whose Eichler order admits the optimal embedding of its CM order (Molina, Ribet bimodules and the
//       specialization of Heegner points, Theorems 5.3, 5.4, 5.8), to a singular point when p ramifies in the CM field
//       and to a component when p is inert, and a singular point whose branches are swapped by an involution of W with
//       p | m becomes a smooth point of the merged component.
// So a torsion D_i - D_j has trivial 3-part by (b) and (c), hence lies in the 2-part of J(Y)(Q)_tors, of order at most
// two by (a); among three divisors with pairwise differences in a group of order at most two, two are linearly
// equivalent, and Y is hyperelliptic.
// Orders of ideal classes are not reproducible between Magma runs, so everything is built from one list.

procedure test_hyperelliptic_6_649()
    printf "X_0^6(649)*: the torsion of J(Y)(Q) divides 6, its 3-part injects into Phi_3 = Z/102, and the three fibre differences reduce to zero at 3 and at 2...";
    D := 6; N := 649;
    W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
    Y := CreateShimuraQuot(D, N, W); Y`g := GenusShimuraCurveQuotient(D, N, W); Y`CurveID := 0;
    assert Y`g eq 3;
    curves := GetQuotientsAndGenera([Y]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};
    discs := [-24, -132, -40, -187, -52, -660];
    assert forall{d : d in discs | DegreeOfFieldOfDefinitionOfCMPoint(star, d) eq 1};
    fibres := [ [<-24, 1>, <-132, 1>], [<-40, 1>, <-187, 1>], [<-52, 1>, <-660, 1>] ];

    // (a) the torsion bound
    g := 0;
    for p in [q : q in PrimesInInterval(5, 120) | (D*N) mod q ne 0] do
        s := [ p^i + 1 - ComputePointsViaTrace(Y, p, i) : i in [1..3] ];
        e := [Rationals() | 1];
        for k in [1..3] do Append(~e, (&+[ (-1)^(i-1) * e[k-i+1] * s[i] : i in [1..k] ]) / k); end for;
        // L_p(1) = 1 - e1 + e2 - e3 + p e2 - p^2 e1 + p^3
        g := GCD(g, Integers()!(1 - e[2] + e[3] - e[4] + p*e[3] - p^2*e[2] + p^3));
    end for;
    assert g eq 6;

    // (b), (c) the component groups at 3 and 2 and the classes of the differences
    function OptEmb(d, Ov)   // optimal embeddings of the order of discriminant d into the definite order Ov, up to O^*
        t0 := d mod 2; n := (t0^2 - d) div 4;
        Zb := ZBasis(Ov);
        Gm := Matrix(Rationals(), 4, 4, [ Trace(Zb[i] * Conjugate(Zb[j])) : i, j in [1..4] ]);
        S := ShortVectors(LatticeWithGram(Gm), 2*n, 2*n);
        d0 := FundamentalDiscriminant(d); _, f := IsSquare(d div d0);
        opt := [];
        for s in S do for sg in [1, -1] do
            x := &+[ sg * s[1][i] * Zb[i] : i in [1..4] ];
            if Trace(x) ne t0 then continue; end if;
            if exists{<q, c> : c in [0..q-1], q in PrimeDivisors(f) | (x - c)/q in Ov} then continue; end if;
            Append(~opt, x);
        end for; end for;
        U := Units(Ov);
        return #{ { u * x * u^-1 : u in U } : x in opt };
    end function;
    function Stab(d)
        e := d eq -3 select 3 else (d eq -4 select 2 else 1);
        if exists{mm : mm in Divisors(D*N) | mm gt 1 and GCD(mm, (D*N) div mm) eq 1 and (d eq -4*mm or (mm mod 4 eq 3 and d eq -mm))} then e *:= 2; end if;
        return e;
    end function;
    function Find(J, reps) hits := [i : i in [1..#reps] | IsIsomorphic(J, reps[i])]; assert #hits eq 1; return hits[1]; end function;
    for p in [3, 2] do
        Q := QuaternionAlgebra(D div p); O := QuaternionOrder(Q, N); OO := Order(O, p);
        L := LeftIdealClasses(O); LL := LeftIdealClasses(OO);
        h := #L; hp := #LL;
        org := [Find(lideal< O | Basis(II) >, L) : II in LL];
        alE := []; alV := [];
        for l in PrimeDivisors(D*N) do
            P := PrimeIdeal(OO, l); assert Norm(P) eq l and IsTwoSidedIdeal(P);
            Append(~alE, <l, [Find(P*II, LL) : II in LL]>);
            if l ne p then
                P0 := PrimeIdeal(O, l); assert Norm(P0) eq l and IsTwoSidedIdeal(P0);
                Append(~alV, <l, [Find(P0*I, L) : I in L]>);
            end if;
        end for;
        uV := [#UnitGroup(RightOrder(I)) : I in L]; uE := [#UnitGroup(RightOrder(II)) : II in LL];
        ValidateDualGraphData(D, N, p, <h, hp, org, alE, alV, uV, uE>);
        img := func< list, q | [t[2] : t in list | t[1] eq q][1] >;
        function vact(mm, x)
            lab := ((x - 1) mod h) + 1; sgn := x le h;
            for q in PrimeDivisors(mm) do if q eq p then sgn := not sgn; else lab := img(alV, q)[lab]; end if; end for;
            return sgn select lab else lab + h;
        end function;
        function eact(mm, e) for q in PrimeDivisors(mm) do e := img(alE, q)[e]; end for; return e; end function;
        wp := img(alE, p);
        vorb := [0 : x in [1..2*h]]; nv := 0;
        for x in [1..2*h] do if vorb[x] eq 0 then nv +:= 1; for mm in W do vorb[vact(mm, x)] := nv; end for; end if; end for;
        eorb := [0 : e in [1..hp]]; edges := []; thick := []; nedge := 0;
        for e in [1..hp] do
            if eorb[e] ne 0 then continue; end if;
            S := { eact(mm, e) : mm in W };
            if exists{mm : mm in W | mm mod p eq 0 and eact(mm, e) eq e} then
                assert vorb[org[e]] eq vorb[org[wp[e]] + h];      // the ends are swapped: a smooth point of Y
                for x in S do eorb[x] := -1; end for; continue;
            end if;
            nedge +:= 1; for x in S do eorb[x] := nedge; end for;
            Append(~edges, <vorb[org[e]], vorb[org[wp[e]] + h]>);
            Append(~thick, (uE[e] div 2) * (#W div #S));
        end for;
        // the regular model subdivides a singular point of thickness t into t - 1 components
        V := nv; Ed := []; chain := [];
        for i->e in edges do
            prev := e[1]; inner := [];
            for k in [1..thick[i]-1] do V +:= 1; Append(~inner, V); Append(~Ed, <prev, V>); prev := V; end for;
            Append(~Ed, <prev, e[2]>); Append(~chain, inner);
        end for;
        Lap := ZeroMatrix(Integers(), V, V);
        for e in Ed do Lap[e[1],e[1]] +:= 1; Lap[e[2],e[2]] +:= 1; Lap[e[1],e[2]] -:= 1; Lap[e[2],e[1]] -:= 1; end for;
        Sm, Pm := SmithForm(Lap);
        inv := [Sm[i,i] : i in [1..V] | Sm[i,i] notin {0, 1}];
        assert inv eq (p eq 3 select [102] else [3, 3, 90]);
        assert #[i : i in [1..V] | Sm[i,i] eq 0] eq 1;          // the graph is connected
        assert #edges - nv + 1 eq 3;                            // the first Betti number is the genus: totally degenerate
        // the star points of discriminant d per vertex of the subdivided graph
        function Place(d)
            d0 := FundamentalDiscriminant(d); sym := KroneckerSymbol(d0, p);
            assert sym ne 1;
            cnt := AssociativeArray(); total := 0;
            if sym eq -1 then
                for i in [1..h] do
                    c := OptEmb(d, RightOrder(L[i]));
                    if c eq 0 then continue; end if;
                    for x in [i, i + h] do v := vorb[x]; cnt[v] := (IsDefined(cnt, v) select cnt[v] else 0) + c; total +:= c; end for;
                end for;
            else
                for e in [1..hp] do
                    c := OptEmb(d, RightOrder(LL[e]));
                    if c eq 0 then continue; end if;
                    o := eorb[e];
                    if o eq -1 then v := vorb[org[e]];
                    else t := thick[o]; assert IsEven(t); v := t eq 2 select chain[o][1] else chain[o][t div 2]; end if;
                    cnt[v] := (IsDefined(cnt, v) select cnt[v] else 0) + c; total +:= c;
                end for;
            end if;
            ratio := #W div Stab(d);                               // CM points of X above one star point
            assert total eq ratio;
            for v in Keys(cnt) do assert cnt[v] mod ratio eq 0; cnt[v] := cnt[v] div ratio; end for;
            return cnt;
        end function;
        places := AssociativeArray(); for d in discs do places[d] := Place(d); end for;
        function ClassOf(dv)
            v := ZeroMatrix(Integers(), V, 1);
            for t in dv do for w in Keys(places[t[1]]) do v[w, 1] +:= t[2] * places[t[1]][w]; end for; end for;
            assert &+Eltseq(v) eq 0;
            w := Pm * v;
            return [ w[i,1] mod Sm[i,i] : i in [1..V] | Sm[i,i] notin {0, 1} ];
        end function;
        for i in [1..3], j in [i+1..3] do
            assert ClassOf(fibres[i] cat [ <t[1], -t[2]> : t in fibres[j] ]) eq [ 0 : k in inv ];
        end for;
    end for;
    printf " Done!\n";
end procedure;

test_hyperelliptic_6_649();
