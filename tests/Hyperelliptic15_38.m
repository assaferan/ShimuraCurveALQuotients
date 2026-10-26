// X_0^15(38)/<w_2, w_15, w_19> IS HYPERELLIPTIC, PART ONE: THE STAR, ITS TORSION, AND THE SPECIALISATION OF THE
// BRANCH POINTS, AND WHICH DOUBLE COVER Y IS.  About a minute.
//
// Y = X_0^15(38)/H, H = <w_2, w_15, w_19>, has genus three and is a double cover of the genus-one star
// E = X_0^15(38)/W branched at the four rational CM points P_-228, P_-760, P_-12, P_-40 (the fixed points of the
// involutions outside H).  Y is hyperelliptic if and only if two of the branch points add up to the other two
// in the group E, and tests/_offline/Hyperelliptic15_38.m certifies that
//     Delta = P_-228 + P_-760 - P_-12 - P_-40
// is TORSION in E(Q).  This test shows it is zero (paper/hyperelliptic-15-38.tex):
//   (a) E is the curve 570e4 of Cremona's tables: its traces of Frobenius agree with the star's point counts at
//       every good prime up to the Sturm bound 240 of level 570, and its component groups at 3 and 5 have
//       orders 8 and 4, which single it out in its isogeny class;
//   (b) E(Q)_tors = Z/4, and neither the generator nor its double reduces to the identity component at 3 or
//       at 5, so the torsion injects into the component group Phi_3 (and into Phi_5);
//   (c) the class of Delta in Phi_3 and in Phi_5 is zero.  The component groups are the critical groups of
//       the quotient by W of the Cerednik-Drinfeld dual graph (vertices: the left ideal classes of the Eichler
//       order of level 38 in the definite quaternion algebra of discriminant 15/p, in two copies swapped by
//       w_p; singular points: the classes of the level-38p order; thickness |O^*/{+-1}| times the order of the
//       stabiliser in W).  A CM point reduces to the component or singular point whose Eichler order admits
//       the optimal embedding of its CM order (Molina, Ribet bimodules and the specialization of Heegner
//       points): to a singular point when p ramifies in the CM field, to a component when p is inert.  A
//       singular point whose branches are swapped by an involution of W with p | m becomes a smooth point of
//       the merged component on E.
// Orders of ideal classes are not reproducible between Magma runs, so everything is built from one list.

procedure test_hyperelliptic_15_38()
    printf "X_0^15(38)/<w_2,w_15,w_19>: the star is 570e4, its torsion injects into Phi_3, and Delta specialises to zero...";
    D := 15; N := 38;
    W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
    Xstar := CreateShimuraQuot(D, N, W);
    Xstar`g := GenusShimuraCurveQuotient(D, N, W); Xstar`CurveID := 0;
    assert Xstar`g eq 1;
    assert GenusShimuraCurveQuotient(D, N, {1, 2, 15, 19, 30, 38, 285, 570}) eq 3;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};
    branch := [-228, -760, -12, -40];
    assert forall{d : d in branch | DegreeOfFieldOfDefinitionOfCMPoint(star, d) eq 1};

    // (a) the isogeny class: a_p for every good prime up to the Sturm bound 240 = (2/12) [SL_2(Z) : Gamma_0(570)]
    E := MinimalModel(EllipticCurve("570e4"));
    assert Conductor(E) eq 570;
    for p in [q : q in PrimesInInterval(2, 240) | 570 mod q ne 0] do
        assert p + 1 - ComputePointsViaTrace(Xstar, p, 1) eq TraceOfFrobenius(E, p);
    end for;
    assert [TamagawaNumber(E, p) : p in [3, 5]] eq [8, 4];   // Cremona: I8 at 3, I4 at 5; 570e1-3 have I2, I4, I2 at 3

    // (b) the torsion and the identity component at 3 and 5
    G, m := TorsionSubgroup(E);
    assert Invariants(G) eq [4];
    T := m(G.1);
    for p in [3, 5] do
        Ep := ChangeRing(E, GF(p));
        sing := SingularPoints(Curve(Ep));
        inE0 := func< P | not (Curve(Ep)!(Ep(GF(p))![GF(p)!c : c in Eltseq(P)]) in sing) >;
        assert not inE0(T) and not inE0(2*T);
    end for;

    // (c) the component groups and the specialisation of the branch points
    function OptEmb(d, Ov)   // optimal embeddings of the order of discriminant d into the definite order Ov
        t0 := d mod 2; n := (t0^2 - d) div 4;
        Zb := ZBasis(Ov);
        Gm := Matrix(Rationals(), 4, 4, [ Trace(Zb[i] * Conjugate(Zb[j])) : i, j in [1..4] ]);
        S := ShortVectors(LatticeWithGram(Gm), 2*n, 2*n);
        d0 := FundamentalDiscriminant(d); _, f := IsSquare(d div d0);
        cnt := 0;
        for s in S do for sg in [1, -1] do
            x := &+[ sg * s[1][i] * Zb[i] : i in [1..4] ];
            if Trace(x) ne t0 then continue; end if;
            if exists{<q, c> : c in [0..q-1], q in PrimeDivisors(f) | (x - c)/q in Ov} then continue; end if;
            cnt +:= 1;
        end for; end for;
        return cnt;
    end function;
    function Find(J, reps) hits := [i : i in [1..#reps] | IsIsomorphic(J, reps[i])]; assert #hits eq 1; return hits[1]; end function;
    for p in [3, 5] do
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
        uE := [#UnitGroup(RightOrder(II)) : II in LL];
        assert forall{u : u in uE | u eq 2} and forall{I : I in L | #UnitGroup(RightOrder(I)) eq 2};
        ValidateDualGraphData(D, N, p, <h, hp, org, alE, alV, [2 : I in L], uE>);
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
                assert vorb[org[e]] eq vorb[org[wp[e]] + h];      // the ends are swapped: a smooth point of E
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
        Sm, Pm, Qm := SmithForm(Lap);
        inv := [Sm[i,i] : i in [1..V] | Sm[i,i] notin {0, 1}];
        assert inv eq [TamagawaNumber(E, p)];                    // Phi_3 = Z/8, Phi_5 = Z/4
        function Place(d)
            d0 := FundamentalDiscriminant(d); sym := KroneckerSymbol(d0, p);
            assert sym ne 1;
            if sym eq -1 then
                vs := { vorb[i] : i in [1..h] | OptEmb(d, RightOrder(L[i])) ne 0 };
                assert #vs eq 1; return Rep(vs);
            end if;
            es := [ e : e in [1..hp] | OptEmb(d, RightOrder(LL[e])) ne 0 ];
            assert #{ eorb[e] : e in es } eq 1 and eorb[es[1]] eq -1;   // on these curves every ramified branch point lands on a swapped singular point
            return vorb[org[es[1]]];
        end function;
        pt := AssociativeArray(); for d in branch do pt[d] := Place(d); end for;
        function Vec(dv) v := ZeroMatrix(Integers(), V, 1); for t in dv do v[t[1], 1] +:= t[2]; end for; return v; end function;
        function ClassOf(v)   // coordinates of the class of a degree-zero divisor in the Smith cokernel
            w := Pm * v;
            return [ w[i,1] mod Sm[i,i] : i in [1..V] | Sm[i,i] notin {0, 1} ];
        end function;
        Delta := Vec([<pt[-228], 1>, <pt[-760], 1>, <pt[-12], -1>, <pt[-40], -1>]);
        assert ClassOf(Delta) eq [0];
        // the other two pairings are not even torsion; at p = 3 both have nonzero image
        if p eq 3 then
            assert ClassOf(Vec([<pt[-228], 1>, <pt[-12], 1>, <pt[-760], -1>, <pt[-40], -1>])) ne [0];
            assert ClassOf(Vec([<pt[-228], 1>, <pt[-40], 1>, <pt[-760], -1>, <pt[-12], -1>])) ne [0];
        end if;

        // (d) which of the two double covers of E branched at the four points is Y.  A double cover of E
        // branched exactly there is E(sqrt(c g)) with g = x, the function with divisor Delta, or g = x f,
        // f the function with divisor 2 S - 2 P_-40 where S - P_-40 is the rational 2-torsion class
        // (c a rational constant).  On the regular model at p the divisor of g is Delta plus a vertical
        // part a (the solution of the Laplace equation), and over a singular point of E of thickness t
        // the cover is etale when the slope of a along it is even and a single node of half the thickness
        // (fixed by the deck involution with its branches) when the slope is odd; a thickness-one node
        // with odd slope makes the fibre non-reduced, so no such cover is semistable at p.  The dual graph
        // of Y, the quotient of the Cerednik-Drinfeld graph by H, must show the first pattern: it does
        // for g = x, and g = x f has odd slope at a thickness-one node.
        function Solve(v)   // integral a with Lap a = v (v of degree zero and in the image)
            w := Pm * v; y := ZeroMatrix(Integers(), V, 1);
            for i in [1..V] do
                if Sm[i,i] eq 0 then assert w[i,1] eq 0; else assert w[i,1] mod Sm[i,i] eq 0; y[i,1] := w[i,1] div Sm[i,i]; end if;
            end for;
            a := Qm * y; assert Lap * a eq v; return a;
        end function;
        ax := Solve(Delta);
        if p eq 3 then
            two := [ v : v in [1..V] | ClassOf(Vec([<v, 1>, <pt[-40], -1>])) eq [TamagawaNumber(E, p) div 2] ];
            assert not IsEmpty(two);
            af := Solve(Vec([<two[1], 2>, <pt[-40], -2>]));
            // the dual graph of Y: the H-orbits of the same edges
            Hgrp := {1, 2, 15, 19, 30, 38, 285, 570};
            eorbH := [0 : e in [1..hp]]; nH := 0;
            for e in [1..hp] do
                if eorbH[e] ne 0 then continue; end if;
                S := { eact(mm, e) : mm in Hgrp };
                if exists{mm : mm in Hgrp | mm mod p eq 0 and eact(mm, e) eq e} then for x in S do eorbH[x] := -1; end for; continue; end if;
                nH +:= 1; for x in S do eorbH[x] := nH; end for;
            end for;
            oddslope := 0;
            for i->e in edges do
                reps := [ ee : ee in [1..hp] | eorb[ee] eq i ];
                above := #{ eorbH[ee] : ee in reps };
                assert above in {1, 2} and 0 notin { eorbH[ee] : ee in reps };
                sx := (ax[e[2],1] - ax[e[1],1]) / thick[i]; sf := ((ax + af)[e[2],1] - (ax + af)[e[1],1]) / thick[i];
                assert IsIntegral(sx) and IsIntegral(sf);
                assert (above eq 1) eq IsOdd(Integers()!sx);       // Y has the reduction of the cover E(sqrt(c x))
                if thick[i] eq 1 and IsOdd(Integers()!sf) then oddslope +:= 1; end if;
            end for;
            assert oddslope gt 0;                                    // the cover E(sqrt(c x f)) is not semistable at 3
        end if;
    end for;
    printf " Done!\n";
end procedure;

test_hyperelliptic_15_38();
