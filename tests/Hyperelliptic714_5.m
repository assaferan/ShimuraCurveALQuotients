// X_0^714(5)/<w_3, w_10, w_14, w_17> IS HYPERELLIPTIC, PART ONE: THE STAR, ITS TORSION, THE SPECIALISATION OF THE BRANCH
// PLACES, AND WHICH DOUBLE COVER Y IS.  About three minutes.
//
// Y = X_0^714(5)/H, H = <w_3, w_10, w_14, w_17>, has genus three and is a double cover of the genus-one star E = X_0^714(5)/W
// branched at the CM places of discriminants -340, -840 (rational points) and -14280 (a place of degree two, Z) (the fixed points of the involutions outside H).  Y is
// hyperelliptic if and only if two branch points add up to the other two in the group E (paper/hyperelliptic-15-38.tex,
// Section 1), and tests/_offline/Hyperelliptic714_5.m certifies that
//     Delta = P_-340 + P_-840 - Z_-14280
// is TORSION in E(Q).  This test shows it is zero and that Y is the hyperelliptic cover (paper/hyperelliptic-714-5.tex):
//   (a) E is the curve 3570x2 of Cremona's tables: its traces of Frobenius agree with the star's point counts at every
//       good prime up to the Sturm bound 1728 of its level, and its component groups at 3, 7 and 17 have orders
//       10, 4 and 2, which single it out in its isogeny class;
//   (b) E(Q)_tors = Z/2 x Z/2 and each nonzero torsion point reduces off the identity component at one of 3, 7, 17, so the torsion injects into Phi_3 x Phi_7 x Phi_17;
//   (c) the class of Delta in each of these component groups is zero.  The component groups are the critical groups
//       of the quotient by W of the Cerednik-Drinfeld dual graph (vertices: the left ideal classes of the Eichler
//       order of level 5 in the definite quaternion algebra of discriminant 714/p, in two copies swapped by w_p;
//       singular points: the classes of the level-5p order; thickness |O^*/{+-1}| times the order of the stabiliser
//       in W).  A CM point reduces to the component or singular point whose Eichler order admits the optimal
//       embedding of its CM order (Molina, Ribet bimodules and the specialization of Heegner points, Theorems 5.3,
//       5.4, 5.8): to a singular point when p ramifies in the CM field, to a component when p is inert; a singular
//       point whose branches are swapped by an involution of W with p | m becomes a smooth point of the merged
//       component.  A place of degree two is treated as its whole CM cycle;
//   (d) Y is the untwisted cover: the Prym of Y -> E has an irreducible L-polynomial at p = 11, so it is simple
//       over Q, whereas the twist of the hyperelliptic cover by the rational 2-torsion class has a Prym that splits
//       into two elliptic curves over Q (paper/hyperelliptic-15-38.tex, Section 5).
// Orders of ideal classes are not reproducible between Magma runs, so everything is built from one list.

procedure test_hyperelliptic_714_5()
    printf "X_0^714(5)/<w_3, w_10, w_14, w_17>: the star is 3570x2, its torsion injects into the component groups, Delta specialises to zero, the Prym is simple...";
    D := 714; N := 5;
    W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
    H := {1, 3, 10, 14, 17, 30, 35, 42, 51, 105, 170, 238, 510, 595, 714, 1785};
    Xstar := CreateShimuraQuot(D, N, W);
    Xstar`g := GenusShimuraCurveQuotient(D, N, W); Xstar`CurveID := 0;
    assert Xstar`g eq 1;
    Y := CreateShimuraQuot(D, N, H); Y`g := GenusShimuraCurveQuotient(D, N, H); Y`CurveID := 0;
    assert Y`g eq 3;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};
    branch := [-340, -840, -14280];
    assert [ DegreeOfFieldOfDefinitionOfCMPoint(star, d) : d in branch ] eq [1, 1, 2];

    // (a) the isogeny class
    E := MinimalModel(EllipticCurve("3570x2"));
    lev := Conductor(E);
    bound := Floor(2 * lev * &*[Rationals() | 1 + 1/q : q in PrimeDivisors(lev)] / 12);
    assert bound eq 1728;
    for p in [q : q in PrimesInInterval(2, bound) | lev mod q ne 0] do
        assert p + 1 - ComputePointsViaTrace(Xstar, p, 1) eq TraceOfFrobenius(E, p);
    end for;
    primes := [3, 7, 17];
    assert [TamagawaNumber(E, p) : p in primes] eq [10, 4, 2];

    // (b) the torsion and the identity components
    function InE0(P, p)   // P reduces to the identity component: to O when p divides a denominator, else to a smooth point
        if exists{c : c in Eltseq(P) | Valuation(Denominator(c), p) gt 0} then return true; end if;
        Ep := ChangeRing(E, GF(p)); sing := SingularPoints(Curve(Ep));
        return not (Curve(Ep)!(Ep(GF(p))![GF(p)!c : c in Eltseq(P)]) in sing);
    end function;
    G, m := TorsionSubgroup(E);
    assert Invariants(G) eq [2, 2];
    tors := [ m(g) : g in G | g ne G!0 ];
    assert forall{T : T in tors | exists{p : p in primes | not InE0(T, p)}};   // every nonzero torsion point is seen by some Phi_p

    // (c) the component groups and the specialisation of the branch places
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
    for p in primes do
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
        Sm, Pm := SmithForm(Lap);
        inv := [Sm[i,i] : i in [1..V] | Sm[i,i] notin {0, 1}];
        assert inv eq [TamagawaNumber(E, p)];
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
            assert total eq ratio * DegreeOfFieldOfDefinitionOfCMPoint(star, d);
            for v in Keys(cnt) do assert cnt[v] mod ratio eq 0; cnt[v] := cnt[v] div ratio; end for;
            return cnt;
        end function;
        places := AssociativeArray(); for d in branch do places[d] := Place(d); end for;
        function ClassOf(dv)
            v := ZeroMatrix(Integers(), V, 1);
            for t in dv do for w in Keys(places[t[1]]) do v[w, 1] +:= t[2] * places[t[1]][w]; end for; end for;
            assert &+Eltseq(v) eq 0;
            w := Pm * v;
            return [ w[i,1] mod Sm[i,i] : i in [1..V] | Sm[i,i] notin {0, 1} ];
        end function;
        assert ClassOf([<-340, 1>, <-840, 1>, <-14280, -1>]) eq [0];
    end for;

    // (d) the Prym is simple over Q: its L-polynomial at 11 is irreducible
    R<x> := PolynomialRing(Rationals());
    p := 11;
    s := [ p^i + 1 - ComputePointsViaTrace(Y, p, i) : i in [1..3] ];
    e := [Rationals() | 1];
    for k in [1..3] do Append(~e, (&+[ (-1)^(i-1) * e[k-i+1] * s[i] : i in [1..k] ]) / k); end for;
    Ek := func< k | k le 3 select e[k+1] else p^(k-3) * e[6-k+1] >;
    LY := &+[ (-1)^k * Ek(k) * x^k : k in [0..6] ];
    ok, Lprym := IsDivisibleBy(LY, 1 - TraceOfFrobenius(E, p) * x + p * x^2);
    assert ok and Degree(Lprym) eq 4 and IsIrreducible(Lprym);
    printf " Done!\n";
end procedure;

test_hyperelliptic_714_5();
