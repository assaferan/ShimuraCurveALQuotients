// X_0^6(323)/<w_6, w_19, w_34> IS HYPERELLIPTIC: THE STAR, THE L-FUNCTIONS BEHIND THE TORSION STATEMENT, THE
// SPECIALISATION OF THE BRANCH PLACE, AND WHICH DOUBLE COVER Y IS.  About two minutes.
//
// Y = X_0^6(323)/H, H = <w_6, w_19, w_34>, has genus three and is a double cover of the genus-one star E = X_0^6(323)/W
// branched at ONE place of degree four, the CM points of discriminant -7752 = -4*1938 (the fixed points of w_1938).
// Its four points P, sP, s^2P, s^3P are permuted cyclically by Cl(-7752)/Cl(-7752)[2] = Z/4, and the only
// Galois-stable pairing is {P, s^2P} | {sP, s^3P}; Y is hyperelliptic if and only if
//     Delta = P + s^2P - sP - s^3P
// is zero in E (paper/hyperelliptic-15-38.tex, Section 1).  Delta is defined over the quadratic subfield K of the
// field of the place, one of the quadratic subfields Q(sqrt -2), Q(sqrt 969) of the genus field of Q(sqrt -1938),
// and lies in the minus part of E(K), the Mordell-Weil group of the twist of E by K.  paper/hyperelliptic-6-323.tex
// proves that Delta is TORSION from the Gross-Zagier formula of Yuan-Zhang-Zhang: Delta is the CM cycle twisted by
// the genus character (-8, 969), and L(E/Q(sqrt -1938), chi, s) = L(E^(-8), s) L(E^(969), s) vanishes to order at
// least three at s = 1.  This test checks the inputs of that argument and shows Delta is zero and Y is the
// hyperelliptic cover:
//   (a) E is the curve 1938h1 of Cremona's tables: its traces of Frobenius agree with the star's point counts at every
//       good prime up to the Sturm bound 720 of its level, and its component groups at 2 and 3 have orders 14 and 3,
//       which single it out in its isogeny class (1938h2 has 7 and 6);
//   (b) the class group of discriminant -7752 is Z/8 x Z/2 x Z/2, its ambiguous classes (through which the
//       Atkin-Lehner group acts on the CM points) are the whole 2-torsion, and the genus character (-8, 969) is
//       trivial on them: the four star points form a torsor under Z/4 and Delta is the (-8, 969)-twisted cycle;
//       E^(-8) has root number -1 and E^(969) has root number +1 and rank two, so L(E/k, chi, s) vanishes to
//       order at least three; the minus part of E(K)_tors is {0, T} for both candidate K, T the rational
//       2-torsion point of E, so Delta is 0 or T;
//   (c) T does not reduce to the identity component at 2 (at 3 it does), and at 2 the four points of the place
//       reduce to two components of the quotient Cerednik-Drinfeld graph, two points each, whose difference has
//       order 14 in Phi_2 = Z/14: over K_p (p above 2; ramified for Q(sqrt -2), where thicknesses double and
//       Phi = Z/28; split for Q(sqrt 969)) twice the difference of the two components is neither 0 nor the image
//       of T, so the K-conjugate pair {P, s^2P} cannot share a component, every other distribution of the four
//       points gives Delta = 0 in Phi, and T does not: Delta = 0.  The component groups are the critical groups of
//       the quotient by W of the Cerednik-Drinfeld dual graph (vertices: the left ideal classes of the Eichler
//       order of level 323 in the definite quaternion algebra of discriminant 6/p, in two copies swapped by w_p;
//       singular points: the classes of the level-323p order; thickness |O^*/{+-1}| times the order of the
//       stabiliser in W).  A CM point reduces to the component or singular point whose Eichler order admits the
//       optimal embedding of its CM order (Molina, Ribet bimodules and the specialization of Heegner points,
//       Theorems 5.3, 5.4, 5.8): to a singular point when p ramifies in the CM field; a singular point whose
//       branches are swapped by an involution of W with p | m becomes a smooth point of the merged component;
//   (d) Y is the untwisted cover: the Prym of Y -> E has an irreducible L-polynomial at p = 7, so it is simple
//       over Q, whereas the twist of the hyperelliptic cover by the rational 2-torsion class has a Prym that splits
//       into two elliptic curves over Q (paper/hyperelliptic-6-385.tex).
// Orders of ideal classes are not reproducible between Magma runs, so everything is built from one list.

procedure test_hyperelliptic_6_323()
    printf "X_0^6(323)/<w_6, w_19, w_34>: the star is 1938h1, the twisted L-function vanishes to order three, Delta specialises to zero at 2 for every labelling, the Prym is simple...";
    D := 6; N := 323;
    W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
    H := {1, 6, 19, 34, 51, 114, 646, 969};
    Xstar := CreateShimuraQuot(D, N, W);
    Xstar`g := GenusShimuraCurveQuotient(D, N, W); Xstar`CurveID := 0;
    assert Xstar`g eq 1;
    Y := CreateShimuraQuot(D, N, H); Y`g := GenusShimuraCurveQuotient(D, N, H); Y`CurveID := 0;
    assert Y`g eq 3;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};
    d := -7752;
    assert d eq -4 * D * N;
    assert DegreeOfFieldOfDefinitionOfCMPoint(star, d) eq 4;

    // (a) the isogeny class
    E := MinimalModel(EllipticCurve("1938h1"));
    lev := Conductor(E);
    bound := Floor(2 * lev * &*[Rationals() | 1 + 1/q : q in PrimeDivisors(lev)] / 12);
    assert bound eq 720;
    for p in [q : q in PrimesInInterval(2, bound) | lev mod q ne 0] do
        assert p + 1 - ComputePointsViaTrace(Xstar, p, 1) eq TraceOfFrobenius(E, p);
    end for;
    primes := [2, 3];
    assert [TamagawaNumber(E, p) : p in primes] eq [14, 3];
    assert #[F : F in IsogenousCurves(E) | [TamagawaNumber(F, p) : p in primes] eq [14, 3]] eq 1;

    // (b) the class group, the genus character, the L-functions and the torsion
    Cl, mCl := ClassGroup(QuadraticForms(d));
    assert Invariants(Cl) eq [2, 2, 8];
    amb := sub< Cl | [ g : g in Cl | 2*g eq Cl!0 ] >;
    assert #amb eq 8;                                                   // the 2-torsion, of order 2^(4-1)
    forms := [ Reduction(QuadraticForms(d)![q, 0, (-d div 4) div q]) : q in [2, 3, 17, 19] ];   // the ramified prime ideals
    ramified := [ [ g : g in Cl | Reduction(mCl(g)) eq F ][1] : F in forms ];                   // their classes generate the 2-torsion
    assert sub< Cl | ramified > eq amb;
    assert IsCyclic(quo< Cl | amb >) and #quo< Cl | amb > eq 4;
    assert forall{q : q in [2, 3, 17, 19] | KroneckerSymbol(q eq 2 select 969 else -8, q) eq 1};   // chi_(-8,969) on the ramified classes
    for twist in [<-8, -1, 1>, <969, 1, 2>] do
        Et := MinimalModel(QuadraticTwist(E, twist[1]));
        assert RootNumber(Et) eq twist[2];
        r, proved := Rank(Et); assert proved and r eq twist[3];
    end for;
    G, m := TorsionSubgroup(E);
    assert Invariants(G) eq [2];
    T := m(G.1);
    for dK in [-2, 969] do
        K := QuadraticField(dK); EK := BaseChange(E, K);
        GK, mK := TorsionSubgroup(EK);
        minus := [ g : g in GK | g ne GK!0 and EK![ Conjugate(c) : c in Eltseq(mK(g)) ] eq -mK(g) ];
        assert #minus eq 1 and mK(minus[1]) eq EK!T;    // the minus part of the torsion is {0, T}
    end for;

    // (c) the torsion against the identity components, and the specialisation of the branch place
    function InE0(P, p)   // P reduces to the identity component: to O when p divides a denominator, else to a smooth point
        if exists{c : c in Eltseq(P) | Valuation(Denominator(c), p) gt 0} then return true; end if;
        Ep := ChangeRing(E, GF(p)); sing := SingularPoints(Curve(Ep));
        return not (Curve(Ep)!(Ep(GF(p))![GF(p)!c : c in Eltseq(P)]) in sing);
    end function;
    assert not InE0(T, 2) and InE0(T, 3);
    function OptEmb(dd, Ov)   // optimal embeddings of the order of discriminant dd into the definite order Ov, up to O^*
        t0 := dd mod 2; n := (t0^2 - dd) div 4;
        Zb := ZBasis(Ov);
        Gm := Matrix(Rationals(), 4, 4, [ Trace(Zb[i] * Conjugate(Zb[j])) : i, j in [1..4] ]);
        S := ShortVectors(LatticeWithGram(Gm), 2*n, 2*n);
        d0 := FundamentalDiscriminant(dd); _, f := IsSquare(dd div d0);
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
    function Stab(dd)
        e := dd eq -3 select 3 else (dd eq -4 select 2 else 1);
        if exists{mm : mm in Divisors(D*N) | mm gt 1 and GCD(mm, (D*N) div mm) eq 1 and (dd eq -4*mm or (mm mod 4 eq 3 and dd eq -mm))} then e *:= 2; end if;
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
        // the critical group of the graph whose singular points of thickness t are chains of t - 1 components,
        // for the given thicknesses (over Q_p) and for the doubled ones (over a ramified quadratic extension)
        function Critical(th)
            V := nv; Ed := [];
            for i->e in edges do
                prev := e[1];
                for k in [1..th[i]-1] do V +:= 1; Append(~Ed, <prev, V>); prev := V; end for;
                Append(~Ed, <prev, e[2]>);
            end for;
            Lap := ZeroMatrix(Integers(), V, V);
            for e in Ed do Lap[e[1],e[1]] +:= 1; Lap[e[2],e[2]] +:= 1; Lap[e[1],e[2]] -:= 1; Lap[e[2],e[1]] -:= 1; end for;
            Sm, Pm := SmithForm(Lap);
            inv := [Sm[i,i] : i in [1..V] | Sm[i,i] notin {0, 1}];
            assert #inv eq 1;                                        // cyclic: the graph has one cycle
            cls := func< v | (Pm * v)[[i : i in [1..V] | Sm[i,i] notin {0, 1}][1], 1] mod inv[1] >;
            return inv[1], cls, V;
        end function;
        n0, cls0, V0 := Critical(thick);
        assert n0 eq TamagawaNumber(E, p);
        if p ne 2 then continue; end if;
        // the four star points of discriminant d: all reduce to singular points of X (2 ramifies in Q(sqrt -1938))
        // whose branches are swapped, hence to smooth points of components of E
        cnt := AssociativeArray(); total := 0;
        for e in [1..hp] do
            c := OptEmb(d, RightOrder(LL[e]));
            if c eq 0 then continue; end if;
            assert eorb[e] eq -1;
            v := vorb[org[e]];
            cnt[v] := (IsDefined(cnt, v) select cnt[v] else 0) + c; total +:= c;
        end for;
        ratio := #W div Stab(d);                               // CM points of X above one star point
        assert total eq ratio * 4;
        verts := Sort(Setseq(Keys(cnt)));
        assert #verts eq 2 and forall{v : v in verts | cnt[v] eq 2 * ratio};
        unit := func< v, V | Matrix(Integers(), V, 1, [ i eq v select 1 else 0 : i in [1..V] ]) >;
        // twice the difference of the two components is neither zero nor the image of T, in Phi over K_p
        // for K = Q(sqrt 969) (2 splits: the same graph) and for K = Q(sqrt -2) (2 ramifies: thicknesses double)
        delta0 := cls0(unit(verts[1], V0) - unit(verts[2], V0));
        assert (2 * delta0) mod n0 notin {0, n0 div 2};
        n1, cls1, V1 := Critical([2 * t : t in thick]);
        assert n1 eq 2 * n0;
        K := QuadraticField(-2); OK := Integers(K);
        loc := LocalInformation(BaseChange(E, K), Factorization(2 * OK)[1][1]);
        assert loc[4] eq n1;                                     // the Tamagawa number of E over Q_2(sqrt -2)
        delta1 := cls1(unit(verts[1], V1) - unit(verts[2], V1));
        assert (2 * delta1) mod n1 notin {0, n1 div 2};
    end for;

    // (d) the Prym is simple over Q: its L-polynomial at 7 is irreducible
    R<x> := PolynomialRing(Rationals());
    p := 7;
    s := [ p^i + 1 - ComputePointsViaTrace(Y, p, i) : i in [1..3] ];
    e := [Rationals() | 1];
    for k in [1..3] do Append(~e, (&+[ (-1)^(i-1) * e[k-i+1] * s[i] : i in [1..k] ]) / k); end for;
    Ek := func< k | k le 3 select e[k+1] else p^(k-3) * e[6-k+1] >;
    LY := &+[ (-1)^k * Ek(k) * x^k : k in [0..6] ];
    ok, Lprym := IsDivisibleBy(LY, 1 - TraceOfFrobenius(E, p) * x + p * x^2);
    assert ok and Degree(Lprym) eq 4 and IsIrreducible(Lprym);
    printf " Done!\n";
end procedure;

test_hyperelliptic_6_323();
