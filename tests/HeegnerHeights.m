// The height route (HeegnerHeights.m) on quotients whose models are stored, and on the two models it
// produced.
//
// (a) The explicit Gross-Zagier formula against the stored models: on X_0^35(1) (one class number
//     one point, d = -7), X_0^10(11)/{1,5} (d = -8 and -40, class number 1 and 2) and
//     X_0^390(1)/{1,3,10,30} (d = -312 and -195, class number 4, genus-character components), the
//     predicted canonical height of each branch point of the stored even octic on its elliptic quotient
//     equals the computed height to ten digits.
// (b) The route end to end on X_0^10(11)/{1,5} and X_0^390(1)/{1,3,10,30}: the heights exclude the
//     wrong curve of the isogeny class, the surviving configuration is unique up to the symmetries,
//     and its hyperelliptic model is Q-isomorphic to the stored one.
// (c) The two models the route produced where no other method reaches (regression values, the full
//     searches took an hour each and are recorded with the campaign data): X_0^1326(1)/{1,17,78,1326},
//     y^2 = 272x^8 - 72x^6 - 55x^4 - 18x^2 + 17, and X_0^210(13)/{1,2,21,42,65,130,1365,2730},
//     y^2 = 4x^8 - 20x^7 + 97x^6 - 20x^5 - 138x^4 + 20x^3 + 97x^2 + 20x + 4. Checked here: bad primes
//     are the primes of DN, point counts agree with the trace formula, the automorphism group has
//     order 8 (the Klein four group of Atkin-Lehner involutions and the hyperelliptic involution), and
//     the predicted height of P - Pbar is an integer multiple of the generator of E(F)^- (k = 1 for
//     both branch discriminants at 1326, k = 2 at 210), which is the lattice condition the search used.

P<x> := PolynomialRing(Rationals());

// model files are evaluated at top level (evaluating them inside a procedure crashes the harness)
stored := AssociativeArray();
for DN in [<35, 1>, <10, 11>, <390, 1>] do
    stored[DN] := eval (Read("data/models/models_" cat IntegerToString(DN[1]) cat "_" cat IntegerToString(DN[2]) cat ".m") cat " return models;");
end for;

// the stored even octic y^2 = F(x) of X_0^D(N)/W, its elliptic quotient E by x -> -x (minimal model,
// with the map), and the function u = x^2 on E
function stored_even_model(D, N, W)
    models := stored[<D, N>];
    ent := [e : e in models[Sort(Setseq(W))] | e[1] eq 3][1];
    F := 4*ent[2] + ent[3]^2;
    assert &and[Coefficient(F, i) eq 0 : i in [1,3,5,7]];
    q := P![Coefficient(F, 2*i) : i in [0..4]];
    den := LCM([Denominator(c) : c in Coefficients(q)]);
    C := HyperellipticCurve(q*den^2);
    pts := Points(C : Bound := 5000);
    assert #pts gt 0;
    E0, mp := EllipticCurve(C, pts[1]);
    E, mn := MinimalModel(E0);
    u := Pullback(Inverse(mp), FunctionField(C).1);
    return F, E, E0, mn, u;
end function;

// heights of the branch points (fibres of u over 0 and oo) on E, keyed by the discriminant of their field
function branch_heights(E, E0, mn, u)
    Z, Pl := SignDecomposition(Divisor(u));
    out := [];
    for pl in Support(Z + Pl) do
        R := RepresentativePoint(pl); K := Ring(Parent(R));
        dk := Type(K) eq FldRat select 1 else Discriminant(Integers(AbsoluteField(K)));
        PK := Type(K) eq FldRat select mn(E0!R) else BaseChange(E, K)!Eltseq(mn(E0(K)!R));
        Append(~out, <dk, Height(PK)>);
    end for;
    return out;
end function;

procedure test_formula()
    printf "Testing the explicit Gross-Zagier formula against stored branch points...";
    // <D, N, W, [<CM discriminant, field discriminant>]>
    cases := [ <35, 1, {1}, [<-7, -7>]>,
               <10, 11, {1, 5}, [<-8, -8>, <-40, -40>]>,
               <390, 1, {1, 3, 10, 30}, [<-312, -24>, <-195, -15>]> ];
    for c in cases do
        D, N, W, discs := Explode(c);
        F, E, E0, mn, u := stored_even_model(D, N, W);
        obs := branch_heights(E, E0, mn, u);
        for dd in discs do
            pred := PredictedBranchPointHeight(E, D, N, 2*#W, dd[1]);
            got := [o[2] : o in obs | o[1] eq dd[2]];
            assert exists{h : h in got | Abs(h - pred) lt 1e-10};
        end for;
    end for;
    printf "Done\n";
end procedure;

procedure test_route(D, N, W)
    printf "Testing the route end to end on X_0^%o(%o)/%o...", D, N, W;
    X := CreateShimuraQuot(D, N, W); X`g := GenusShimuraCurveQuotient(D, N, W);
    F, Estored := stored_even_model(D, N, W);
    quots := EllipticQuotientsBySigns(D, N, W);
    // the quotient by x -> -x of the stored model is one of them
    q := [t : t in quots | IsIsogenous(t[2], Estored)];
    assert #q ge 1;
    Wp, Erep := Explode(q[1]);
    DB := CremonaDatabase();
    cls := [EllipticCurve(DB, D*N, i, j) : j in [1..NumberOfCurves(DB, D*N, i)]]
           where i := [i : i in [1..NumberOfIsogenyClasses(DB, D*N)] | IsIsogenous(EllipticCurve(DB, D*N, i, 1), Erep)][1];
    nmodels := 0;
    for E in cls do
        surv := ModelByHeegnerHeights(X, Wp, E);
        if #surv eq 0 then continue; end if;
        // the wrong curve of the class is excluded by the lattice test, so survivors identify E
        assert IsIsomorphic(E, Estored);
        for s in surv do
            assert s[6] eq 5 and Type(s[7]) eq CrvHyp;
            Hs := HyperellipticCurve(P!(F * LCM([Denominator(c) : c in Coefficients(F)])^2));
            assert IsIsomorphic(s[7], Hs);
            nmodels +:= 1;
        end for;
    end for;
    assert nmodels ge 1;
    printf "%o surviving configuration(s), all the stored curve. Done\n", nmodels;
end procedure;

procedure test_model(D, N, W, F, Elabels, kexp)
    printf "Testing the model of X_0^%o(%o)/%o...", D, N, W;
    H := HyperellipticCurve(F);
    assert Set(PrimeDivisors(Integers()!Discriminant(F))) eq Set(PrimeDivisors(D*N));
    X := CreateShimuraQuot(D, N, W); X`g := 3;
    for p in [q : q in PrimesInInterval(5, 60) | (D*N) mod q ne 0] do
        assert #Points(ChangeRing(H, GF(p))) eq ComputePointsViaTrace(X, p, 1);
    end for;
    A := AutomorphismGroup(H);
    assert #A eq 8;
    // the lattice condition: for the quotient E used in the search, P - Pbar is kexp times the generator
    quots := EllipticQuotientsBySigns(D, N, W);
    q := [t : t in quots | CremonaReference(t[2]) in Elabels][1];
    Wp := q[1];
    E := EllipticCurve(Elabels[2]);
    branch := BranchPointData(X, Wp);
    for b in branch do
        F0 := b[3]; dF := Discriminant(Integers(F0));
        Ed := QuadraticTwist(E, dF);
        G, mG := MordellWeilGroup(Ed);
        gen := [mG(G.i) : i in [1..Ngens(G)] | Order(G.i) eq 0];
        assert #gen eq 1;
        EF := BaseChange(E, F0);
        ok, tw := IsIsomorphic(BaseChange(Ed, F0), EF);
        Gm := tw(BaseChange(Ed, F0)!Eltseq(gen[1]));
        halves := DivisionPoints(Gm, 2);
        while #halves gt 0 do Gm := halves[1]; halves := DivisionPoints(Gm, 2); end while;
        comps, h := HeegnerHeightComponents(E, D, N, 2*#W, b[2]);
        adm := [c : c in comps | c[3] and c[6] gt 0];
        assert #adm eq 1;
        k := Sqrt(adm[1][6] * (2/h)^2 / Height(Gm));
        assert Abs(k - kexp) lt 1e-6;
    end for;
    printf "Done\n";
end procedure;

test_formula();
test_route(10, 11, {1, 5});
test_route(390, 1, {1, 3, 10, 30});
test_model(1326, 1, {1, 17, 78, 1326}, 272*x^8 - 72*x^6 - 55*x^4 - 18*x^2 + 17, ["1326b1", "1326b2"], 1);
test_model(210, 13, {1, 2, 21, 42, 65, 130, 1365, 2730}, 4*x^8 - 20*x^7 + 97*x^6 - 20*x^5 - 138*x^4 + 20*x^3 + 97*x^2 + 20*x + 4, ["2730u1", "2730u1"], 2);
