// Cai-Shu-Tian prediction with the TRUE CM discriminants of the branch points (fixed points of the elements of
// W' \ W on X), compared with the observed canonical heights of the branch points on E, matched by field of
// definition. Fundamental CM discriminants with 2-elementary class group only.
AttachSpec("ShimuraQuotients.spec"); SetColumns(0); SetDefaultRealFieldPrecision(20);
P<x> := PolynomialRing(Rationals());
load "/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/heeg/rowsW.m";
star := func<a, b | a*b div GCD(a, b)^2>;
area := func<E | Abs(Imaginary(ComplexConjugate(Periods(E)[1])*Periods(E)[2]))>;
function Lstar(E)
    L := LSeries(E); v := Evaluate(L, 1);
    if Abs(v) gt 1e-8 then return v, 0; end if;
    v1 := Evaluate(L, 1 : Derivative := 1);
    if Abs(v1) gt 1e-8 then return v1, 1; end if;
    return 0, 2;
end function;
function predicted_height(E, Eopt, D, N, nW, d)   // absolute canonical height of f(P), P a CM point of disc d
    SigD := [p : p in PrimeDivisors(Conductor(E)) | d mod p eq 0];
    eps := AssociativeArray();
    for v in SigD do eps[v] := (D mod v eq 0) select TraceOfFrobenius(E, v) else -TraceOfFrobenius(E, v); end for;
    facs := [<1, d>];
    for d1 in Divisors(Abs(d)) do for s in [1, -1] do
        e1 := s*d1; e2 := d div e1;
        if e1 lt e2 and e1 ne 1 and IsFundamentalDiscriminant(e1) and IsFundamentalDiscriminant(e2) then Append(~facs, <e1, e2>); end if;
    end for; end for;
    h := #facs;
    u2 := d eq -3 select 9 else (d eq -4 select 4 else 1);
    tot := 0; used := "";
    for f in facs do
        d1 := f[1]; d2 := f[2]; adm := true;
        for v in SigD do
            chiv := (d1 mod v ne 0) select KroneckerSymbol(d1, v) else KroneckerSymbol(d2, v);
            if chiv ne eps[v] then adm := false; end if;
        end for;
        if not adm then continue; end if;
        E1 := d1 eq 1 select E else QuadraticTwist(E, d1); E2 := QuadraticTwist(E, d2);
        v1, o1 := Lstar(E1); v2, o2 := Lstar(E2);
        if o1 + o2 eq 1 then tot +:= v1*v2; used cat:= Sprintf("(%o,%o)", d1, d2); end if;
    end for;
    return (1/2) * (1/h^2) * u2 * Sqrt(Abs(d)) * 2*nW * 2^#SigD * tot / (2 * ModularDegree(Eopt) * area(Eopt)), used, h;
end function;
nok := 0; nbad := 0;
for c in cases do
    D, N, W, lab := Explode(c); W := {Integers() | w : w in W};
    E := EllipticCurve(lab);
    if Conductor(E) ne D*N then continue; end if;
    Eopt := EllipticCurve(Substring(lab, 1, #lab-1) cat "1");
    eps := func<v | (D mod v eq 0) select TraceOfFrobenius(E, v) else -TraceOfFrobenius(E, v)>;
    epsm := func<m | &*[Integers() | eps(v) : v in PrimeDivisors(m)]>;
    als := [d : d in Divisors(D*N) | GCD(d, (D*N) div d) eq 1 and d notin W];
    Wp := 0;
    for u in als do Wu := W join {star(u, w) : w in W};
        if GenusShimuraCurveQuotient(D, N, Wu) eq 1 and &and[epsm(m) eq 1 : m in Wu] then Wp := Wu; break; end if; end for;
    X := CreateShimuraQuot(D, N, W); X`g := GenusShimuraCurveQuotient(D, N, W);
    // observed: rebuild E from the stored model and compute branch-point heights by field
    models := eval (Read("data/models/models_" cat IntegerToString(D) cat "_" cat IntegerToString(N) cat ".m") cat " return models;");
    ent := [e : e in models[Sort(Setseq(W))] | e[1] eq 3][1];
    F := 4*ent[2] + ent[3]^2; q := P![Coefficient(F, 2*i) : i in [0..4]];
    den := LCM([Denominator(cc) : cc in Coefficients(q)]); C := HyperellipticCurve(q*den^2);
    E0, mp := EllipticCurve(C, Points(C : Bound := 5000)[1]); Em, mn := MinimalModel(E0);
    assert IsIsomorphic(Em, E);
    _, iso := IsIsomorphic(Em, E);
    u := Pullback(Inverse(mp), FunctionField(C).1); Z, Pl := SignDecomposition(Divisor(u));
    obs := [];   // <field disc, height>
    for pl in Support(Z + Pl) do
        R := RepresentativePoint(pl); K := Ring(Parent(R));
        if Type(K) eq FldRat then Append(~obs, <1, Height(iso(mn(E0!R)))>); continue; end if;
        dk := Discriminant(Integers(AbsoluteField(K)));
        PK := BaseChange(E, K)!Eltseq(iso(mn(E0(K)!R)));
        Append(~obs, <dk, Height(PK)>);
    end for;
    printf "\n%o_%o W=%o E=%o W'\\W=%o\n", D, N, Setseq(W), lab, Setseq(Wp diff W);
    for m in Wp diff W do
        fp := NumFixedPointsByCMOrder(D, N, m);
        for d in Keys(fp) do
            if fp[d] eq 0 then continue; end if;
            K := FieldOfDefinitionOfCMPoint(X, d); dk := Type(K) eq FldRat select 1 else Discriminant(Integers(AbsoluteField(K)));
            fund := IsFundamentalDiscriminant(d);
            cl2 := fund and Exponent(ClassGroup(QuadraticField(d))) le 2;
            if cl2 then hp, used, h := predicted_height(E, Eopt, D, N, #W, d); else hp := -1; used := ""; h := 0; end if;
            cand := [o[2] : o in obs | o[1] eq dk];
            match := exists{o : o in cand | Abs(o - hp) lt 1e-6*(1+Abs(hp))};
            if cl2 then if match then nok +:= 1; else nbad +:= 1; end if; end if;
            printf "   w_%o: disc %o (h=%o, fund %o) field disc %o: observed %o  predicted %o  %o %o\n", m, d, h, fund, dk, [Sprintf("%.6o", o) : o in cand], cl2 select Sprintf("%.6o", hp) else "n/a", cl2 select (match select "OK" else "MISMATCH") else "", used;
        end for;
    end for;
end for;
printf "\nmatched %o, mismatched %o\n", nok, nbad;
