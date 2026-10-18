// The height route for X = X_0^D(N)/W with E = X/W' of rank 1 over Q and both branch points over the same
// quadratic field F (the Klein-four cases).  For a branch point P of CM discriminant d:
//   P - Pbar lies in E(F)^- = E^{dF}(Q), with height (2/h)^2 times that of the admissible character
//   component P^0_chi (Cai-Shu-Tian), so P - Pbar = +-k G' + torsion for a known integer k;
//   P + Pbar = 2 f(xi) + torsion, f(xi) the image of the Hodge class, the SAME for both branch points, an
//   unknown small multiple b of the generator G of E(Q) (the trivial character is inadmissible here, so its
//   Gross-Zagier component is torsion).
// Candidates: halves in E(F) of b G +- k G' + t; configurations with a common b; a mod-p prefilter on the
// point counts of the covers, then the exact test against the trace formula.
AttachSpec("ShimuraQuotients.spec"); SetColumns(0); SetDefaultRealFieldPrecision(20);
P<x> := PolynomialRing(Rationals());
D := StringToInteger(DD); N := StringToInteger(NN); W := {Integers() | StringToInteger(s) : s in Split(WW, ",")};
BMAX := assigned BM select StringToInteger(BM) else 6;
star := func<a, b | a*b div GCD(a, b)^2>;
area := func<E | Abs(Imaginary(ComplexConjugate(Periods(E)[1])*Periods(E)[2]))>;
function Lstar(E)
    L := LSeries(E); v := Evaluate(L, 1);
    if Abs(v) gt 1e-8 then return v, 0; end if;
    v1 := Evaluate(L, 1 : Derivative := 1);
    if Abs(v1) gt 1e-8 then return v1, 1; end if;
    return 0, 2;
end function;
// absolute height of P^0_chi(f) for each genus character chi = (d1,d2), with admissibility
function character_heights(E, Eopt, D, N, nW, d)
    SigD := [p : p in PrimeDivisors(Conductor(E)) | d mod p eq 0];
    eps := AssociativeArray();
    for v in SigD do eps[v] := (D mod v eq 0) select TraceOfFrobenius(E, v) else -TraceOfFrobenius(E, v); end for;
    facs := [<1, d>];
    for d1 in Divisors(Abs(d)) do for s in [1, -1] do
        e1 := s*d1; e2 := d div e1;
        if e1 lt e2 and e1 ne 1 and IsFundamentalDiscriminant(e1) and IsFundamentalDiscriminant(e2) then Append(~facs, <e1, e2>); end if;
    end for; end for;
    u2 := d eq -3 select 9 else (d eq -4 select 4 else 1);
    out := [];
    for f in facs do
        d1 := f[1]; d2 := f[2];
        adm := not exists{v : v in SigD | ((d1 mod v ne 0) select KroneckerSymbol(d1, v) else KroneckerSymbol(d2, v)) ne eps[v]};
        E1 := d1 eq 1 select E else QuadraticTwist(E, d1); E2 := QuadraticTwist(E, d2);
        v1, o1 := Lstar(E1); v2, o2 := Lstar(E2);
        Lp := (o1 + o2 eq 1) select v1*v2 else 0;
        hchi := (1/2) * u2 * Sqrt(Abs(d)) * 2*nW * 2^#SigD * Lp / (2 * ModularDegree(Eopt) * area(Eopt));
        Append(~out, <d1, d2, adm, o1, o2, hchi>);
    end for;
    return out, #facs;
end function;

X := CreateShimuraQuot(D, N, W); X`g := GenusShimuraCurveQuotient(D, N, W);
assert X`g eq 3;
als := [d : d in Divisors(D*N) | GCD(d, (D*N) div d) eq 1 and d notin W];
DB := CremonaDatabase();
found := [];
for u in als do
    Wu := W join {star(u, w) : w in W};
    if GenusShimuraCurveQuotient(D, N, Wu) ne 1 then continue; end if;
    for i in [1..NumberOfIsogenyClasses(DB, D*N)] do
        E := EllipticCurve(DB, D*N, i, 1);
        eps := func<v | (D mod v eq 0) select TraceOfFrobenius(E, v) else -TraceOfFrobenius(E, v)>;
        if &and[&*[Integers() | eps(v) : v in PrimeDivisors(m)] eq 1 : m in Wu] then Append(~found, <Wu, E>); end if;
    end for;
end for;
printf "elliptic quotients consistent with the signs: %o\n", [<Setseq(f[1] diff W), CremonaReference(f[2])> : f in found];
sel := assigned UU select [f : f in found | StringToInteger(UU) in f[1]] else found;
Wp, Eopt := Explode(sel[1]);
cls := [EllipticCurve(DB, D*N, i, j) : j in [1..NumberOfCurves(DB, D*N, i)]] where i := [i : i in [1..NumberOfIsogenyClasses(DB, D*N)] | IsIsogenous(EllipticCurve(DB, D*N, i, 1), Eopt)][1];
printf "W' \\ W = %o, isogeny class %o\n", Setseq(Wp diff W), [CremonaReference(E) : E in cls];

branch := [];
for m in Wp diff W do
    fp := NumFixedPointsByCMOrder(D, N, m);
    for d in Keys(fp) do
        if fp[d] eq 0 then continue; end if;
        K := FieldOfDefinitionOfCMPoint(X, d);
        Append(~branch, <m, d, AbsoluteField(K), fp[d] div #W>);
    end for;
end for;
printf "branch points: %o\n", [<b[1], b[2], Discriminant(Integers(b[3])), b[4]> : b in branch];
require #branch eq 2 and &and[b[4] eq 2 : b in branch]: "two degree-2 branch places expected";
dF := Discriminant(Integers(branch[1][3]));
require Discriminant(Integers(branch[2][3])) eq dF: "both branch points must lie over the same field";
F := branch[1][3];

primes := [p : p in PrimesInInterval(3, 130) | (D*N) mod p ne 0];
target := [ComputePointsViaTrace(X, p, 1) : p in primes];
printf "trace-formula counts: %o\n", target;

// mod-p prefilter: point counts of all double covers of E/F_p branched along the reduction of B
function fp_counts(E, pts, p)
    // pts: points of E(F); returns the set of possible #X(F_p) over the F_p-rational covers, or {} if unusable
    Ep := ChangeRing(E, GF(p));
    OF := Integers(F); pp := Decomposition(OF, p)[1][1]; k, red := ResidueClassField(pp);
    Epk := ChangeRing(Ep, k);
    places := {@ @};
    for Q in pts do
        c := Eltseq(Q);
        if exists{a : a in c | Valuation(a, pp) lt 0} then return {}; end if;
        Include(~places, Place(Epk![red(a) : a in c]));
    end for;
    B := &+[Divisor(pl) : pl in places];
    if Degree(B) ne 4 then return {}; end if;
    O := Divisor(Place(Ep!0));
    D1 := B - 3*O;
    V, mp := RiemannRochSpace(D1);
    if Dimension(V) ne 1 then return {}; end if;
    Deff := D1 + Divisor(mp(V.1));
    supp, mults := Support(Deff);
    PB := Ep!RepresentativePoint(supp[1]);
    AF, m := AlgorithmicFunctionField(FunctionField(Ep));
    counts := {};
    nE := #Ep(GF(p));
    for Q in DivisionPoints(PB, 2) do
        ok, f := IsPrincipal(B - 2*(Divisor(Place(Q)) + O));
        if not ok then continue; end if;
        Fc := FunctionField(Polynomial([AF | -m(f), 0, 1]));
        if Genus(Fc) ne 3 then continue; end if;
        n := NumberOfPlacesOfDegreeOneECF(RationalExtensionRepresentation(Fc));
        counts join:= {n, 2*nE - n};
    end for;
    return counts;
end function;

for E in cls do
    printf "\n==== E = %o (rank %o) ====\n", CremonaReference(E), Rank(E);
    EF := BaseChange(E, F);
    // E(Q) and E(F)^-
    G0, mG0 := MordellWeilGroup(E); gQ := [mG0(G0.i) : i in [1..Ngens(G0)] | Order(G0.i) eq 0];
    Ed := QuadraticTwist(E, dF); Gd, mGd := MordellWeilGroup(Ed); gd := [mGd(Gd.i) : i in [1..Ngens(Gd)] | Order(Gd.i) eq 0];
    require #gQ le 1 and #gd eq 1: "this driver expects rank(E(Q)) <= 1 and rank(E^dF(Q)) = 1";
    ok, tw := IsIsomorphic(BaseChange(Ed, F), EF);
    Gm := tw(BaseChange(Ed, F)!Eltseq(gd[1]));
    halves := DivisionPoints(Gm, 2); while #halves gt 0 do Gm := halves[1]; halves := DivisionPoints(Gm, 2); end while;
    GQ := #gQ eq 1 select EF!Eltseq(gQ[1]) else EF!0;
    if #gQ eq 1 then halves := DivisionPoints(GQ, 2); while #halves gt 0 do GQ := halves[1]; halves := DivisionPoints(GQ, 2); end while; end if;
    T, mT := TorsionSubgroup(EF); tors := [mT(t) : t in T];
    printf "E(F): generator of the minus part height %.6o, generator of E(Q) height %.6o, torsion %o\n", Height(Gm), #gQ eq 1 select Height(GQ) else 0, #tors;
    cands := [* *]; good := true;
    for b in branch do
        d := b[2];
        chs, h := character_heights(E, Eopt, D, N, #W, d);
        printf "disc %o (h = %o): characters %o\n", d, h, [<c[1], c[2], c[3], c[4], c[5], Sprintf("%.6o", c[6])> : c in chs];
        adm := [c : c in chs | c[3] and c[6] gt 0];
        require #adm eq 1: "expected exactly one admissible character with nonzero L-value";
        hminus := adm[1][6] * (2/h)^2;
        k := Sqrt(hminus / Height(Gm));
        printf "   predicted height of P - Pbar: %.6o = k^2 h(G'), k = %.6o\n", hminus, k;
        kk := Round(k);
        if Abs(k - kk) gt 1e-4 then printf "   not a lattice multiple: this E is not the quotient\n"; good := false; break; end if;
        Append(~cands, <d, kk>);
    end for;
    if not good then continue; end if;
    // candidates per b
    nsurv := 0;
    for b in [-BMAX..BMAX] do
        if #gQ eq 0 and b ne 0 then continue; end if;
        lists := [* *];
        for c in cands do
            pts := {@ @};
            for s in [1, -1], t in tors do
                twoP := b*GQ + s*c[2]*Gm + t;
                for Q in DivisionPoints(twoP, 2) do

                    // must be defined over F but not over Q
                    if &and[a in Rationals() : a in Eltseq(Q)] then continue; end if;
                    Include(~pts, E(F)!Eltseq(Q));
                end for;
            end for;
            Append(~lists, pts);
        end for;
        nconf := 0; npass := 0;
        for P1 in lists[1], P2 in lists[2] do
            if Place(P1) eq Place(P2) then continue; end if;
            nconf +:= 1;
            // prefilter at three primes
            rej := false;
            for p in [q : q in primes | q le 60][1..6] do
                S := fp_counts(E, [P1, P2], p);
                if #S eq 0 then continue; end if;
                i := Index(primes, p);
                if target[i] notin S then rej := true; break; end if;
            end for;
            if rej then continue; end if;
            npass +:= 1;
            B := Divisor(Place(P1)) + Divisor(Place(P2));
            if Degree(B) ne 4 then continue; end if;
            covers := DoubleCoversWithBranchDivisor(E, B);
            hits := MatchCoverToCounts(E, covers, target, primes : ExtraPrimes := PrimeDivisors(D*N));
            if #hits eq 0 then continue; end if;
            nsurv +:= 1;
            printf "SURVIVOR b=%o: P1 = %o, P2 = %o; covers noether %o; consistent %o\n", b, P1, P2, [c`noether : c in covers], hits;
            for hh in hits do
                c := TwistedCover(covers[hh[1]], hh[2]);
                if c`noether eq 5 then
                    try
                        H := HyperellipticModelOfCover(c);
                        printf "   HYPERELLIPTIC, model y^2 = %o\n", HyperellipticPolynomials(H);
                    catch err
                        printf "   HYPERELLIPTIC, but %o\n", err`Object;
                    end try;
                else
                    printf "   NOT hyperelliptic\n";
                end if;
            end for;
        end for;
        printf "b = %o: %o configurations, %o passed the prefilter\n", b, nconf, npass;
    end for;
    printf "survivors for E = %o: %o\n", CremonaReference(E), nsurv;
end for;
printf "DONE\n";
