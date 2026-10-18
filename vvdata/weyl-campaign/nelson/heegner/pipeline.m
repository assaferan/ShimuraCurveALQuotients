// The height route end to end, on a quotient X = X_0^D(N)/W of genus 3 with an elliptic quotient E = X/W':
//   1. E from the elliptic-curve tables by its Atkin-Lehner signs (a_v at v | D, -a_v at v | N), W' the
//      index-2 overgroup of W with genus 1 on which those signs are +1;
//   2. the branch points of X -> E: the fixed points of the elements of W' \ W, with their CM discriminants
//      and fields of definition F;
//   3. their canonical heights from the explicit Gross-Zagier formula (Cai-Shu-Tian, Theorem 1.5) summed over
//      the admissible characters of the class group;
//   4. the points of E(F) of that height, from the Mordell-Weil groups of E and of its twist;
//   5. every configuration tested exactly: the double covers of E branched there against the trace formula.
// If the stored model is available it is compared with the survivor.
AttachSpec("ShimuraQuotients.spec"); SetColumns(0); SetDefaultRealFieldPrecision(20);
P<x> := PolynomialRing(Rationals());
D := StringToInteger(DD); N := StringToInteger(NN); W := {Integers() | StringToInteger(s) : s in Split(WW, ",")};
star := func<a, b | a*b div GCD(a, b)^2>;
area := func<E | Abs(Imaginary(ComplexConjugate(Periods(E)[1])*Periods(E)[2]))>;
function Lstar(E)
    L := LSeries(E); v := Evaluate(L, 1);
    if Abs(v) gt 1e-8 then return v, 0; end if;
    v1 := Evaluate(L, 1 : Derivative := 1);
    if Abs(v1) gt 1e-8 then return v1, 1; end if;
    return 0, 2;
end function;
function predicted_height(E, Eopt, D, N, nW, d)
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
    tot := 0;
    for f in facs do
        d1 := f[1]; d2 := f[2];
        if exists{v : v in SigD | ((d1 mod v ne 0) select KroneckerSymbol(d1, v) else KroneckerSymbol(d2, v)) ne eps[v]} then continue; end if;
        E1 := d1 eq 1 select E else QuadraticTwist(E, d1); E2 := QuadraticTwist(E, d2);
        v1, o1 := Lstar(E1); v2, o2 := Lstar(E2);
        if o1 + o2 eq 1 then tot +:= v1*v2; end if;
    end for;
    return (1/2) * (1/h^2) * u2 * Sqrt(Abs(d)) * 2*nW * 2^#SigD * tot / (2 * ModularDegree(Eopt) * area(Eopt));
end function;

// 1. E and W'
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
Wp, Eopt := Explode(found[1]);
// The branch points lie on the quotient curve E = X/W', which is SOME curve in the isogeny class; the heights
// depend on it through Area(E). Try every curve in the class.
cls := [EllipticCurve(DB, D*N, i, j) : j in [1..NumberOfCurves(DB, D*N, i)]] where i := [i : i in [1..NumberOfIsogenyClasses(DB, D*N)] | IsIsogenous(EllipticCurve(DB, D*N, i, 1), Eopt)][1];
printf "isogeny class: %o\n", [CremonaReference(E) : E in cls];

// 2. branch points: discriminants and fields
branch := [];   // <m, d, F>
for m in Wp diff W do
    fp := NumFixedPointsByCMOrder(D, N, m);
    for d in Keys(fp) do
        if fp[d] eq 0 then continue; end if;
        K := FieldOfDefinitionOfCMPoint(X, d);
        Append(~branch, <m, d, K, fp[d] div #W>);
    end for;
end for;
printf "branch points: %o\n", [<b[1], b[2], Type(b[3]) eq FldRat select 1 else Discriminant(Integers(AbsoluteField(b[3]))), b[4]> : b in branch];
require #branch eq 2 and &and[b[4] eq 2 : b in branch]: "this driver handles two degree-2 branch places";

// trace-formula target
primes := [p : p in PrimesInInterval(3, 130) | (D*N) mod p ne 0];
target := [ComputePointsViaTrace(X, p, 1) : p in primes];

for E in cls do
    printf "\n==== E = %o (rank %o) ====\n", CremonaReference(E), Rank(E);
    if Rank(E) gt 0 then printf "rank > 0: the trivial-character component needs the Hodge-class shift; skipped\n"; continue; end if;
    cands := [* *];
    for b in branch do
        d := b[2]; K := AbsoluteField(b[3]); dK := Discriminant(Integers(K));
        hpred := predicted_height(E, Eopt, D, N, #W, d);
        // 4. E(K) = E(Q) + E^dK(Q) up to 2-torsion; E(Q) is torsion here
        EK := BaseChange(E, K);
        Ed := QuadraticTwist(E, dK);
        G, mG := MordellWeilGroup(Ed);
        gens := [mG(G.i) : i in [1..Ngens(G)] | Order(G.i) eq 0];
        ok, tw := IsIsomorphic(BaseChange(Ed, K), EK);
        gensK := [tw(BaseChange(Ed, K)!Eltseq(g)) : g in gens];
        // saturate at 2 inside E(K)
        sat := [];
        for g in gensK do
            g2 := g; halves := DivisionPoints(g2, 2);
            while #halves gt 0 do g2 := halves[1]; halves := DivisionPoints(g2, 2); end while;
            Append(~sat, g2);
        end for;
        T, mT := TorsionSubgroup(EK);
        tors := [mT(t) : t in T];
        printf "disc %o over Q(sqrt %o): predicted height %.6o; E^d rank %o, generator heights %o, torsion %o\n", d, SquarefreeFactorization(dK), hpred, #sat, [Sprintf("%.6o", Height(g)) : g in sat], #tors;
        require #sat eq 1: "rank of E(K) is not 1; the enumeration below assumes a single generator";
        g := sat[1]; hg := Height(g);
        k := Sqrt(hpred/hg);
        printf "   height ratio gives multiple k = %.6o\n", k;
        kk := Round(k);
        if Abs(k - kk) gt 1e-4 then
            printf "   not a lattice point on this curve of the isogeny class; this E is not the quotient\n";
            cands := [* *]; break;
        end if;
        pts := [E(K)!Eltseq(s*kk*g + t) : s in [1, -1], t in tors];
        Append(~cands, <d, pts>);
    end for;
    if #cands ne 2 then continue; end if;
    // 5. configurations
    places := [* [Place(p) : p in c[2]] : c in cands *];
    // distinct places only (a point and its conjugate give the same place)
    uniq := [SetToIndexedSet({@ pl : pl in L @}) : L in places];
    printf "candidate places: %o x %o\n", #uniq[1], #uniq[2];
    survivors := [];
    for i in [1..#uniq[1]], j in [1..#uniq[2]] do
        B := Divisor(uniq[1][i]) + Divisor(uniq[2][j]);
        if Degree(B) ne 4 then continue; end if;
        covers := DoubleCoversWithBranchDivisor(E, B);
        hits := MatchCoverToCounts(E, covers, target, primes : ExtraPrimes := PrimeDivisors(D*N));
        if #hits gt 0 then
            Append(~survivors, <i, j, hits, [c`noether : c in covers]>);
            printf "   configuration (%o,%o): covers noether %o, consistent: %o\n", i, j, [c`noether : c in covers], hits;
            for h in hits do
                c := TwistedCover(covers[h[1]], h[2]);
                if c`noether eq 5 then
                    H := HyperellipticModelOfCover(c);
                    printf "      hyperelliptic model: y^2 = %o\n", HyperellipticPolynomials(H);
                    if assigned MODELS then
                        models := eval (Read(MODELS) cat " return models;");
                        ent := [e : e in models[Sort(Setseq(W))] | e[1] eq 3][1];
                        Hs := HyperellipticCurve(P!((4*ent[2] + ent[3]^2) * LCM([Denominator(cc) : cc in Coefficients(4*ent[2] + ent[3]^2)])^2));
                        printf "      Q-isomorphic to the stored model: %o\n", IsIsomorphic(H, Hs);
                    end if;
                else
                    printf "      not hyperelliptic\n";
                end if;
            end for;
        end if;
    end for;
    printf "surviving configurations: %o\n", #survivors;
end for;
printf "DONE\n";
