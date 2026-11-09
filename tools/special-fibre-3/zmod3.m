// The genus-3 quotients Z = X_0^D(N)/W0 below (and the three X = X_0^D(3N)/W whose special fibre at 3 is Z),
// each a double cover of its elliptic star E: the configurations of branch points allowed by the certified
// Heegner heights (ModelByHeegnerHeights), and for every double cover of E with such a branch divisor its
// reduction mod p, its genus and the span of the products of its differentials (5 hyperelliptic, 6 not), with
// its L-polynomial, which for the surviving cover must be the trace formula's.  For 14_17 and 38_5 the
// surviving cover is not hyperelliptic over Q but is mod 3, so the special fibre of 14_51 and 38_15 at 3 does
// not rule them out; on X (OLD:=3, the star's form old at 3) the surviving cover is hyperelliptic over Q.
//   magma -b LAB:=14_17 EL:=238a2 tools/special-fibre-3/zmod3.m < /dev/null
//   magma -b LAB:=38_5 EL:=190a1 tools/special-fibre-3/zmod3.m < /dev/null
//   magma -b LAB:=14_51 EL:=238a1 OLD:=3 tools/special-fibre-3/zmod3.m < /dev/null
//   magma -b LAB:=38_15 EL:=190a1 OLD:=3 tools/special-fibre-3/zmod3.m < /dev/null
AttachSpec("ShimuraQuotients.spec");
SetVerbose("HeegnerHeights", 1);
W0s := AssociativeArray(); W0s["14_17"] := [1,14,34,119]; W0s["38_5"] := [1,10,38,95]; W0s["158_1"] := [1,2]; W0s["14_51"] := [1,6,14,21,34,51,119,714]; W0s["38_15"] := [1,3,10,30,38,95,114,285]; W0s["158_3"] := [1,2,237,474];
D := StringToInteger(Split(LAB, "_")[1]); N := StringToInteger(Split(LAB, "_")[2]);
p := assigned PP select StringToInteger(PP) else 3;
H := {Integers() | m : m in W0s[LAB]}; DN := D*N;
full := {Integers() | m : m in Divisors(DN) | GCD(m, DN div m) eq 1};
Z := CreateShimuraQuot(D, N, H); Z`g := GenusShimuraCurveQuotient(D, N, H); Z`CurveID := 1;
E := MinimalModel(EllipticCurve(EL));
function noether_mod_p(E, f, p, kappa)
    Ep := ChangeRing(E, GF(p)); FEp := FunctionField(Ep);
    // numerator and denominator made primitive separately; f mod p is then p^v times their quotient,
    // and the cover y^2 = f has good reduction only when v is even (y rescaled by p^(v/2))
    pf := ProjectiveFunction(f);
    Rq := Parent(Numerator(pf)); R3Q := Type(Rq) eq RngMPolRes select PreimageRing(Rq) else Rq;
    num := R3Q!Numerator(pf); den := R3Q!Denominator(pf);
    prim := func< P | LCM([Denominator(c) : c in Coefficients(P)]) / GCD([Numerator(c) : c in Coefficients(P)]) >;
    sn := prim(num); sd := prim(den);
    v := Valuation(sd, p) - Valuation(sn, p);
    R3 := PolynomialRing(GF(p), 3);
    nump := R3!(sn*num); denp := R3!(sd*den);
    eqp := R3!DefiningPolynomial(Ep);
    if IsDivisibleBy(nump, eqp) or IsDivisibleBy(denp, eqp) then return -4, v; end if;
    if IsOdd(v) then return -5, v; end if;
    u := kappa * (sd/sn) / p^v;
    fp := (GF(p)!u) * Evaluate(nump, [FEp.1, FEp.2, 1]) / Evaluate(denp, [FEp.1, FEp.2, 1]);
    AFp, mp := AlgorithmicFunctionField(FEp);
    try
        Fp := FunctionField(Polynomial([AFp | -mp(fp), 0, 1]));
        g := Genus(Fp);
    catch err
        return -2, -2;
    end try;
    if g ne 3 then return g, -3; end if;
    diffs := BasisOfDifferentialsFirstKind(Fp);
    gs := [d / diffs[1] : d in diffs];
    prods := [gs[i]*gs[j] : j in [i..#gs], i in [1..#gs]];
    // coordinates down the tower, over GF(p)(x), cleared of denominators
    coords := [[c : c in Eltseq(a), a in Eltseq(e)] : e in prods];
    den := LCM([Denominator(c) : c in cs, cs in coords]);
    polys := [[Numerator(c*den) : c in cs] : cs in coords];
    deg := Max([Degree(q) : q in ps, ps in polys] cat [0]);
    rows := [&cat[[Coefficient(q, i) : i in [0..deg]] : q in ps] : ps in polys];
    M := Matrix(GF(p), #prods, #polys[1]*(deg+1), rows);
    printf "   L-polynomial mod %o: %o\n", p, LPolynomial(Fp);
    return g, Rank(M);
end function;
if N mod p ne 0 then printf "trace formula for Z at %o: %o\n", p, ZetaNumeratorViaTrace(Z, p : Check := false); end if;
t0 := Realtime();
old := assigned OLD select StringToInteger(OLD) else 0;
res := ModelByHeegnerHeights(Z, full, E : Certified := true, OldformPrime := old);
printf "%o on %o: %o survivor(s) (%o s)\n", LAB, EL, #res, Realtime() - t0;
for r in res do
    if Type(r[6]) eq CrvHyp then
        H := SmallHyperellipticModel(r[6]);
        printf "MODEL %o: %o, conductor-free check: genus %o, discriminant %o\n", LAB, DefiningPolynomial(H), Genus(H), Factorization(Integers()!Discriminant(H));
    end if;
end for;
done := {};
for r in res do
    B := &+[Divisor(Place(P)) : P in r[2]];
    if B in done then continue; end if; Include(~done, B);
    covers := DoubleCoversWithBranchDivisor(E, B);
    for i in [1..#covers] do
        C := covers[i];
        g, nq := noether_mod_p(E, C`f, p, i eq r[3] select r[4] else 1);
        printf "RESULT %o: config %o, cover %o of %o (survivor uses %o, kappa %o; Noether over Q %o): mod %o genus %o, span %o\n",
            LAB, r[1], i, #covers, r[3], r[4], C`noether, p, g, nq;
    end for;
end for;
quit;
