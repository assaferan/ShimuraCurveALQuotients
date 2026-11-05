// Models of genus-3 Atkin-Lehner quotients from the heights of their branch points.
//
// Let X = X_0^D(N)/W have genus 3 and an elliptic quotient E = X/W', W' = <W, w> of index 2. Then X -> E
// is a double cover branched at the fixed points of the elements of W' \ W, which are CM points of known
// discriminants over known quadratic fields F. The explicit Gross-Zagier formula of Cai, Shu and Tian
// (Algebra & Number Theory 8 (2014), Theorem 1.5) gives the canonical heights of their images on E in
// terms of L-values, up to one component that only depends on the choice of origin; the points of E(F)
// of those heights are a finite set; and each configuration of four such points is tested exactly by
// ExactDoubleCover.m against the trace formula. Nothing here is numerical except the L-values, which
// enter only through the recognition of an integer.
//
// Conventions. The Atkin-Lehner sign of E at a prime v is a_v(E) for v | D and -a_v(E) for v | N (the
// sign of the quaternionic newform). A character chi = (d1, d2) of the class group of a fundamental
// discriminant d = d1 d2 with 2-elementary class group (a genus character) is admissible when, at every
// prime v dividing both the conductor of E and d, chi at the ramified prime above v equals that sign;
// otherwise the chi-component of the Heegner divisor is torsion. Heights are the absolute canonical
// heights of Magma; the height over the quadratic field used by Gross-Zagier is twice that.

declare verbose HeegnerHeights, 2;

al_star := func<a, b | a*b div GCD(a, b)^2>;

function lattice_area(E)
    per := Periods(E);
    return Abs(Imaginary(ComplexConjugate(per[1])*per[2]));
end function;

// L(E,1) if nonzero, else L'(E,1); the second value is the order used (2 means neither, within precision).
function L_value_or_derivative(E)
    L := LSeries(E); v := Evaluate(L, 1);
    if Abs(v) gt 1e-8 then return v, 0; end if;
    v1 := Evaluate(L, 1 : Derivative := 1);
    if Abs(v1) gt 1e-8 then return v1, 1; end if;
    return 0, 2;
end function;

intrinsic QuaternionicAtkinLehnerSign(E::CrvEll, v::RngIntElt, D::RngIntElt, N::RngIntElt) -> RngIntElt
{The eigenvalue of the Atkin-Lehner involution w_v on the newform of X_0^D(N) attached to E, for a
prime v dividing DN: a_v(E) when v divides D and -a_v(E) when v divides N.}
    require IsPrime(v) and (D*N) mod v eq 0: "v must be a prime dividing DN";
    return (D mod v eq 0) select TraceOfFrobenius(E, v) else -TraceOfFrobenius(E, v);
end intrinsic;

intrinsic EllipticQuotientsBySigns(D::RngIntElt, N::RngIntElt, W::SetEnum) -> SeqEnum
{The pairs <W', E> with W' an index-2 overgroup of W in the Atkin-Lehner group such that X_0^D(N)/W' has
genus 1, and E the first curve of the isogeny class of conductor DN whose Atkin-Lehner signs are +1 on
W'. Only newforms of level exactly DN are considered, so quotients whose Jacobian has lower conductor
are not listed.}
    W := {Integers() | w : w in W};
    DB := CremonaDatabase();
    out := [];
    for u in [d : d in Divisors(D*N) | GCD(d, (D*N) div d) eq 1 and d notin W] do
        Wu := W join {al_star(u, w) : w in W};
        if exists{o : o in out | o[1] eq Wu} then continue; end if;
        if GenusShimuraCurveQuotient(D, N, Wu) ne 1 then continue; end if;
        for i in [1..NumberOfIsogenyClasses(DB, D*N)] do
            E := EllipticCurve(DB, D*N, i, 1);
            if &and[&*[Integers() | QuaternionicAtkinLehnerSign(E, v, D, N) : v in PrimeDivisors(m)] eq 1 : m in Wu] then
                Append(~out, <Wu, E>);
            end if;
        end for;
    end for;
    return out;
end intrinsic;

intrinsic BranchPointData(X::ShimuraQuot, Wp::SetEnum) -> List
{For the double cover X = X_0^D(N)/W -> X/W', the branch points grouped by CM order: a list of
tuples <m, d, F, n> with m in W' \ W, d the discriminant of the CM order of the fixed points of
w_m, F the field of definition of those points on X (the rationals or an absolute number field),
and n their number on X.}
    D := X`D; N := X`N; W := {Integers() | w : w in X`W};
    out := [* *];
    for m in Wp diff W do
        fp := NumFixedPointsByCMOrder(D, N, m);
        for d in Keys(fp) do
            if fp[d] eq 0 then continue; end if;
            K := FieldOfDefinitionOfCMPoint(X, d);
            Append(~out, <m, d, Type(K) eq FldRat select K else AbsoluteField(K), fp[d] div #W>);
        end for;
    end for;
    return out;
end intrinsic;

intrinsic GenusCharacters(d::RngIntElt) -> SeqEnum
{The factorisations d = d1 d2 into fundamental discriminants with d1 < d2, the trivial character
<1, d> first: the genus characters, i.e. the characters of order dividing 2 of the class group of d.
For a 2-elementary class group these are all its characters.}
    require IsFundamentalDiscriminant(d): "d must be a fundamental discriminant";
    facs := [<1, d>];
    for d1 in Divisors(Abs(d)), s in [1, -1] do
        e1 := s*d1; e2 := d div e1;
        if e1 lt e2 and e1 ne 1 and IsFundamentalDiscriminant(e1) and IsFundamentalDiscriminant(e2) then
            Append(~facs, <e1, e2>);
        end if;
    end for;
    return facs;
end intrinsic;

intrinsic HeegnerHeightComponents(E::CrvEll, D::RngIntElt, N::RngIntElt, degf::RngIntElt, d::RngIntElt) -> SeqEnum, RngIntElt
{For a CM point P of fundamental discriminant d on X_0^D(N) and f : X_0^D(N) -> E of degree degf, the
genus-character components of the Heegner divisor of f(P): for each genus character chi = (d1, d2) a
tuple <d1, d2, admissible, ord L(E^d1,1), ord L(E^d2,1), h_chi>, where h_chi is the absolute canonical
height of sum_t chi(t) f(P)^t over the whole class group when chi is admissible, by Cai-Shu-Tian's
formula with 8 pi^2 (phi,phi) = 2 deg(pi_opt) Area(E_opt). The second value is the class number.
Characters of higher order are not computed: a point defined over a quadratic field F has nonzero
components only for the characters trivial on Gal(H/F), the trivial one and the genus character
attached to F.}
    require IsFundamentalDiscriminant(d) and d lt 0: "d must be a negative fundamental discriminant";
    // for a squarefree conductor no Euler factor is removed from L(s, E, chi) in the formula, and the
    // Manin constant of the optimal curve is 1
    require IsSquarefree(Conductor(E)): "the conductor of E must be squarefree";
    Cl := ClassGroup(QuadraticField(d));
    // the first curve of its class in Cremona's table is the optimal one, except in class 990h
    require Conductor(E) ne 990: "the optimal curve of class 990h is not the first of its class";
    lab := CremonaReference(E);
    Eopt := EllipticCurve(Substring(lab, 1, #lab - 1) cat "1");
    SigD := [p : p in PrimeDivisors(Conductor(E)) | d mod p eq 0];
    eps := AssociativeArray();
    for v in SigD do eps[v] := QuaternionicAtkinLehnerSign(E, v, D, N); end for;
    u2 := d eq -3 select 9 else (d eq -4 select 4 else 1);
    norm := u2 * Sqrt(RealField(30)!Abs(d)) * degf * 2^#SigD / (2 * ModularDegree(Eopt) * lattice_area(Eopt));
    out := [];
    for f in GenusCharacters(d) do
        d1 := f[1]; d2 := f[2];
        adm := not exists{v : v in SigD | ((d1 mod v ne 0) select KroneckerSymbol(d1, v) else KroneckerSymbol(d2, v)) ne eps[v]};
        E1 := d1 eq 1 select E else QuadraticTwist(E, d1); E2 := QuadraticTwist(E, d2);
        v1, o1 := L_value_or_derivative(E1); v2, o2 := L_value_or_derivative(E2);
        Lp := (o1 + o2 eq 1) select v1*v2 else 0;
        Append(~out, <d1, d2, adm, o1, o2, (1/2) * norm * Lp>);
    end for;
    return out, #Cl;
end intrinsic;

intrinsic PredictedBranchPointHeight(E::CrvEll, D::RngIntElt, N::RngIntElt, degf::RngIntElt, d::RngIntElt) -> FldReElt
{The absolute canonical height of f(P) for a CM point P of discriminant d, when the trivial-character
component is given by the formula (E of rank 0 or the trivial character inadmissible) and the class
group is 2-elementary, so that the genus characters are all the characters: the sum of the admissible
components divided by the square of the class number.}
    comps, h := HeegnerHeightComponents(E, D, N, degf, d);
    require Exponent(ClassGroup(QuadraticField(d))) le 2: "the class group must be 2-elementary";
    return &+[RealField(30) | c[6] : c in comps | c[3]] / h^2;
end intrinsic;

// Context for counting double covers of E over F_p: the reduction of E and its function field.
function reduction_context(E, p)
    Ep := ChangeRing(E, GF(p));
    AFp, mp := AlgorithmicFunctionField(FunctionField(Ep));
    return rec< recformat< p, Ep, AFp, mp, nE > | p := p, Ep := Ep, AFp := AFp, mp := mp, nE := #Ep(GF(p)) >;
end function;

// The set of #X(F_p) over the double covers of E/F_p branched along the reduction of the points pts
// (each over a number field), with their twists; empty when the reduction is unusable.
function counts_mod_p(ctx, pts)
    places := {@ @};
    for Q in pts do
        c := Eltseq(Q); F := Ring(Parent(Q));
        if Type(F) eq FldRat then
            if exists{a : a in c | Valuation(a, ctx`p) lt 0} then return {}; end if;
            Include(~places, Place(ctx`Ep![GF(ctx`p) | a : a in c]));
            continue;
        end if;
        pp := Decomposition(Integers(F), ctx`p)[1][1];
        if exists{a : a in c | Valuation(a, pp) lt 0} then return {}; end if;
        k, red := ResidueClassField(pp);
        Include(~places, Place(ctx`Ep(k)![red(a) : a in c]));
    end for;
    B := &+[Divisor(pl) : pl in places];
    if Degree(B) ne 4 then return {}; end if;
    O := Divisor(Place(ctx`Ep!0));
    V, mV := RiemannRochSpace(B - 3*O);
    if Dimension(V) ne 1 then return {}; end if;
    supp := Support(B - 3*O + Divisor(mV(V.1)));
    PB := ctx`Ep!RepresentativePoint(supp[1]);
    counts := {};
    for Q in DivisionPoints(PB, 2) do
        ok, f := IsPrincipal(B - 2*(Divisor(Place(Q)) + O));
        if not ok then continue; end if;
        // two branch points meeting modulo p make f a square; then there is no information
        try
            Fc := FunctionField(Polynomial([ctx`AFp | -ctx`mp(f), 0, 1]));
            g := Genus(Fc);
        catch err
            return {};
        end try;
        if g ne 3 then continue; end if;
        n := NumberOfPlacesOfDegreeOneECF(RationalExtensionRepresentation(Fc));
        counts join:= {n, 2*ctx`nE - n};
    end for;
    return counts;
end function;

intrinsic ModelByHeegnerHeights(X::ShimuraQuot, Wp::SetEnum, E::CrvEll : Primes := 130, TwistPrimes := 800, Effort := 1) -> List
{The configurations of branch points of X = X_0^D(N)/W -> E = X/W' allowed by the Heegner heights
whose double cover is consistent with the trace formula. E must be a curve of the isogeny class of
X/W' (every curve of the class should be tried: on the wrong one the predicted heights are not
heights of lattice points, and the list is empty). The branch places must have degree 1 or 2, E must
have rank 0 or 1 over Q, and its twist by the field of a place of degree 2 must have rank 1. Each
entry is a tuple <ns, pts, cover, kappa, noether, H>: ns the multiples of the generator of E(Q) in
the sums of the points over the branch places, pts one point per place, cover and kappa the index and
the twist identified by MatchCoverToCounts, noether 5 or 6, and H the hyperelliptic model when
noether is 5 (0 when the canonical conic has no rational point). Configurations that are translates of
an earlier one by a rational point of E, or their negatives, are omitted. Effort is passed to
MordellWeilGroup, which must prove the rank and the saturation of E and of the twists; raise it when
it cannot.}
    D := X`D; N := X`N; W := {Integers() | w : w in X`W};
    require assigned X`g and X`g eq 3: "X must have genus 3 (assign X`g)";
    branch := BranchPointData(X, Wp);
    require &and[b[4] in {1, 2} : b in branch] and &+[b[4] : b in branch] eq 4:
        "branch places of degree 1 or 2, of total degree 4, are expected";
    degf := 2*#W;
    G0, mG0, rk0, sat0 := MordellWeilGroup(E : Effort := Effort);
    require rk0 and sat0: "the Mordell-Weil group of E is not proved";
    gQ := [mG0(G0.i) : i in [1..Ngens(G0)] | Order(G0.i) eq 0];
    require #gQ le 1: "E must have rank at most 1";
    hG := #gQ eq 1 select Height(gQ[1]) else 0;
    // Per branch place: its degree e, its field F (Q when e = 1), the group E(F), which is
    // E(Q) + E^dF(Q) up to 2-torsion with generators saturated at 2, the multiple k of the generator
    // of E(F)^- in P - Pbar, and the multiple m of the generator of E(Q) in the trivial character's
    // component of the Heegner point.
    data := [* *];
    for b in branch do
        e := b[4]; F := e eq 2 select b[3] else Rationals();
        EF := e eq 2 select BaseChange(E, F) else E;
        Gm := EF!0; k := 0;
        comps, h := HeegnerHeightComponents(E, D, N, degf, b[2]);
        if e eq 2 then
            dF := Discriminant(Integers(F));
            Ed := QuadraticTwist(E, dF);
            Gd, mGd, rkd, satd := MordellWeilGroup(Ed : Effort := Effort);
            require rkd and satd: "the Mordell-Weil group of the twist of E is not proved";
            gd := [mGd(Gd.i) : i in [1..Ngens(Gd)] | Order(Gd.i) eq 0];
            require #gd eq 1: "the twist of E by the field of a branch point must have rank 1";
            _, tw := IsIsomorphic(BaseChange(Ed, F), EF);
            Gm := tw(BaseChange(Ed, F)!Eltseq(gd[1]));
            halves := DivisionPoints(Gm, 2);
            while #halves gt 0 do Gm := halves[1]; halves := DivisionPoints(Gm, 2); end while;
            // P - Pbar is (2/h) times the component of the character attached to the field F of the
            // branch point: the trivial character when F is the CM field itself (then Gal(H/K) fixes P
            // and the other components vanish), otherwise the genus character (d1, d2) with Q(sqrt d1)
            // or Q(sqrt d2) equal to F; characters attached to other fields vanish on the orbit of P.
            // That character must be admissible with a nonzero value, or P would be rational up to
            // torsion.
            adm := IsSquare(dF*b[2]) select [c : c in comps | c[1] eq 1]
                   else [c : c in comps | c[1] ne 1 and (IsSquare(c[1]*dF) or IsSquare(c[2]*dF))];
            require #adm eq 1: "no genus character cuts out the field of the branch point";
            require adm[1][3] and adm[1][6] gt 0:
                "the character attached to the field of the branch point is inadmissible or has vanishing value";
            hminus := adm[1][6] * (2/h)^2;
            kr := Sqrt(hminus / Height(Gm));
            vprintf HeegnerHeights: "disc %o over Q(sqrt %o): P - Pbar has predicted height %o = k^2 h(G'), k = %o\n",
                b[2], SquarefreeFactorization(dF), hminus, kr;
            if Abs(kr - Round(kr)) gt 1e-4 then
                vprintf HeegnerHeights: "not a lattice point: E is not the quotient\n";
                return [* *];
            end if;
            k := Round(kr);
        else
            // A branch place of degree 1: every conjugate of the CM point has the same image in E, so
            // the component of every nontrivial character vanishes.
            require &and[c[6] lt 1e-10 : c in comps | c[1] ne 1 and c[3]]:
                "a branch place of degree 1 whose nontrivial characters do not all have vanishing components";
            vprintf HeegnerHeights: "disc %o: the branch place is rational\n", b[2];
        end if;
        // The trivial character's component is (h/e) times the sum of the points over the place minus
        // e*xi_E, with xi_E the image in E of the Hodge class of X_0^D(N). The sum lies in E(Q) and
        // xi_E in E(Q) tensor Q, so that component's height gives |n - e*x| with n the multiple of the
        // generator of E(Q) in the sum and x the multiple in xi_E: a rational number, whose denominator
        // is that of x and need not be 1.
        triv := [c : c in comps | c[1] eq 1][1];
        bp := Rationals()!0;
        if #gQ eq 1 and triv[3] and triv[6] gt 0 then
            br := (e/h) * Sqrt(triv[6] / hG);
            bp := BestApproximation(br, 10^4);
            vprintf HeegnerHeights: "disc %o: the sum of the points over the place is e*xi_E +- %o*G, measured %o\n",
                b[2], bp, br;
            if Abs(br - bp) gt 1e-15 then
                vprintf HeegnerHeights: "not a rational multiple of the generator: E is not the quotient\n";
                return [* *];
            end if;
        end if;
        T, mT := TorsionSubgroup(EF);
        Append(~data, rec< recformat< e, F, EF, Gm, k, bp, h, tors, d > | e := e, F := F, EF := EF,
            Gm := Gm, k := k, bp := bp, h := h, tors := [mT(t) : t in T], d := b[2] >);
    end for;
    // The sum of the points over place i is n_i*G + torsion with n_i = e_i*y + s_i*b_i, where y is the
    // multiple of G in xi_E: it is common to all the places and matters only modulo 1, because
    // translating every branch point by a rational point of E gives an isomorphic cover. Each sign s_i
    // is free; changing all of them at once is [-1] on E, whose configurations are removed as
    // duplicates below. Any y making every n_i an integer has denominator dividing the bound below, so
    // the enumeration over y is finite and complete.
    nvecs := {};
    L := LCM([Integers() | 2*Denominator(dt`bp) : dt in data]);
    for sg in CartesianPower([1, -1], #data) do
        if #gQ eq 0 then Include(~nvecs, [0 : i in [1..#data]]); continue; end if;
        for u in [0..L-1] do
            ns := [Rationals() | data[i]`e * u/L + sg[i] * data[i]`bp : i in [1..#data]];
            if &and[IsIntegral(n) : n in ns] then Include(~nvecs, [Integers() | n : n in ns]); end if;
        end for;
    end for;
    if #nvecs eq 0 then
        vprintf HeegnerHeights: "the rational parts of the branch places are not compatible: E is not the quotient\n";
        return [* *];
    end if;
    // trace-formula counts
    primes := [p : p in PrimesInInterval(3, Primes) | (D*N) mod p ne 0];
    target := [ComputePointsViaTrace(X, p, 1) : p in primes];
    small := [q : q in primes | q le 60];
    ctxs := [reduction_context(E, p) : p in small[1..Min(6, #small)]];
    out := [* *];
    seen := [* *];
    for ns in Sort(SetToSequence(nvecs)) do
        lists := [* *];
        for i in [1..#data] do
            dt := data[i]; pts := {@ @};
            nGF := #gQ eq 1 select dt`EF!Eltseq(ns[i]*gQ[1]) else dt`EF!0;
            if dt`e eq 1 then
                for t in dt`tors do Include(~pts, E(dt`F)!Eltseq(nGF + t)); end for;
            else
                for sg in [1, -1], t in dt`tors do
                    for Q in DivisionPoints(nGF + sg*dt`k*dt`Gm + t, 2) do
                        if &and[a in Rationals() : a in Eltseq(Q)] then continue; end if;
                        Include(~pts, E(dt`F)!Eltseq(Q));
                    end for;
                end for;
            end if;
            Append(~lists, pts);
        end for;
        vprintf HeegnerHeights: "rational parts %o: %o candidate points\n", ns, [#l : l in lists];
        for idx in CartesianProduct([[1..#l] : l in lists]) do
            cfg := [* lists[i][idx[i]] : i in [1..#data] *];
            // the branch divisor must be reduced of degree 4: two of the chosen points may lie over
            // the same place
            B := &+[Divisor(Place(P)) : P in cfg];
            _, mlt := Support(B);
            if Degree(B) ne 4 or not &and[mm eq 1 : mm in mlt] then continue; end if;
            // translates by a rational point, or negatives, of an earlier configuration give an
            // isomorphic cover (translations and [-1] are automorphisms of E over Q)
            dup := false;
            for sn in seen, sg in [1, -1] do
                Tc := sg*cfg[1] - sn[1];
                if not &and[a in Rationals() : a in Eltseq(Tc)] then continue; end if;
                Tq := E(Rationals())!Eltseq(Tc);
                if &and[Place(sg*cfg[i]) eq Place(sn[i] + E(data[i]`F)!Eltseq(Tq)) : i in [1..#data]] then
                    dup := true; break;
                end if;
            end for;
            if dup then continue; end if;
            rejected := false;
            for ctx in ctxs do
                S := counts_mod_p(ctx, cfg);
                if #S eq 0 then continue; end if;
                if target[Index(primes, ctx`p)] notin S then rejected := true; break; end if;
            end for;
            if rejected then continue; end if;
            covers := DoubleCoversWithBranchDivisor(E, B);
            hits := MatchCoverToCounts(E, covers, target, primes : ExtraPrimes := PrimeDivisors(D*N));
            if #hits eq 0 then continue; end if;
            Append(~seen, cfg);
            // pin the twist with more primes
            primes2 := [p : p in PrimesInInterval(3, TwistPrimes) | (D*N) mod p ne 0];
            target2 := [ComputePointsViaTrace(X, p, 1) : p in primes2];
            hits := MatchCoverToCounts(E, covers, target2, primes2 : ExtraPrimes := PrimeDivisors(D*N));
            for hh in hits do
                c := TwistedCover(covers[hh[1]], hh[2]);
                H := 0;
                if c`noether eq 5 then
                    try H := HyperellipticModelOfCover(c); catch err H := 0; end try;
                end if;
                Append(~out, <ns, cfg, hh[1], hh[2], c`noether, H>);
                vprintf HeegnerHeights: "survivor %o: cover %o twisted by %o, noether %o\n", ns, hh[1], hh[2], c`noether;
            end for;
        end for;
    end for;
    return out;
end intrinsic;

intrinsic SmallHyperellipticModel(H::CrvHyp : Shifts := 6) -> CrvHyp
{A model of the genus-3 hyperelliptic curve H over Q with smaller coefficients, obtained by rational
scalings x -> p^e x at the primes of the discriminant, shifts x -> x + r with small rational r,
x -> 1/x, and removal of square content. Every step is a Q-isomorphism.}
    f, h := HyperellipticPolynomials(H);
    require h eq 0: "a model y^2 = f(x) is expected";
    P := Parent(f);
    size := func<Q | &+[#IntegerToString(Numerator(c)) + #IntegerToString(Denominator(c)) : c in Coefficients(Q)]>;
    function clean(G)
        G *:= LCM([Denominator(c) : c in Coefficients(G)])^2;
        gg := GCD([Integers() | c : c in Coefficients(G)]);
        sq := &*[Integers() | fac[1]^(fac[2] div 2) : fac in Factorisation(gg)];
        return G div sq^2;
    end function;
    best := clean(f);
    ps := PrimeDivisors(Integers()!Discriminant(best));
    changed := true;
    while changed do
        changed := false;
        for p in ps, e in [-10..10] do
            if e eq 0 then continue; end if;
            G := clean(P![Coefficient(best, k)*p^(e*k) : k in [0..8]]);
            if size(G) lt size(best) then best := G; changed := true; end if;
        end for;
        for r in [a/b : a in [-Shifts..Shifts], b in [1..Shifts] | a ne 0] do
            G := clean(Evaluate(best, P.1 + r));
            if size(G) lt size(best) then best := G; changed := true; end if;
        end for;
        G := clean(P!Reverse(Coefficients(best) cat [0 : i in [Degree(best)+1..8]]));
        if size(G) lt size(best) then best := G; changed := true; end if;
    end while;
    return HyperellipticCurve(best);
end intrinsic;

intrinsic KleinFourEvenModel(H::CrvHyp) -> BoolElt, CrvHyp
{For a genus-3 hyperelliptic curve with an involution x -> (ax+b)/(cx+d) of the projective line having
rational fixed points, the model y^2 = G(x^2) in which that involution is x -> -x, reduced by the
scalings x -> p^e x and removal of square content, and true; false if no involution has rational fixed points. The change of
coordinate is a Q-isomorphism: only scalings by squares of x^2 are used.}
    f, h := HyperellipticPolynomials(H);
    require h eq 0: "a model y^2 = f(x) is expected";
    P := Parent(f); X := P.1;
    A, mA := AutomorphismGroup(H);
    for a in A do
        if Order(a) ne 2 then continue; end if;
        eqs := DefiningPolynomials(mA(a)); R := Parent(eqs[1]);
        a1 := MonomialCoefficient(eqs[1], R.1); b1 := MonomialCoefficient(eqs[1], R.3);
        c1 := MonomialCoefficient(eqs[3], R.1); d1 := MonomialCoefficient(eqs[3], R.3);
        if [a1, b1, c1, d1] eq [1, 0, 0, 1] then continue; end if;
        if c1 eq 0 then
            if a1/d1 ne -1 then continue; end if;
            Fn := Evaluate(f, X + b1/(2*d1));
        else
            ok, s := IsSquare((d1 - a1)^2 + 4*c1*b1);
            if not ok then continue; end if;
            r1 := (a1 - d1 + s)/(2*c1); r2 := (a1 - d1 - s)/(2*c1);
            Fn := &+[Coefficient(f, k) * (r2*X - r1)^k * (X - 1)^(8-k) : k in [0..8]];
        end if;
        if not &and[Coefficient(Fn, k) eq 0 : k in [1, 3, 5, 7]] then continue; end if;
        // reduce by x -> p^e x (even-preserving) and square content
        function clean(G)
            G *:= LCM([Denominator(c) : c in Coefficients(G)])^2;
            gg := GCD([Integers() | c : c in Coefficients(G)]);
            sq := &*[Integers() | fac[1]^(fac[2] div 2) : fac in Factorisation(gg)];
            return G div sq^2;
        end function;
        size := func<Q | &+[#IntegerToString(Integers()!c) : c in Coefficients(Q)]>;
        best := clean(Fn); changed := true;
        ps := PrimeDivisors(Integers()!Discriminant(best));
        while changed do
            changed := false;
            for p in ps, e in [-10..10] do
                if e eq 0 then continue; end if;
                G := clean(P![Coefficient(best, k)*p^(e*k) : k in [0..8]]);
                if size(G) lt size(best) then best := G; changed := true; end if;
            end for;
        end while;
        return true, HyperellipticCurve(best);
    end for;
    return false, _;
end intrinsic;
