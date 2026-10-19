// Double covers of an elliptic curve with a prescribed branch divisor, decided exactly.
//
// A double cover X -> E of an elliptic curve E/Q branched along a reduced Q-rational divisor B of
// degree 4 has genus 3 and function field Q(E)(sqrt f) with div f = B - 2D for a Q-rational
// divisor D of degree 2.  Writing D ~ (Q) + (O), the condition 2D ~ B says 2Q = P_B for the point
// P_B with B ~ (P_B) + 3(O), so the covers over Q are indexed by the rational halves of P_B: at
// most four, a torsor under E(Q)[2], each determined up to quadratic twist by the half Q.
//
// Everything here is exact: divisors and Riemann-Roch spaces on E, the genus and the holomorphic
// differentials of the cover from its function field, hyperellipticity by the number of quadrics
// through the canonical image (one quadric in genus 3 if and only if hyperelliptic), and point
// counts of each candidate from its reduction, which single out the quotient among the candidates.

declare verbose ExactDoubleCover, 2;

// Q-linear relations among elements of a function field F / A / Q(x), F and A of finite degree.
// Coordinates are flattened down the tower and cleared of denominators, so a relation among the
// coefficient vectors is a relation among the elements.
function q_relations(elts)
    coords := [[c : c in Eltseq(a), a in Eltseq(e)] : e in elts];
    den := LCM([Denominator(c) : c in cs, cs in coords]);
    polys := [[Numerator(c*den) : c in cs] : cs in coords];
    nslots := #polys[1];
    deg := Max([Degree(p) : p in ps, ps in polys] cat [0]);
    rows := [&cat[[Coefficient(p, i) : i in [0..deg]] : p in ps] : ps in polys];
    M := Matrix(Rationals(), #elts, nslots*(deg+1), rows);
    return Basis(Kernel(M));
end function;

// The rational point P with D ~ (P) - deg(D)*(O) for a Q-rational divisor D on E.
function point_of_class(E, D)
    O := Divisor(Place(E!0));
    D1 := D - (Degree(D) - 1)*O;
    V, mp := RiemannRochSpace(D1);
    assert Dimension(V) eq 1;
    Deff := D1 + Divisor(mp(V.1));
    supp, mults := Support(Deff);
    assert #supp eq 1 and mults[1] eq 1 and Degree(supp[1]) eq 1;
    return E!RepresentativePoint(supp[1]);
end function;

// The reduction of a function f on E modulo p as an element of the function field of E mod p,
// or false when the coefficients are not p-integral after clearing denominators.
function reduce_function(E, f, p, FEp)
    pf := ProjectiveFunction(f);
    Rq := Parent(Numerator(pf));
    R3Q := Type(Rq) eq RngMPolRes select PreimageRing(Rq) else Rq;
    num := R3Q!Numerator(pf); den := R3Q!Denominator(pf);
    cs := Coefficients(num) cat Coefficients(den);
    scale := LCM([Denominator(c) : c in cs]) / GCD([Numerator(c) : c in cs]);
    num *:= scale; den *:= scale;
    R3 := PolynomialRing(GF(p), 3);
    nump := R3!num; denp := R3!den;
    if denp eq 0 or nump eq 0 then return false, _; end if;
    return true, Evaluate(nump, [FEp.1, FEp.2, 1]) / Evaluate(denp, [FEp.1, FEp.2, 1]);
end function;

CoverRF := recformat< Q : PtEll, f : FldFunFracSchElt, F : FldFun, genus : RngIntElt,
                      noether : RngIntElt, diffs : SeqEnum, afmap : Map, B : DivCrvElt >;

// The cover y^2 = f of E as a record: its function field, genus, and the span of the pairwise
// products of a basis of holomorphic differentials.
function cover_record(m, f, Q, B)
    AF := Codomain(m);
    F := FunctionField(Polynomial([AF | -m(f), 0, 1]));
    g := Genus(F);
    diffs := BasisOfDifferentialsFirstKind(F);
    gs := [d / diffs[1] : d in diffs];
    prods := [gs[i]*gs[j] : j in [i..#gs], i in [1..#gs]];
    nrel := #q_relations(prods);
    return rec< CoverRF | Q := Q, f := f, F := F, genus := g, noether := #prods - nrel,
                          diffs := gs, afmap := m, B := B >;
end function;

intrinsic DoubleCoversWithBranchDivisor(E::CrvEll, B::DivCrvElt) -> List
{Every double cover of E over Q branched exactly along the reduced Q-rational divisor B of degree
4, up to quadratic twist. Each entry is a record with the half Q of the class point, the function
f with div f = B - 2((Q) + (O)), the function field of y^2 = f, its genus, and the dimension of
the span of the pairwise products of a basis of holomorphic differentials (5 if hyperelliptic,
6 if not).}
    _, mults := Support(B);
    require Degree(B) eq 4 and IsEffective(B) and &and[m eq 1 : m in mults]:
        "the branch divisor must be reduced, effective, of degree 4";
    FE := FunctionField(E);
    O := Divisor(Place(E!0));
    PB := point_of_class(E, B);
    halves := DivisionPoints(PB, 2);
    vprintf ExactDoubleCover: "class point %o has %o rational halves\n", PB, #halves;
    AF, m := AlgorithmicFunctionField(FE);
    out := [* *];
    for Q in halves do
        D := Divisor(Place(Q)) + O;
        ok, f := IsPrincipal(B - 2*D);
        assert ok;
        C := cover_record(m, f, Q, B);
        Append(~out, C);
        vprintf ExactDoubleCover: "  half %o: genus %o, products span %o\n", Q, C`genus, C`noether;
    end for;
    return out;
end intrinsic;

intrinsic DoubleCoversWithBranchDivisor(E::CrvEll, Bpts::List) -> List
{The same for the branch divisor spanned by the points in Bpts, points of E over number fields
of which conjugates may be listed or omitted.}
    places := {@ Place(P) : P in Bpts @};
    return DoubleCoversWithBranchDivisor(E, &+[Divisor(pl) : pl in places]);
end intrinsic;

intrinsic TwistedCover(C::Rec, kappa::RngIntElt) -> Rec
{The cover y^2 = kappa*f: the quadratic twist of the cover C along its deck involution, which is
not the twist of a hyperelliptic model of C along its hyperelliptic involution.}
    return cover_record(C`afmap, kappa*C`f, C`Q, C`B);
end intrinsic;

intrinsic HyperellipticModelOfCover(C::Rec) -> CrvHyp, FldFunElt, FldFunElt
{For a genus-3 cover from DoubleCoversWithBranchDivisor whose canonical image is a conic with a
rational point, a hyperelliptic model y^2 = G(t) over Q together with the functions t, y realising
it in the function field of the cover.}
    require C`genus eq 3 and C`noether eq 5: "the cover is not a hyperelliptic curve of genus 3";
    F := C`F; gs := C`diffs;
    prods := [gs[i]*gs[j] : j in [i..3], i in [1..3]];
    idx := [[i,j] : j in [i..3], i in [1..3]];
    rel := q_relations(prods)[1];
    P2<X1,X2,X3> := ProjectiveSpace(Rationals(), 2);
    quadric := &+[rel[k]*P2.idx[k][1]*P2.idx[k][2] : k in [1..#idx]];
    conic := Conic(P2, quadric);
    ok, pt := HasRationalPoint(conic);
    require ok: "the canonical conic has no rational point";
    // Projection from the rational point: the pencil of lines through it has degree one on the
    // conic, hence degree two on the curve.
    v := Vector(Rationals(), Eltseq(pt));
    lines := Basis(Kernel(Transpose(Matrix(v))));
    l1 := &+[lines[1][i]*gs[i] : i in [1..3]];
    l2 := &+[lines[2][i]*gs[i] : i in [1..3]];
    t := l1 / l2;
    // A second generator z and its minimal polynomial over Q(t), which must be quadratic: a
    // generator lying in Q(t) is skipped, and a relation of lower degree in z is looked for first
    // so that the quadratic one is the minimal polynomial and not a multiple of a linear relation.
    Pt<T> := PolynomialRing(Rationals());
    for z in [F.1, F!BaseField(F).1, F.1*BaseField(F).1] do
        for n in [2..40] do
            if #q_relations([t^i*z^j : i in [0..n], j in [0..1]]) gt 0 then break; end if;
            rels := q_relations([t^i*z^j : i in [0..n], j in [0..2]]);
            if #rels eq 0 then continue; end if;
            r := rels[1];
            cf := [Pt![r[j*(n+1) + i + 1] : i in [0..n]] : j in [0..2]];
            A := cf[3]; Bc := cf[2]; Cc := cf[1];
            G := Bc^2 - 4*A*Cc;
            vprintf ExactDoubleCover, 2: "relation at n = %o: A = %o, B = %o, C = %o\n", n, A, Bc, Cc;
            // y = 2Az + B has y^2 = G(t); strip square factors from G.
            sq := Pt!1;
            for fac in Factorisation(G) do
                sq *:= fac[1]^(fac[2] div 2);
            end for;
            G0 := G div sq^2;
            // An integral model with no square factor in the content.
            lcm := LCM([Denominator(c) : c in Coefficients(G0)]);
            G0 *:= lcm^2;
            cont := GCD([Integers() | c : c in Coefficients(G0)]);
            csq := &*[Integers() | fac[1]^(fac[2] div 2) : fac in Factorisation(cont)];
            G0 := G0 div csq^2;
            H := HyperellipticCurve(G0);
            y := (2*Evaluate(A, t)*z + Evaluate(Bc, t)) * (lcm/csq) / Evaluate(sq, t);
            return H, t, y;
        end for;
    end for;
    require false: "no quadratic relation found";
end intrinsic;

intrinsic HyperellipticModelOverConic(C::Rec) -> RngMPolElt, RngMPolElt, SeqEnum, FldFunElt
{For a genus-3 hyperelliptic cover from DoubleCoversWithBranchDivisor, a model
    a X1^2 + b X2^2 + c X3^2 = 0,   y^2 = G(X1, X2, X3)
over Q with a, b, c integers, G a quartic form with integral coefficients, no square content and
no monomial divisible by X3^2 (it is the representative of the class of y^2 modulo the conic in
that complement); X1, X2, X3 have weight 1 and y weight 2. The model exists whether or not the
canonical conic has a rational point. Also returned are the functions X1, X2, X3, y of the
function field of the cover realising it.}
    require C`genus eq 3 and C`noether eq 5: "the cover is not a hyperelliptic curve of genus 3";
    F := C`F; gs := C`diffs;
    prods := [gs[i]*gs[j] : j in [i..3], i in [1..3]];
    idx := [[i,j] : j in [i..3], i in [1..3]];
    rel := q_relations(prods)[1];
    R<X1,X2,X3> := PolynomialRing(Rationals(), 3);
    quadric := &+[rel[k]*R.idx[k][1]*R.idx[k][2] : k in [1..#idx]];
    // H^0(2K) = span of the products g_i g_j plus one more function b; the function y = b + sum c_k g_i g_j
    // anti-invariant under the hyperelliptic involution is the one whose square lies in the span of the
    // quartic monomials in the g_i, which is a linear condition on the c_k
    omega := BasisOfDifferentialsFirstKind(F)[1];
    L2, m2 := RiemannRochSpace(2*Divisor(omega));
    assert Dimension(L2) eq 6;
    b := 0;
    for e in [m2(L2.i) : i in [1..6]] do
        if #q_relations(prods cat [e]) eq #q_relations(prods) then b := e; break; end if;
    end for;
    assert b ne 0;
    mons4 := MonomialsOfDegree(R, 4);
    vals := [Evaluate(m, gs) : m in mons4];
    rels := q_relations([b^2] cat [b*q : q in prods] cat vals);
    r := [v : v in rels | v[1] ne 0][1];
    y := b + &+[r[1 + k] / (2*r[1]) * prods[k] : k in [1..#prods]];
    rels := q_relations(vals cat [y^2]);
    r := [v : v in rels | v[#mons4 + 1] ne 0][1];
    G := -&+[r[k]*mons4[k] : k in [1..#mons4]] / r[#mons4 + 1];
    // diagonalise the conic, x = (s_j w_j)_j T with T Gram T^t diagonal, and make the diagonal
    // entries squarefree integers with no common factor
    Gram := Matrix(Rationals(), 3, 3, [MonomialCoefficient(quadric, R.i*R.j) * (i eq j select 1 else 1/2) : i, j in [1..3]]);
    Dg, T := OrthogonalizeGram(Gram);
    scal := [];
    for i in [1..3] do
        dd := Dg[i,i]; nn := Numerator(dd)*Denominator(dd);
        sq := &*[Integers() | pf[1]^(pf[2] div 2) : pf in Factorization(nn)];
        Append(~scal, Denominator(dd)/sq);
    end for;
    diag := [Dg[i,i]*scal[i]^2 : i in [1..3]];
    g := GCD([Integers() | d : d in diag]);
    diag := [Integers() | d/g : d in diag];
    subst := [&+[T[j,i]*scal[j]*R.j : j in [1..3]] : i in [1..3]];
    conic := &+[diag[i]*R.i^2 : i in [1..3]];
    assert Evaluate(quadric, subst) eq g*conic;
    Ti := T^-1;
    ws := [(&+[gs[i]*Ti[i,j] : i in [1..3]]) / scal[j] : j in [1..3]];
    Gw := Evaluate(G, subst);
    // the representative without X3^2: on the conic X3^2 = -(a X1^2 + b X2^2)/c
    while Degree(Gw, 3) ge 2 do
        d3 := Degree(Gw, 3); co := Coefficient(Gw, 3, d3);
        Gw := Gw - co*X3^d3 - co*X3^(d3-2)*(diag[1]*X1^2 + diag[2]*X2^2)/diag[3];
    end while;
    den := LCM([Denominator(co) : co in Coefficients(Gw)]);
    Gw *:= den^2;
    cont := GCD([Integers() | co : co in Coefficients(Gw)]);
    csq := &*[Integers() | pf[1]^(pf[2] div 2) : pf in Factorization(cont)];
    Gw := Gw / csq^2;
    y *:= den/csq;
    assert Evaluate(conic, ws) eq 0 and y^2 eq Evaluate(Gw, ws);
    return conic, Gw, ws, y;
end intrinsic;

intrinsic PointCountOverConic(conic::RngMPolElt, G::RngMPolElt, q::RngIntElt) -> RngIntElt
{The number of points over F_q of the curve  conic = 0, y^2 = G  in the weighted projective space
with X1, X2, X3 of weight 1 and y of weight 2, for q a prime power coprime to the coefficients of
the conic: each of the q + 1 points of the conic carries 1 + chi(G) points, chi the quadratic
character of F_q.}
    Fq := GF(q);
    P2 := ProjectiveSpace(Fq, 2);
    Cq := Conic(P2, Evaluate(conic, [P2.1, P2.2, P2.3]));
    Gq := ChangeRing(G, Fq);
    n := 0;
    for pt in Points(Cq) do
        v := Evaluate(Gq, Eltseq(pt));
        n +:= v eq 0 select 1 else (IsSquare(v) select 2 else 0);
    end for;
    return n;
end intrinsic;

intrinsic CoverPointCounts(E::CrvEll, C::Rec, primes::SeqEnum) -> SeqEnum
{For each prime p in primes, the pair <#X(F_p), #X'(F_p)> for the cover y^2 = f and its
quadratic twist by a non-square, or <0,0> when the reduction is not available (p | 2, bad
reduction or a non-integral coefficient of the model of E, non-integral coefficients of f, two
branch points meeting modulo p, or the reduced cover has a different genus).}
    out := [];
    for p in primes do
        if p eq 2 or Valuation(Discriminant(E), p) ne 0 or exists{a : a in aInvariants(E) | Valuation(a, p) lt 0} then
            Append(~out, <0,0>); continue;
        end if;
        Ep := ChangeRing(E, GF(p));
        FEp := FunctionField(Ep);
        ok, fp := reduce_function(E, C`f, p, FEp);
        if not ok then Append(~out, <0,0>); continue; end if;
        AFp, mp := AlgorithmicFunctionField(FEp);
        // the reduced function can be a square when two branch points meet modulo p
        try
            Fp := FunctionField(Polynomial([AFp | -mp(fp), 0, 1]));
            g := Genus(Fp);
        catch err
            Append(~out, <0,0>); continue;
        end try;
        if g ne C`genus then Append(~out, <0,0>); continue; end if;
        n := NumberOfPlacesOfDegreeOneECF(RationalExtensionRepresentation(Fp));
        Append(~out, <n, 2*#Ep(GF(p)) - n>);
    end for;
    return out;
end intrinsic;

// Primes dividing a numerator or denominator of a coefficient of the function f on E.
function coefficient_primes(f)
    pf := ProjectiveFunction(f);
    cs := Coefficients(Numerator(pf)) cat Coefficients(Denominator(pf));
    return &join[Set(PrimeDivisors(Numerator(c)*Denominator(c))) : c in cs];
end function;

// Primes modulo which two of the geometric points of the divisor B can coincide: the prime
// divisors of the norms of the ideals of 2x2 minors of pairs of integral projective coordinate
// vectors, over a splitting field of the residue fields of B.
function branch_collision_primes(E, B)
    supp := Support(B);
    flds := [* Ring(Parent(RepresentativePoint(pl))) : pl in supp *];
    polys := [DefiningPolynomial(AbsoluteField(K)) : K in flds | Type(K) ne FldRat];
    L := #polys eq 0 select RationalsAsNumberField() else SplittingField(polys);
    OL := MaximalOrder(L);
    pts := [];
    for i in [1..#supp] do
        P := Eltseq(RepresentativePoint(supp[i]));
        K := flds[i];
        if Type(K) eq FldRat then
            embs := [map< K -> L | a :-> L!a >];
        else
            KA := AbsoluteField(K);
            embs := [hom< KA -> L | r[1] > : r in Roots(DefiningPolynomial(KA), L)];
            P := [KA!c : c in P];
        end if;
        for e in embs do
            v := [e(c) : c in P];
            den := LCM([Integers() | Denominator(c*OL) : c in v | c ne 0]);
            Append(~pts, [OL | c*den : c in v]);
        end for;
    end for;
    assert #pts eq Degree(B);
    primes := {};
    for i in [1..#pts], j in [i+1..#pts] do
        a := pts[i]; b := pts[j];
        minors := [a[1]*b[2] - a[2]*b[1], a[1]*b[3] - a[3]*b[1], a[2]*b[3] - a[3]*b[2]];
        I := ideal< OL | minors >;
        if I ne 1*OL then primes join:= Set(PrimeDivisors(Norm(I))); end if;
    end for;
    return primes;
end function;

intrinsic MatchCoverToCounts(E::CrvEll, covers::List, target::SeqEnum, primes::SeqEnum :
                             ExtraPrimes := []) -> SeqEnum
{Which of the covers, twisted by which squarefree kappa, has the point counts target[i] at
primes[i] (a target of 0 means no count is available at that prime). Returns the <cover index,
kappa> pairs consistent with every usable prime. A cover is excluded as soon as one prime matches
neither it nor its twist, whatever kappa. The sign of the twist at each remaining prime then
determines the Kronecker symbols of kappa, and kappa is solved for as a product of -1, 2, the
primes in ExtraPrimes (which must contain the primes of bad reduction of the curve being
identified) and the primes where the untwisted cover can have bad reduction: those dividing the
discriminant of E or a coefficient of f, and those modulo which two branch points coincide. At
any other odd prime the cover y^2 = f is smooth, so a twist ramified there would have bad
reduction where the identified curve has good reduction. Two pairs with the same cover index
mean the primes used do not separate the twists; kappa = 0 means no usable prime distinguished
the cover from its twist. kappa refers to the function f of the cover as returned, which is only
determined up to a constant: a kappa computed in one session goes with that session's f.}
    disc := Discriminant(E);
    hits := [];
    for ci in [1..#covers] do
        counts := CoverPointCounts(E, covers[ci], primes);
        S := Set(PrimeDivisors(2*Numerator(disc)*Denominator(disc))) join Set(ExtraPrimes)
             join &join[Set(PrimeDivisors(Denominator(a))) : a in aInvariants(E)]
             join coefficient_primes(covers[ci]`f) join branch_collision_primes(E, covers[ci]`B);
        vprintf ExactDoubleCover: "cover %o: twist supported on %o\n", ci, Sort(Setseq(S));
        gens := [-1] cat Sort(Setseq(S));
        rows := []; rhs := [];
        good := true;
        for i in [1..#primes] do
            p := primes[i];
            if target[i] eq 0 or counts[i][1] eq 0 or p in S then continue; end if;
            if target[i] eq counts[i][1] then sgn := 1;
            elif target[i] eq counts[i][2] then sgn := -1;
            else good := false; break; end if;
            if counts[i][1] eq counts[i][2] then continue; end if;   // no information on the sign
            Append(~rows, [GF(2) | KroneckerSymbol(q, p) eq -1 select 1 else 0 : q in gens]);
            Append(~rhs, GF(2)!(sgn eq -1 select 1 else 0));
        end for;
        if not good then continue; end if;
        if #rows eq 0 then Append(~hits, <ci, 0>); continue; end if;   // consistent, twist undetermined
        M := Matrix(GF(2), rows);
        ok, sol, ker := IsConsistent(Transpose(M), Vector(GF(2), rhs));
        if not ok then
            vprintf ExactDoubleCover: "cover %o: the signs of the twist are not those of a character on the support\n", ci;
            continue;
        end if;
        for k in ker do
            v := sol + k;
            kappa := &*[Integers() | gens[j] : j in [1..#gens] | v[j] eq 1];
            Append(~hits, <ci, kappa>);
        end for;
    end for;
    return hits;
end intrinsic;
