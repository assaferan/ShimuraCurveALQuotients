// Component groups of the Jacobian of X_0^D(N)/W at a prime p | D, the reduction of CM points to them, and what
// they say about rational torsion.  paper/component-groups.tex has the statements and proofs.
//
// Y = X_0^D(N)/W (N squarefree, gcd(D, N) = 1, p | D) is a quotient of the Cerednik-Drinfeld model, so its stable
// reduction at p is a union of rational curves meeting at ordinary double points: the dual graph is the quotient by W
// of the graph of DualGraph.m, with a singular point of thickness t for an edge orbit (thickness |O^*/{+-1}| times
// the order of the stabiliser of the edge in W), and an edge fixed by an involution of W divisible by p, whose two
// ends it swaps, becoming a smooth point of the merged component.  The component group Phi_p of the Neron model is
// the critical group of the graph in which a singular point of thickness t is a chain of t - 1 components (the
// regular model), i.e. the torsion of the cokernel of its Laplacian (Raynaud; Bosch-Lutkebohmert-Raynaud 9.6).
//
// A CM point of discriminant d (p not split in its field) reduces to the singular point or component whose Eichler
// order admits an optimal embedding of the order of discriminant d (Molina, Ribet bimodules and the specialization
// of Heegner points, Theorems 5.3, 5.4, 5.8): to a singular point when p ramifies in Q(sqrt d), to a component when
// p is inert; and the number of CM points of X above a singular point or component is the number of optimal
// embeddings into its order up to conjugation by the units.  So the class in Phi_p of a degree-zero combination of
// CM cycles of the star is computable, and a nonzero class shows the divisor is not principal.
//
// Rational torsion.  Prime-to-p torsion of J(Q) injects into the special fibre of the Neron model, which is an
// extension of Phi_p by the F_p-points of the identity component; when the Hall divisor p lies in W, Frobenius
// (which acts on the Cerednik-Drinfeld graph as w_p) acts trivially on the quotient graph, so the identity component
// is a split torus with F_p-points (F_p^*)^g, of order prime to every l not dividing p - 1.  For l = p the kernel
// of reduction is the formal group of the torus, torsion-free over Z_p when p >= 3.  Hence the l-part of J(Q)_tors
// injects into Phi_p for every l not dividing p - 1, p = l included when p >= 3, and in particular ALL odd torsion
// injects into Phi_2.  A torsion CM divisor whose class in Phi_p vanishes therefore has trivial l-part for those l.
//
// Orders of ideal classes are not reproducible between Magma runs, so a record carries the ideals it was built from.

declare verbose ComponentGroups, 1;

CDQuotFmt := recformat< D : RngIntElt, N : RngIntElt, W : SetEnum, p : RngIntElt,
    L : SeqEnum, LL : SeqEnum, h : RngIntElt, hp : RngIntElt, org : SeqEnum, alE : SeqEnum, alV : SeqEnum,
    uV : SeqEnum, uE : SeqEnum,
    vorb : SeqEnum, nv : RngIntElt, eorb : SeqEnum, edges : SeqEnum, thick : SeqEnum,
    V : RngIntElt, chain : SeqEnum, Sm : Mtrx, Pm : Mtrx, inv : SeqEnum >;

function cg_find(J, reps)
    hits := [i : i in [1..#reps] | IsIsomorphic(J, reps[i])];
    assert #hits eq 1;
    return hits[1];
end function;

function cg_img(list, q)
    return [t[2] : t in list | t[1] eq q][1];
end function;

// optimal embeddings of the order of discriminant d into the definite order Ov, up to conjugation by Ov^*
function cg_optemb(d, Ov)
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

// the order of the stabiliser in W of a CM point of discriminant d of X, times |Aut| of the point: the e_d of the
// dictionary of Guo-Yang (3 at -3, 2 at -4, doubled at the fixed points of an Atkin-Lehner involution)
function cg_stab(D, N, W, d)
    e := d eq -3 select 3 else (d eq -4 select 2 else 1);
    if exists{m : m in W | m gt 1 and (d eq -4*m or (m mod 4 eq 3 and d eq -m))} then e *:= 2; end if;
    return e;
end function;

intrinsic ComponentGroupData(D::RngIntElt, N::RngIntElt, W::SetEnum, p::RngIntElt) -> Rec
{The special fibre at p | D of X_0^D(N)/W as a record: the Cerednik-Drinfeld graph built from one list of ideal
classes (vertices L, edges LL, with the Atkin-Lehner actions and unit orders, validated by ValidateDualGraphData),
its quotient by W with the thicknesses of the singular points, and the Smith form of the Laplacian of the regular
model's graph, whose torsion is the component group.  W is a group of Hall divisors of DN containing 1.}
    require DualGraphApplicable(D, N) : "needs D > 1, N squarefree, gcd(D, N) = 1";
    require IsPrime(p) and D mod p eq 0 : "p must be a prime dividing D";
    require 1 in W and forall{m : m in W | (D*N) mod m eq 0 and GCD(m, (D*N) div m) eq 1} : "W must consist of Hall divisors of DN and contain 1";
    Q := QuaternionAlgebra(D div p); O := QuaternionOrder(Q, N); OO := Order(O, p);
    L := LeftIdealClasses(O); LL := LeftIdealClasses(OO);
    h := #L; hp := #LL;
    org := [cg_find(lideal< O | Basis(II) >, L) : II in LL];
    alE := []; alV := [];
    for l in PrimeDivisors(D*N) do
        P := PrimeIdeal(OO, l); assert Norm(P) eq l and IsTwoSidedIdeal(P);
        Append(~alE, <l, [cg_find(P*II, LL) : II in LL]>);
        if l ne p then
            P0 := PrimeIdeal(O, l); assert Norm(P0) eq l and IsTwoSidedIdeal(P0);
            Append(~alV, <l, [cg_find(P0*I, L) : I in L]>);
        end if;
    end for;
    uV := [#UnitGroup(RightOrder(I)) : I in L]; uE := [#UnitGroup(RightOrder(II)) : II in LL];
    ValidateDualGraphData(D, N, p, <h, hp, org, alE, alV, uV, uE>);
    vprintf ComponentGroups : "p = %o: h = %o, h' = %o, graph data validated\n", p, h, hp;
    // vertex x in [1..2h]: class ((x-1) mod h) + 1 in the copy + if x <= h, - otherwise; w_q acts on classes for
    // q | DN/p and swaps the copies for q = p
    function vact(m, x)
        lab := ((x - 1) mod h) + 1; sgn := x le h;
        for q in PrimeDivisors(m) do if q eq p then sgn := not sgn; else lab := cg_img(alV, q)[lab]; end if; end for;
        return sgn select lab else lab + h;
    end function;
    function eact(m, e) for q in PrimeDivisors(m) do e := cg_img(alE, q)[e]; end for; return e; end function;
    wp := cg_img(alE, p);
    vorb := [0 : x in [1..2*h]]; nv := 0;
    for x in [1..2*h] do if vorb[x] eq 0 then nv +:= 1; for m in W do vorb[vact(m, x)] := nv; end for; end if; end for;
    // eorb[e]: the edge orbit of e in the quotient, or -1 when an involution of W divisible by p fixes e and so
    // swaps its ends, making it a smooth point of the merged component
    eorb := [0 : e in [1..hp]]; edges := []; thick := []; nedge := 0;
    for e in [1..hp] do
        if eorb[e] ne 0 then continue; end if;
        S := { eact(m, e) : m in W };
        if exists{m : m in W | m mod p eq 0 and eact(m, e) eq e} then
            assert vorb[org[e]] eq vorb[org[wp[e]] + h];
            for x in S do eorb[x] := -1; end for; continue;
        end if;
        nedge +:= 1; for x in S do eorb[x] := nedge; end for;
        Append(~edges, <vorb[org[e]], vorb[org[wp[e]] + h]>);
        Append(~thick, (uE[e] div 2) * (#W div #S));
    end for;
    // the regular model: a singular point of thickness t is a chain of t - 1 components
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
    vprintf ComponentGroups : "quotient graph: %o vertices, %o edges of thicknesses %o, b1 = %o; Phi_%o = %o\n", nv, #edges, thick, #edges - nv + 1, p, inv;
    return rec< CDQuotFmt | D := D, N := N, W := W, p := p, L := L, LL := LL, h := h, hp := hp, org := org, alE := alE,
        alV := alV, uV := uV, uE := uE, vorb := vorb, nv := nv, eorb := eorb, edges := edges, thick := thick,
        V := V, chain := chain, Sm := Sm, Pm := Pm, inv := inv >;
end intrinsic;

intrinsic ComponentGroup(G::Rec) -> SeqEnum
{The invariants of the component group Phi_p of the Jacobian of X_0^D(N)/W at p, from ComponentGroupData.}
    return G`inv;
end intrinsic;

intrinsic QuotientGraph(G::Rec) -> RngIntElt, SeqEnum, SeqEnum
{The dual graph of the stable reduction of X_0^D(N)/W at p, from ComponentGroupData: the number of vertices, the
edges as pairs of vertices, and their thicknesses.}
    return G`nv, G`edges, G`thick;
end intrinsic;

intrinsic IsTotallyDegenerate(G::Rec, g::RngIntElt) -> BoolElt
{Whether the reduction described by ComponentGroupData is totally degenerate for a curve of genus g: the quotient
graph is connected with first Betti number g, so the identity component of the Neron model is a torus.}
    return #[i : i in [1..G`V] | G`Sm[i,i] eq 0] eq 1 and #G`edges - G`nv + 1 eq g;
end intrinsic;

intrinsic CMSpecialisation(G::Rec, d::RngIntElt) -> Assoc, RngIntElt
{Where the CM points of discriminant d of the star X_0^D(N)/W reduce at p, from ComponentGroupData: an associative
array vertex -> number of star points on that component of the regular model's graph (a singular point of
thickness t is its chain of t - 1 components, numbered after the vertices of the quotient graph), and the number
of star points of discriminant d.  Requires p not split in Q(sqrt d).  A CM point reduces to a singular point or
component whose Eichler order admits an optimal embedding of the order of discriminant d (Molina, Theorems 5.3,
5.4, 5.8); a singular point whose ends are swapped by W is a smooth point of the merged component; a singular
point that survives in the quotient carries the point on the middle component of its chain.}
    D := G`D; N := G`N; W := G`W; p := G`p;
    d0 := FundamentalDiscriminant(d); sym := KroneckerSymbol(d0, p);
    require sym ne 1 : "p splits in the CM field: ordinary reduction, not handled";
    cnt := AssociativeArray(); total := 0;
    if sym eq -1 then
        for i in [1..G`h] do
            c := cg_optemb(d, RightOrder(G`L[i]));
            if c eq 0 then continue; end if;
            for x in [i, i + G`h] do v := G`vorb[x]; cnt[v] := (IsDefined(cnt, v) select cnt[v] else 0) + c; total +:= c; end for;
        end for;
    else
        for e in [1..G`hp] do
            c := cg_optemb(d, RightOrder(G`LL[e]));
            if c eq 0 then continue; end if;
            o := G`eorb[e];
            if o eq -1 then v := G`vorb[G`org[e]];
            else t := G`thick[o]; require IsEven(t) : "a CM point at a singular point of odd thickness"; v := t eq 2 select G`chain[o][1] else G`chain[o][t div 2]; end if;
            cnt[v] := (IsDefined(cnt, v) select cnt[v] else 0) + c; total +:= c;
        end for;
    end if;
    ratio := #W div cg_stab(D, N, W, d);                 // CM points of X above one star point
    require total mod ratio eq 0 and forall{v : v in Keys(cnt) | cnt[v] mod ratio eq 0} : "embedding counts are not multiples of the orbit size";
    for v in Keys(cnt) do cnt[v] := cnt[v] div ratio; end for;
    return cnt, total div ratio;
end intrinsic;

intrinsic ClassInComponentGroup(G::Rec, divisor::SeqEnum) -> SeqEnum
{The class in Phi_p of the divisor sum m_d Z_d of the star, Z_d the full CM cycle of discriminant d (the sum of the
star points of discriminant d), given as pairs <d, m_d> of total degree zero; the answer is the vector of
coordinates with respect to the invariants of ComponentGroup(G), zero meaning the divisor reduces to the identity
component.}
    v := ZeroMatrix(Integers(), G`V, 1);
    for t in divisor do
        cnt := CMSpecialisation(G, t[1]);
        for w in Keys(cnt) do v[w, 1] +:= t[2] * cnt[w]; end for;
    end for;
    require &+Eltseq(v) eq 0 : "the divisor must have degree zero";
    w := G`Pm * v;
    return [ w[i,1] mod G`Sm[i,i] : i in [1..G`V] | G`Sm[i,i] notin {0, 1} ];
end intrinsic;

intrinsic IsTorsionPartInjective(G::Rec, l::RngIntElt) -> BoolElt
{Whether the l-part of the rational torsion of the Jacobian of X_0^D(N)/W injects into the component group at p
described by G, which requires the Hall divisor p to lie in W (so that Frobenius acts trivially on the quotient
graph and the identity component is a split torus): true exactly when l does not divide p - 1, and, for l = p,
when p is at least 3.}
    require IsPrime(l) : "l must be prime";
    require G`p in G`W : "the torus is known to be split only when the Hall divisor p lies in W";
    if l eq G`p then return G`p ge 3; end if;
    return (G`p - 1) mod l ne 0;
end intrinsic;

intrinsic TorsionOrderBoundOfCMDivisor(D::RngIntElt, N::RngIntElt, W::SetEnum, divisor::SeqEnum : Bound := 0, Primes := []) -> RngIntElt, BoolElt, SeqEnum
{For a degree-zero combination of CM cycles of the star X_0^D(N)/W (pairs <d, m_d>) that is known to be torsion
in the Jacobian, a bound on its order: the bound on the rational torsion (JacobianTorsionBound unless Bound is
given) with the l-parts removed for every l whose part injects into some component group Phi_p, p in Primes
(default: the primes of D at which p is unsplit in every Q(sqrt d) of the divisor), in which the divisor has class
zero.  Also returns false if some class is nonzero, which proves the divisor is not principal, and the list of
<p, invariants of Phi_p, class>.}
    if Bound eq 0 then
        X := CreateShimuraQuot(D, N, W); X`g := GenusShimuraCurveQuotient(D, N, W); X`CurveID := 0;
        B := JacobianTorsionBound(X);
    else
        B := Bound;
    end if;
    ps := #Primes eq 0 select [p : p in PrimeDivisors(D) | forall{t : t in divisor | KroneckerSymbol(FundamentalDiscriminant(t[1]), p) ne 1}] else Primes;
    principal := true; record := [];
    for p in ps do
        G := ComponentGroupData(D, N, W, p);
        c := ClassInComponentGroup(G, divisor);
        Append(~record, <p, G`inv, c>);
        if exists{x : x in c | x ne 0} then principal := false; continue; end if;
        for l in PrimeDivisors(B) do
            if p in W and IsTorsionPartInjective(G, l) then
                while B mod l eq 0 do B := B div l; end while;
            end if;
        end for;
    end for;
    return B, principal, record;
end intrinsic;
