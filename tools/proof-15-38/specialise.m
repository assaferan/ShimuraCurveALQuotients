// THE CLASS OF Delta = P_-228 + P_-760 - P_-12 - P_-40 IN THE COMPONENT GROUP Phi_p OF THE STAR OF X_0^15(38), p | D.
// Everything is built from ONE set of ideal classes (the orderings of QuaternionOrder / LeftIdealClasses are
// not reproducible), following dg_build of DualGraph.m: vertices (class of the level-N Eichler order O in
// B_{D/p}) x {+,-}, edges = classes of the level-Np order OO, edge e joining (org[e], +) to (org[w_p e], -).
// CM points (Molina, Thm 4.4 and its square-free-N form): p inert in K -> a component, the vertex class whose
// right order admits the optimal embedding; p ramified -> a singular point, the edge class likewise.  Under the
// quotient by W, an edge fixed by an m with p | m (ends swapped) becomes a smooth point of the merged component;
// an edge orbit otherwise is a singular point of thickness |O_R(I_e)^*/+-1| * |Stab_W(e)| and a point on it sits
// at its midpoint.  Phi_p = critical group of the subdivided quotient graph (Raynaud).
//   magma -b PP:=3 spec1538b.m          (run from wt-metric)
AttachSpec("ShimuraQuotients.spec");
D := 15; N := 38; p := StringToInteger(PP);
W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
Dp := D div p;
Q := QuaternionAlgebra(Dp); O := QuaternionOrder(Q, N); OO := Order(O, p);
L := LeftIdealClasses(O); LL := LeftIdealClasses(OO);
h := #L; hp := #LL;
printf "p = %o: h = %o, h' = %o\n", p, h, hp;
function Find(J, reps) hits := [i : i in [1..#reps] | IsIsomorphic(J, reps[i])]; assert #hits eq 1; return hits[1]; end function;
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
data := <h, hp, org, alE, alV, uV, uE>;
ValidateDualGraphData(D, N, p, data);
printf "graph data validated (unit orders: vertices %o, edges %o)\n", Set(uV), Set(uE);
function img(list, q) for t in list do if t[1] mod q eq 0 then return t[2]; end if; end for; error "no such involution"; end function;
function vact(m, x)
    lab := ((x - 1) mod h) + 1; sgn := x le h;
    for q in PrimeDivisors(m) do if q eq p then sgn := not sgn; else lab := img(alV, q)[lab]; end if; end for;
    return sgn select lab else lab + h;
end function;
function eact(m, e) for q in PrimeDivisors(m) do e := img(alE, q)[e]; end for; return e; end function;
wp := img(alE, p);
vorb := [0 : x in [1..2*h]]; nv := 0;
for x in [1..2*h] do if vorb[x] eq 0 then nv +:= 1; for m in W do vorb[vact(m, x)] := nv; end for; end if; end for;
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
printf "quotient graph: %o vertex orbits, edges %o, thicknesses %o (total %o)\n", nv, edges, thick, &+thick;
// cross-check with the library routine (its lengths are twice the thicknesses)
nvm, edges2, lens2, _ := DualGraphMetricQuotientFromData(data, p, W);
assert nvm eq nv and Multiset(lens2) eq Multiset([2*t : t in thick]);
function OptEmb(d, Ov)
    t0 := d mod 2; n := (t0^2 - d) div 4;
    Zb := ZBasis(Ov);
    G := Matrix(Rationals(), 4, 4, [ Trace(Zb[i] * Conjugate(Zb[j])) : i, j in [1..4] ]);
    S := ShortVectors(LatticeWithGram(G), 2*n, 2*n);
    d0 := FundamentalDiscriminant(d); _, f := IsSquare(d div d0);
    cnt := 0;
    for s in S do for sg in [1, -1] do
        x := &+[ sg * s[1][i] * Zb[i] : i in [1..4] ];
        if Trace(x) ne t0 then continue; end if;
        if exists{<q, c> : c in [0..q-1], q in PrimeDivisors(f) | (x - c)/q in Ov} then continue; end if;
        cnt +:= 1;
    end for; end for;
    return cnt;
end function;
// the subdivided graph: orbit vertices 1..nv, then thick[i]-1 interior vertices per edge
V := nv; sub := [];   // sub[i] = the list of vertices along edge i from end 1 to end 2 (interior ones)
E := [];
for i->e in edges do
    chain := []; prev := e[1];
    for k in [1..thick[i]-1] do V +:= 1; Append(~chain, V); Append(~E, <prev, V>); prev := V; end for;
    Append(~E, <prev, e[2]>); Append(~sub, chain);
end for;
Lap := ZeroMatrix(Integers(), V, V);
for e in E do Lap[e[1],e[1]] +:= 1; Lap[e[2],e[2]] +:= 1; Lap[e[1],e[2]] -:= 1; Lap[e[2],e[1]] -:= 1; end for;
Sm, Pm, Qm := SmithForm(Lap);
inv := [Sm[i,i] : i in [1..V] | Sm[i,i] notin {0,1}];
printf "Phi_%o = critical group: invariants %o, order %o\n", p, inv, &*(inv cat [1]);
// placement of each CM point of the star
function Place(d)
    d0 := FundamentalDiscriminant(d); sym := KroneckerSymbol(d0, p);
    if sym eq 1 then error "split prime: ordinary reduction, not handled"; end if;
    if sym eq -1 then
        vs := { vorb[i] : i in [1..h] | OptEmb(d, RightOrder(L[i])) ne 0 };
        assert #vs eq 1;
        return <"component", Rep(vs)>;
    end if;
    es := [ e : e in [1..hp] | OptEmb(d, RightOrder(LL[e])) ne 0 ];
    orbs := { eorb[e] : e in es };
    assert #orbs eq 1;
    o := Rep(orbs);
    if o eq -1 then return <"smooth point of the merged component", vorb[org[es[1]]]>; end if;
    t := thick[o];
    assert IsEven(t);
    return <"midpoint of singular point (edge orbit)", o, "thickness", t, "vertex", t eq 2 select sub[o][1] else sub[o][t div 2]>;
end function;
pts := AssociativeArray();
for d in [-228, -760, -12, -40, -532] do
    pl := Place(d);
    pts[d] := pl[1] eq "component" select pl[2] else (pl[1] eq "smooth point of the merged component" select pl[2] else pl[5]);
    printf "P_%o: %o  -> vertex %o of the subdivided graph\n", d, pl, pts[d];
end for;
function ClassOf(dv)   // dv: list of <vertex, mult>; returns the coordinates of the class in the Smith cokernel
    v := ZeroMatrix(Integers(), V, 1);
    for t in dv do v[t[1], 1] +:= t[2]; end for;
    assert &+Eltseq(v) eq 0;
    w := Pm * v;
    return [ Sm[i,i] eq 0 select w[i,1] else (w[i,1] mod Sm[i,i]) : i in [1..V] | Sm[i,i] ne 1 ];
end function;
Delta := [<pts[-228], 1>, <pts[-760], 1>, <pts[-12], -1>, <pts[-40], -1>];
printf "class of Delta = P_-228 + P_-760 - P_-12 - P_-40 in Phi_%o: %o  (zero means trivial)\n", p, ClassOf(Delta);
for t in [ <"P_-228 - P_-760", [<pts[-228],1>,<pts[-760],-1>]>, <"P_-12 - P_-40", [<pts[-12],1>,<pts[-40],-1>]>,
           <"P_-228 + P_-12 - P_-760 - P_-40", [<pts[-228],1>,<pts[-12],1>,<pts[-760],-1>,<pts[-40],-1>]>,
           <"P_-228 + P_-40 - P_-760 - P_-12", [<pts[-228],1>,<pts[-40],1>,<pts[-760],-1>,<pts[-12],-1>]>,
           <"P_-532 - P_-40", [<pts[-532],1>,<pts[-40],-1>]>, <"P_-532 - P_-228", [<pts[-532],1>,<pts[-228],-1>]> ] do
    printf "  %o: %o\n", t[1], ClassOf(t[2]);
end for;
quit;
