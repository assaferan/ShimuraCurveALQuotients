// SPECIALISATION OF CM CYCLES OF THE STAR X_0^D(N)* AT p | D, GENERAL VERSION OF spec1538b.m.
// For each discriminant d in DS: the number of star points of discriminant d on each vertex orbit of the
// quotient Cerednik-Drinfeld graph (components; p inert) or on each singular point / merged component
// (p ramified), from the optimal-embedding counts (Molina).  Then the class in the critical group Phi_p
// of the divisor sum_d mult_d Z_d given by TGT, Z_d the full CM cycle of discriminant d on the star.
//   magma -b DD:=6 NN:=323 PP:=2 DS:="[-7752]" TGT:="[<-7752,0>]" tools/proof-6-323/specialise.m
AttachSpec("ShimuraQuotients.spec");
D := StringToInteger(DD); N := StringToInteger(NN); p := StringToInteger(PP);
ds := eval DS; tgt := eval TGT;
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
printf "graph data validated; unit orders vertices %o, edges %o\n", Set(uV), Set(uE);
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
printf "quotient graph: %o vertex orbits, edges %o, thicknesses %o, b1 = %o\n", nv, edges, thick, #edges - #({e[1] : e in edges} join {e[2] : e in edges}) + 1;
function OptEmb(d, Ov)
    t0 := d mod 2; n := (t0^2 - d) div 4;
    Zb := ZBasis(Ov);
    G := Matrix(Rationals(), 4, 4, [ Trace(Zb[i] * Conjugate(Zb[j])) : i, j in [1..4] ]);
    S := ShortVectors(LatticeWithGram(G), 2*n, 2*n);
    d0 := FundamentalDiscriminant(d); _, f := IsSquare(d div d0);
    opt := [];
    for s in S do for sg in [1, -1] do
        x := &+[ sg * s[1][i] * Zb[i] : i in [1..4] ];
        if Trace(x) ne t0 then continue; end if;
        if exists{<q, c> : c in [0..q-1], q in PrimeDivisors(f) | (x - c)/q in Ov} then continue; end if;
        Append(~opt, x);
    end for; end for;
    U := Units(Ov);
    return #{ { u * x * u^-1 : u in U } : x in opt };   // embeddings up to conjugation by O^*
end function;
// subdivided graph and its critical group
V := nv; E := []; sub := [];
for i->e in edges do
    chain := []; prev := e[1];
    for k in [1..thick[i]-1] do V +:= 1; Append(~chain, V); Append(~E, <prev, V>); prev := V; end for;
    Append(~E, <prev, e[2]>); Append(~sub, chain);
end for;
Lap := ZeroMatrix(Integers(), V, V);
for e in E do Lap[e[1],e[1]] +:= 1; Lap[e[2],e[2]] +:= 1; Lap[e[1],e[2]] -:= 1; Lap[e[2],e[1]] -:= 1; end for;
Sm, Pm, Qm := SmithForm(Lap);
inv := [Sm[i,i] : i in [1..V] | Sm[i,i] notin {0,1}];
printf "Phi_%o: invariants %o, order %o\n", p, inv, &*(inv cat [1]);
function Stab(d)
    e := d eq -3 select 3 else (d eq -4 select 2 else 1);
    if exists{m : m in Divisors(D*N) | m gt 1 and GCD(m, (D*N) div m) eq 1 and (d eq -4*m or (m mod 4 eq 3 and d eq -m))} then e *:= 2; end if;
    return e;
end function;
// the star points of discriminant d per vertex of the subdivided graph (an associative array vertex -> count)
function Place(d)
    d0 := FundamentalDiscriminant(d); sym := KroneckerSymbol(d0, p);
    if sym eq 1 then error "split prime: ordinary reduction, not handled"; end if;
    cnt := AssociativeArray(); total := 0;
    if sym eq -1 then
        for i in [1..h] do
            c := OptEmb(d, RightOrder(L[i]));
            if c eq 0 then continue; end if;
            for x in [i, i + h] do v := vorb[x]; cnt[v] := (IsDefined(cnt, v) select cnt[v] else 0) + c; total +:= c; end for;
        end for;
    else
        for e in [1..hp] do
            c := OptEmb(d, RightOrder(LL[e]));
            if c eq 0 then continue; end if;
            o := eorb[e];
            if o eq -1 then v := vorb[org[e]];
            else t := thick[o]; assert IsEven(t); v := t eq 2 select sub[o][1] else sub[o][t div 2]; end if;
            cnt[v] := (IsDefined(cnt, v) select cnt[v] else 0) + c; total +:= c;
        end for;
    end if;
    // CM points of X per star point: |W| / Stab(d); convert counts to star points
    ratio := #W div Stab(d);
    assert total mod ratio eq 0;
    for v in Keys(cnt) do assert cnt[v] mod ratio eq 0; cnt[v] := cnt[v] div ratio; end for;
    return cnt, total div ratio;
end function;
places := AssociativeArray();
for d in ds do
    cnt, npts := Place(d);
    places[d] := cnt;
    printf "d = %o (%o at p, %o star points): vertex -> points %o\n", d, KroneckerSymbol(FundamentalDiscriminant(d), p) eq 0 select "ramified" else "inert", npts, [ <v, cnt[v]> : v in Sort(Setseq(Keys(cnt))) ];
end for;
function ClassOf(dv)
    v := ZeroMatrix(Integers(), V, 1);
    for t in dv do for w in Keys(places[t[1]]) do v[w, 1] +:= t[2] * places[t[1]][w]; end for; end for;
    assert &+Eltseq(v) eq 0;
    w := Pm * v;
    return [ Sm[i,i] eq 0 select w[i,1] else (w[i,1] mod Sm[i,i]) : i in [1..V] | Sm[i,i] ne 1 ];
end function;
printf "class of %o in Phi_%o: %o (zero means trivial)\n", tgt, p, ClassOf(tgt);
quit;
