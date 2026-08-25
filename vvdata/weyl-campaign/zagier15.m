// w_sq readout at 15_2: is the Gross-basis E_eis expressible over ZAGIER
// PULLBACKS V_d E_Z (Hurwitz series, 4d | 120) + cusp?  Which Gross averages
// live in that (square-class-trivial) channel?  Plus genus identification
// among the 13 Gross lattices (the aliasing structure).
DEPTH := 200;
MMAX := DEPTH - 1;
R<q> := PowerSeriesRing(Rationals(), DEPTH);
DN := 30; lev := 120;
function eta_unit(d0, e)
    ser := R!1; n := 1;
    while d0*n lt DEPTH do ser *:= (1 - q^(d0*n))^e; n +:= 1; end while;
    return ser;
end function;
ds60 := [1,2,3,4,5,6,10,12,15,20,30,60];
function mono(rr)
    s1 := &+[ds60[i]*rr[i] : i in [1..#ds60]];
    ser := R!1;
    for i in [1..#ds60] do
        if rr[i] ne 0 then ser *:= eta_unit(ds60[i], rr[i]); end if;
    end for;
    return q^(s1 div 24) * ser;
end function;
E := (4/5)*mono([-6,15,0,-6,0,0,0,0,0,0,0,0]) - 4*mono([-3,7,0,-3,1,0,0,0,0,1,0,0])
     - (4/5)*mono([-1,0,0,4,-1,0,3,0,0,-2,0,0]) + 4*mono([-1,2,0,0,0,0,0,0,-3,0,7,-2]);

M := HalfIntegralWeightForms(lev, 3/2);
S := CuspidalSubspace(M);
BS := Basis(S, DEPTH);
printf "dim M = %o, dim S = %o (Eis %o)\n", Dimension(M), Dimension(S), Dimension(M)-Dimension(S);

// the 13 Gross averages (construction from the validated allgross driver)
pairs := [];
for Dp in [ d : d in Divisors(DN) | d gt 1 and IsSquarefree(d) and IsOdd(#PrimeDivisors(d)) ] do
    for Rl in Divisors(DN div Dp) do Append(~pairs, <Dp, Rl>); end for;
end for;
gcols := [* *]; gnames := []; glats := [* *];
for pr in pairs do
    Bq := QuaternionAlgebra(pr[1]);
    Oq := MaximalOrder(Bq);
    OR := pr[2] eq 1 select Oq else Order(Oq, pr[2]);
    CM := Matrix(Rationals(), 4, 4, &cat[ Eltseq(x) : x in Basis(OR) ]);
    LZ := Lattice(CM);
    Bvecs := [ Bq! Eltseq(v) : v in Basis(LZ) ];
    den0 := Lcm([Denominator(Trace(x)) : x in Bvecs]);
    TrInt := Matrix(Integers(), #Bvecs, 1, [ Integers()!(den0*Trace(x)) : x in Bvecs ]);
    KintZ := KernelMatrix(TrInt);
    S0 := [ &+[ KintZ[i][j]*Bvecs[j] : j in [1..#Bvecs] ] : i in [1..Nrows(KintZ)] ];
    GG := Matrix(Rationals(), 3,3, [ Trace(S0[i]*Conjugate(S0[j])) : j in [1..3], i in [1..3] ]);
    LGr := LatticeWithGram(GG);
    Append(~glats, LGr);
    reps := Representatives(Genus(LGr));
    num := [ Rationals()| 0 : m in [0..MMAX] ]; den := 0;
    for Lr in reps do
        wt := 1/#AutomorphismGroup(Lr);
        T := ThetaSeries(Lr, 2*MMAX+1);
        for m in [0..MMAX] do num[m+1] +:= wt*Coefficient(T, 2*m); end for;
        den +:= wt;
    end for;
    Append(~gcols, [ x/den : x in num ]);
    Append(~gnames, Sprintf("GROSS(%o,%o)", pr[1], pr[2]));
end for;

// genus identification among the 13
printf "genus classes among the 13 Gross lattices:\n";
gen := [* Genus(L) : L in glats *];
done := [ false : p in pairs ];
for i in [1..#pairs] do
    if done[i] then continue; end if;
    cls := [ i ];
    for j in [i+1..#pairs] do
        if not done[j] and gen[i] eq gen[j] then Append(~cls, j); done[j] := true; end if;
    end for;
    printf "  { %o }  disc %o\n", Join([gnames[c] : c in cls], ", "), Determinant(GramMatrix(glats[i]))/2;
end for;

// Zagier pullbacks V_d E_Z, 4d | lev
function H24(k)
    // Hurwitz class number H(k): sum over f^2 | k of h(-k/f^2)/(w/2)
    if k eq 0 then return Rationals()!(-1/12); end if;
    if k mod 4 in {1,2} then return Rationals()!0; end if;
    tot := Rationals()!0;
    for f in [1..Isqrt(k)] do
        if k mod (f*f) ne 0 then continue; end if;
        d := -(k div (f*f));
        if ((d mod 4) + 4) mod 4 notin {0,1} then continue; end if;
        w := d eq -3 select 1/3 else (d eq -4 select 1/2 else 1);
        tot +:= w*ClassNumber(d);
    end for;
    return tot;
end function;
zds := [ d : d in Divisors(DN) ];
zcols := [* *]; znames := [];
for d in zds do
    Append(~zcols, [ (m mod d eq 0) select H24(m div d) else 0 : m in [0..MMAX] ]);
    Append(~znames, Sprintf("Z_%o", d));
end for;

function spantest(cols, names, vec)
    nc := #cols;
    A := Matrix(Rationals(), MMAX+1, nc, [ [cols[j][m+1] : j in [1..nc]] : m in [0..MMAX] ]);
    rA := Rank(A);
    rAug := Rank(HorizontalJoin(A, Matrix(MMAX+1,1,vec)));
    return rA, rAug eq rA, A;
end function;

vE := [ Coefficient(E, m) : m in [0..MMAX] ];
cuspcols := [* *];
for j in [1..#BS] do Append(~cuspcols, [ Coefficient(BS[j], m) : m in [0..MMAX] ]); end for;

// 1) Zagier alone and Zagier+cusp
zc := zcols cat cuspcols; zcn := znames cat [ Sprintf("CUSP_%o", j) : j in [1..#BS] ];
rZ, _, _ := spantest(zcols, znames, vE);
rZC, EinZC, A := spantest(zc, zcn, vE);
printf "rank(Zagier) = %o (of %o), rank(Zagier+cusp) = %o, E in Zagier+cusp: %o\n", rZ, #zds, rZC, EinZC;
if EinZC then
    x := Solution(Transpose(A), Vector(Rationals(), vE));
    printf "E_eis over Zagier+cusp:\n";
    for j in [1..#zcn] do if x[j] ne 0 then printf "  %8o * %o\n", x[j], zcn[j]; end if; end for;
end if;
// 2) which Gross averages live in Zagier+cusp?
for i in [1..#pairs] do
    _, inz, _ := spantest(zc, zcn, [ Rationals()!g : g in gcols[i] ]);
    printf "  %o in Zagier+cusp: %o\n", gnames[i], inz;
end for;
// 3) combined dictionary size
allc := gcols cat zc; alln := gnames cat zcn;
rAll, _, _ := spantest(allc, alln, vE);
printf "rank(Gross+Zagier+cusp) = %o (Gross+cusp alone was 14)\n", rAll;
quit;
