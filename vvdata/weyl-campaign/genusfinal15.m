// FINAL: (a) confirm the (3,10)/(5,6) Gross genera; (b) decompose E_eis over the
// 13 genus averages (Siegel-Weil Eisenstein basis) + cusp forms of M_{3/2}(60).
DEPTH := 200;
MMAX := DEPTH - 1;
R<q> := PowerSeriesRing(Rationals(), DEPTH);
function eta_unit(d0, e)
    ser := R!1; n := 1;
    while d0*n lt DEPTH do ser *:= (1 - q^(d0*n))^e; n +:= 1; end while;
    return ser;
end function;
ds := [1,2,3,4,5,6,10,12,15,20,30,60];
function mono(rr)
    s1 := &+[ds[i]*rr[i] : i in [1..#ds]];
    ser := R!1;
    for i in [1..#ds] do
        if rr[i] ne 0 then ser *:= eta_unit(ds[i], rr[i]); end if;
    end for;
    return q^(s1 div 24) * ser;
end function;
E := (4/5)*mono([-6,15,0,-6,0,0,0,0,0,0,0,0]) - 4*mono([-3,7,0,-3,1,0,0,0,0,1,0,0])
     - (4/5)*mono([-1,0,0,4,-1,0,3,0,0,-2,0,0]) + 4*mono([-1,2,0,0,0,0,0,0,-3,0,7,-2]);

grams := [
    [1,1,1,   0,0,0], [1,1,3,   0,0,1], [1,3,3,   0,0,0], [1,4,15,  0,0,1],
    [1,5,5,   0,0,0], [1,15,15, 0,0,0], [2,2,3,   0,0,2], [2,2,15,  0,0,1],
    [2,3,5,   0,0,2], [2,8,15,  0,0,2], [3,4,4,   4,0,0], [3,5,5,   5,0,0],
    [3,5,15,  0,0,0], [3,10,10, 10,0,0], [4,4,15, 0,0,2], [5,6,9,   6,0,0]
];
Ls := [ LatticeWithGram(Matrix(Rationals(), 3,3,
        [2*g[1], g[6], g[5],  g[6], 2*g[2], g[4],  g[5], g[4], 2*g[3]])) : g in grams ];
genera := [ Genus(L) : L in Ls ];
genusof := [0 : L in Ls];
ng := 0;
for i in [1..#Ls] do
    if genusof[i] ne 0 then continue; end if;
    ng +:= 1; genusof[i] := ng;
    for j in [i+1..#Ls] do
        if genusof[j] eq 0 and genera[i] eq genera[j] then genusof[j] := ng; end if;
    end for;
end for;

// (a) all four definite Gross lattices (trace-zero of O_R, form Nrd)
for pairDR in [ <30,1>, <2,15>, <3,10>, <5,6> ] do
    Bq := QuaternionAlgebra(pairDR[1]);
    Oq := MaximalOrder(Bq);
    OR := pairDR[2] eq 1 select Oq else Order(Oq, pairDR[2]);
    Bvs := Basis(OR);
    CM := Matrix(Rationals(), 4, 4, &cat[ Eltseq(x) : x in Bvs ]);
    LZ := Lattice(CM);
    Bvecs := [ Bq! Eltseq(v) : v in Basis(LZ) ];
    den0 := Lcm([Denominator(Trace(x)) : x in Bvecs]);
    TrInt := Matrix(Integers(), #Bvecs, 1, [ Integers()!(den0*Trace(x)) : x in Bvecs ]);
    KintZ := KernelMatrix(TrInt);
    S0 := [ &+[ KintZ[i][j]*Bvecs[j] : j in [1..#Bvecs] ] : i in [1..Nrows(KintZ)] ];
    GG := Matrix(Rationals(), 3,3, [ Trace(S0[i]*Conjugate(S0[j])) : j in [1..3], i in [1..3] ]);
    LGr := LatticeWithGram(GG);
    GLG := Genus(LGr);
    hit := false;
    for g in [1..ng] do
        if GLG eq genera[Index(genusof, g)] then
            printf "Gross(D'=%o,R=%o) is GENUS %o (disc %o)\n", pairDR[1], pairDR[2], g,
                Determinant(GramMatrix(Ls[Index(genusof,g)]))/2;
            hit := true;
        end if;
    end for;
    if not hit then printf "Gross(D'=%o,R=%o): disc %o, NOT among the 13\n", pairDR[1], pairDR[2], Determinant(GG)/2; end if;
end for;

// (b) genus averages and the decomposition
avgs := [* *];
for g in [1..ng] do
    idx := [ i : i in [1..#Ls] | genusof[i] eq g ];
    reps := Representatives(Genus(Ls[idx[1]]));
    num := [ Rationals()| 0 : m in [0..MMAX] ]; den := 0;
    for Lr in reps do
        wt := 1/#AutomorphismGroup(Lr);
        T := ThetaSeries(Lr, 2*MMAX+1);
        for m in [0..MMAX] do num[m+1] +:= wt*Coefficient(T, 2*m); end for;
        den +:= wt;
    end for;
    Append(~avgs, [ x/den : x in num ]);
end for;
M := HalfIntegralWeightForms(60, 3/2);
BS := Basis(CuspidalSubspace(M), DEPTH);
cols := [* *]; names := [];
for g in [1..ng] do
    Append(~cols, avgs[g]);
    Append(~names, Sprintf("AVG_G%o(disc %o)", g, Determinant(GramMatrix(Ls[Index(genusof,g)]))/2));
end for;
for j in [1..#BS] do Append(~cols, [ Coefficient(BS[j], m) : m in [0..MMAX] ]); Append(~names, Sprintf("CUSP_%o", j)); end for;
nc := #cols;
A := Matrix(Rationals(), MMAX+1, nc, [ [cols[j][m+1] : j in [1..nc]] : m in [0..MMAX] ]);
v := Vector(Rationals(), [ Coefficient(E, m) : m in [0..MMAX] ]);
rA := Rank(A); rAug := Rank(HorizontalJoin(A, Matrix(MMAX+1,1,Eltseq(v))));
printf "genus-average+cusp span rank %o vs %o; E in it: %o\n", rA, rAug, rA eq rAug;
if rA eq rAug then
    x := Solution(Transpose(A), v);
    printf "E_eis over the Siegel-Weil basis:\n";
    for j in [1..nc] do if x[j] ne 0 then printf "  %8o * %o\n", x[j], names[j]; end if; end for;
    K := Kernel(Transpose(A));
    printf "aliasing dim %o\n", Dimension(K);
    for k in Basis(K) do printf "  reln: %o\n", [k[j] : j in [1..nc]]; end for;
end if;
quit;
