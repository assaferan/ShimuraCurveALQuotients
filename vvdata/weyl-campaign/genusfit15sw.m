// Canonical form: decompose E_eis over GENUS AVERAGES (Siegel-Weil Eisenstein basis) + cusp.
// Also identify the Gross lattice of the (15,2) Eichler order among the genera.
DEPTH := 400;
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
printf "%o lattices fall into %o genera\n", #Ls, ng;

avgs := [* *]; gnames := [];
for g in [1..ng] do
    idx := [ i : i in [1..#Ls] | genusof[i] eq g ];
    reps := Representatives(Genus(Ls[idx[1]]));
    d0 := Determinant(GramMatrix(Ls[idx[1]]))/2;
    printf "  genus %o: disc %o, members %o, classes in genus %o\n", g, d0, idx, #reps;
    num := [ Rationals()| 0 : m in [0..MMAX] ];
    den := 0;
    for Lr in reps do
        wt := 1/#AutomorphismGroup(Lr);
        T := ThetaSeries(Lr, 2*MMAX+1);
        for m in [0..MMAX] do num[m+1] +:= wt*Coefficient(T, 2*m); end for;
        den +:= wt;
    end for;
    Append(~avgs, [ x/den : x in num ]);
    Append(~gnames, Sprintf("GENUS%o(disc %o)", g, d0));
end for;

// (Gross-lattice identification excised -- answered by grossid100/grossid900,
// banked 10e4bdb + 92b3039; the Z+2O construction here was broken anyway.)

// decompose E over genus averages + cusp
M := HalfIntegralWeightForms(60, 3/2);
BS := Basis(CuspidalSubspace(M), DEPTH);
cols := [* *]; names := [];
for g in [1..ng] do Append(~cols, avgs[g]); Append(~names, gnames[g]); end for;
for j in [1..#BS] do Append(~cols, [ Coefficient(BS[j], m) : m in [0..MMAX] ]); Append(~names, Sprintf("CUSP_%o", j)); end for;
nc := #cols;
A := Matrix(Rationals(), MMAX+1, nc, [ [cols[j][m+1] : j in [1..nc]] : m in [0..MMAX] ]);
v := Vector(Rationals(), [ Coefficient(E, m) : m in [0..MMAX] ]);
rA := Rank(A); rAug := Rank(HorizontalJoin(A, Matrix(MMAX+1,1,Eltseq(v))));
printf "genus-average span rank %o vs %o; E in it: %o\n", rA, rAug, rA eq rAug;
if rA eq rAug then
    x := Solution(Transpose(A), v);
    printf "E_eis over the Siegel-Weil basis:\n";
    for j in [1..nc] do if x[j] ne 0 then printf "  %8o * %o\n", x[j], names[j]; end if; end for;
    K := Kernel(Transpose(A));
    printf "aliasing dim %o\n", Dimension(K);
end if;
quit;
