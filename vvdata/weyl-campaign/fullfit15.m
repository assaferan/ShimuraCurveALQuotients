// Explicit decomposition of E_eis over the 16 member ternary thetas + cusp basis.
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
names := []; cols := [* *];
for g in grams do
    G := Matrix(Rationals(), 3,3, [2*g[1], g[6], g[5],  g[6], 2*g[2], g[4],  g[5], g[4], 2*g[3]]);
    L := LatticeWithGram(G);
    T := ThetaSeries(L, 2*MMAX+1);
    Append(~cols, [Rationals()| Coefficient(T, 2*m) : m in [0..MMAX] ]);
    Append(~names, Sprintf("disc%o[%o,%o,%o|%o,%o,%o]", Determinant(G)/2, g[1],g[2],g[3],g[4],g[5],g[6]));
end for;
M := HalfIntegralWeightForms(60, 3/2);
BS := Basis(CuspidalSubspace(M), DEPTH);
for j in [1..#BS] do Append(~cols, [ Coefficient(BS[j], m) : m in [0..MMAX] ]); Append(~names, Sprintf("CUSP_%o", j)); end for;
nc := #cols;
A := Matrix(Rationals(), MMAX+1, nc, [ [cols[j][m+1] : j in [1..nc]] : m in [0..MMAX] ]);
v := Vector(Rationals(), [ Coefficient(E, m) : m in [0..MMAX] ]);
bcol := Matrix(Rationals(), MMAX+1, 1, Eltseq(v));
rA := Rank(A); rAug := Rank(HorizontalJoin(A, bcol));
printf "rank %o vs %o; E in span: %o\n", rA, rAug, rA eq rAug;
x := Solution(Transpose(A), v);
printf "one representative decomposition:\n";
for j in [1..nc] do if x[j] ne 0 then printf "  %8o * %o\n", x[j], names[j]; end if; end for;
// kernel: aliasing among thetas (theta linear relations)
K := Kernel(Transpose(A));
printf "aliasing kernel dim %o\n", Dimension(K);
for k in Basis(K) do printf "  reln: %o\n", [k[j] : j in [1..nc]]; end for;
quit;
