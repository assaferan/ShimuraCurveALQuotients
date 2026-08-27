// The 1155_1 Gross-span test (allgross.m with base-4 data).
DEPTH := 200;
MMAX := DEPTH - 1;
R<q> := PowerSeriesRing(Rationals(), DEPTH);
DN := 1155;
dsv := [1,2,3,4,5,6,7,10,11,12,14,15,20,21,22,28,30,33,35,42,44,55,60,66,70,77,84,105,110,132,140,154,165,210,220,231,308,330,385,420,462,660,770,924,1155,1540,2310,4620];
rs := [ [-2,5,0,-2,0,0,-4,0,0,0,10,0,0,0,0,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [-2,5,0,-2,0,0,0,0,-4,0,0,0,0,0,10,0,0,0,0,0,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [-1,2,0,0,0,0,-3,0,0,0,7,0,0,0,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [-1,2,0,0,0,0,-1,0,0,0,0,0,0,2,0,1,0,0,0,1,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [-1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,1,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0],
        [-1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3,0,0,0,0,7,0,0,-2,0,0,0,0],
        [0,-1,0,2,0,0,0,0,-2,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,-1,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,1,0,0,0,3,0,-1,0],
        [0,0,-2,0,0,5,0,0,-3,-2,0,0,0,0,6,0,0,1,0,0,0,0,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,-2,0,0,0,5,0,0,0,2,-2,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,0,0,-1,-1,0,-1,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,0,2,-1,0,0,-1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] ];
cs := [ -1/10, 1/9, -3/5, 8/5, -48/5, 48/5, -8/9, -12, 2/9, 3/5, 3/2, -4/3 ];
lev := 4*DN;
printf "BASE 1155_1: DN = %o, level %o\n", DN, lev;

function eta_unit(d0, e)
    ser := R!1; n := 1;
    while d0*n lt DEPTH do ser *:= (1 - q^(d0*n))^e; n +:= 1; end while;
    return ser;
end function;
function mono(rr, dsx)
    s1 := &+[dsx[i]*rr[i] : i in [1..#dsx]];
    error if s1 mod 24 ne 0, "bad monomial";
    ser := R!1;
    for i in [1..#dsx] do
        if rr[i] ne 0 then ser *:= eta_unit(dsx[i], rr[i]); end if;
    end for;
    return q^(s1 div 24) * ser;
end function;
E := &+[ cs[i]*mono(rs[i], dsv) : i in [1..#rs] ];

M := HalfIntegralWeightForms(lev, 3/2);
S := CuspidalSubspace(M);
printf "dim M = %o, dim S = %o (Eis %o)\n", Dimension(M), Dimension(S), Dimension(M)-Dimension(S);
BM := Basis(M, DEPTH);
BS := Basis(S, DEPTH);
AM := Matrix(Rationals(), MMAX+1, #BM, [ [Coefficient(BM[j], m) : j in [1..#BM]] : m in [0..MMAX] ]);
rkM := Rank(AM);
bE := Matrix(Rationals(), MMAX+1, 1, [ Coefficient(E, m) : m in [0..MMAX] ]);
printf "E in M: %o\n", Rank(HorizontalJoin(AM, bE)) eq rkM;

pairs := [ <231, 1>, <231, 5> ];
printf "definite structures: %o\n", pairs;

cols := [* *]; names := [];
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
    reps := Representatives(Genus(LGr));
    num := [ Rationals()| 0 : m in [0..MMAX] ]; den := 0;
    for Lr in reps do
        wt := 1/#AutomorphismGroup(Lr);
        T := ThetaSeries(Lr, 2*MMAX+1);
        for m in [0..MMAX] do num[m+1] +:= wt*Coefficient(T, 2*m); end for;
        den +:= wt;
    end for;
    Append(~cols, [ x/den : x in num ]);
    Append(~names, Sprintf("GROSS(%o,%o) disc %o, %o cls", pr[1], pr[2], Determinant(GG)/2, #reps));
    printf "  built %o\n", names[#names];
end for;
ngross := #cols;
for j in [1..#BS] do Append(~cols, [ Coefficient(BS[j], m) : m in [0..MMAX] ]); Append(~names, Sprintf("CUSP_%o", j)); end for;
nc := #cols;
A := Matrix(Rationals(), MMAX+1, nc, [ [cols[j][m+1] : j in [1..nc]] : m in [0..MMAX] ]);
v := Vector(Rationals(), [ Coefficient(E, m) : m in [0..MMAX] ]);
rA := Rank(A); rAug := Rank(HorizontalJoin(A, Matrix(MMAX+1,1,Eltseq(v))));
printf "Gross-average+cusp span: rank %o vs %o -> E in span: %o\n", rA, rAug, rA eq rAug;
if rA eq rAug then
    x := Solution(Transpose(A), v);
    printf "E_eis =\n";
    for j in [1..nc] do if x[j] ne 0 then printf "  %8o * %o\n", x[j], names[j]; end if; end for;
    K := Kernel(Transpose(A));
    printf "aliasing dim %o\n", Dimension(K);
end if;
quit;
