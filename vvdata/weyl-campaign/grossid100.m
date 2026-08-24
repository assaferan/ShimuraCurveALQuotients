// Which definite structure gives the disc-100 genera G5 ([1,5,5]) and G9 ([2,3,5|xy=2])?
gramsL := [ [1,5,5, 0,0,0], [2,3,5, 0,0,2], [1,3,3, 0,0,0], [2,2,3, 0,0,2], [1,1,3, 0,0,1], [1,1,1, 0,0,0], [1,4,15,0,0,1], [2,2,15,0,0,1], [3,5,5,5,0,0], [3,4,4,4,0,0] ];
Lnames := [ "G5[1,5,5](d100)", "G9[2,3,5|2](d100)", "G3[1,3,3](d36)", "G7[2,2,3|2](d36)", "G2[1,1,3|1](d9)", "G1[1,1,1](d4)", "G4[1,4,15|1](d225)", "G8[2,2,15|1](d225)", "G8b[3,5,5|5](d225)", "G11[3,4,4|4](d144)" ];
Ls := [ LatticeWithGram(Matrix(Rationals(), 3,3,
        [2*g[1], g[6], g[5],  g[6], 2*g[2], g[4],  g[5], g[4], 2*g[3]])) : g in gramsL ];
genera := [ Genus(L) : L in Ls ];
for pairDR in [ <2,5>, <5,2>, <2,3>, <3,2>, <3,1>, <5,1>, <2,1>, <30,2>, <2,15>, <5,3>, <3,5>, <2,2>, <2,4>, <2,6>, <2,10>, <3,4>, <5,4> ] do
    Bq := QuaternionAlgebra(pairDR[1]);
    Oq := MaximalOrder(Bq);
    OR := pairDR[2] eq 1 select Oq else Order(Oq, pairDR[2]);
    CM := Matrix(Rationals(), 4, 4, &cat[ Eltseq(x) : x in Basis(OR) ]);
    LZ := Lattice(CM);
    Bvecs := [ Bq! Eltseq(v) : v in Basis(LZ) ];
    den0 := Lcm([Denominator(Trace(x)) : x in Bvecs]);
    TrInt := Matrix(Integers(), #Bvecs, 1, [ Integers()!(den0*Trace(x)) : x in Bvecs ]);
    KintZ := KernelMatrix(TrInt);
    S0 := [ &+[ KintZ[i][j]*Bvecs[j] : j in [1..#Bvecs] ] : i in [1..Nrows(KintZ)] ];
    GG := Matrix(Rationals(), 3,3, [ Trace(S0[i]*Conjugate(S0[j])) : j in [1..3], i in [1..3] ]);
    GLG := Genus(LatticeWithGram(GG));
    hit := "-(disc " cat Sprintf("%o", Determinant(GG)/2) cat ")";
    for g in [1..#Ls] do
        if GLG eq genera[g] then hit := Lnames[g]; end if;
    end for;
    printf "Gross(D'=%o,R=%o) -> %o\n", pairDR[1], pairDR[2], hit;
end for;
quit;
