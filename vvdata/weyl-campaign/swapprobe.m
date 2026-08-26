// Which 900-genus do the (5,6) and (30,1) trace-zero lattices REALLY lie in,
// and how many classes does each have? (Resolves the fit-vs-dictionary swap.)
gramsL := [ [1,15,15,0,0,0], [4,4,15,0,0,2], [2,8,15,0,0,2], [3,10,10,10,0,0], [3,5,15,0,0,0], [5,6,9,6,0,0] ];
Lnames := [ "G6a[1,15,15]", "G6b[4,4,15|2]", "G10a[2,8,15|2]", "G10b[3,10,10|10]", "G12[3,5,15]", "G13[5,6,9|6]" ];
Ls := [ LatticeWithGram(Matrix(Rationals(), 3,3,
        [2*g[1], g[6], g[5],  g[6], 2*g[2], g[4],  g[5], g[4], 2*g[3]])) : g in gramsL ];
genera := [ Genus(L) : L in Ls ];
printf "reference class counts: %o\n", [ #Representatives(g) : g in genera ];
for trial in [1..3] do
for pairDR in [ <5,6>, <30,1>, <2,15>, <3,10> ] do
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
    LGr := LatticeWithGram(GG);
    GLG := Genus(LGr);
    nreps := #Representatives(GLG);
    hit := "NONE";
    for g in [1..#Ls] do if GLG eq genera[g] then hit := Lnames[g]; end if; end for;
    printf "trial %o: (D'=%o,R=%o) disc %o, %o cls -> %o\n", trial, pairDR[1], pairDR[2], Determinant(GG)/2, nreps, hit;
end for;
end for;
quit;
