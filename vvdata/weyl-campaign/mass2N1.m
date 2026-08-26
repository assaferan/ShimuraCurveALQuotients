grossmass := function(Dp, Rl)
    Bq := QuaternionAlgebra(Dp);
    Oq := MaximalOrder(Bq);
    OR := Rl eq 1 select Oq else Order(Oq, Rl);
    CM := Matrix(Rationals(), 4, 4, &cat[ Eltseq(x) : x in Basis(OR) ]);
    LZ := Lattice(CM);
    Bvecs := [ Bq! Eltseq(v) : v in Basis(LZ) ];
    den0 := Lcm([Denominator(Trace(x)) : x in Bvecs]);
    TrInt := Matrix(Integers(), #Bvecs, 1, [ Integers()!(den0*Trace(x)) : x in Bvecs ]);
    KintZ := KernelMatrix(TrInt);
    S0 := [ &+[ KintZ[i][j]*Bvecs[j] : j in [1..#Bvecs] ] : i in [1..Nrows(KintZ)] ];
    GG := Matrix(Rationals(), 3,3, [ Trace(S0[i]*Conjugate(S0[j])) : j in [1..3], i in [1..3] ]);
    reps := Representatives(Genus(LatticeWithGram(GG)));
    return Determinant(GG)/2, #reps, &+[ 1/#AutomorphismGroup(Lr) : Lr in reps ];
end function;
for pr in [ <2,1>, <2,5>, <2,7>, <2,11>, <2,13> ] do
    d0, ncls, mass := grossmass(pr[1], pr[2]);
    printf "(%o,%o): disc %o, cls %o, mass %o\n", pr[1], pr[2], d0, ncls, mass;
end for;
quit;
