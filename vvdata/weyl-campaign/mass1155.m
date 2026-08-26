// masses of the candidate supports at DN = 1155, computed IN ADVANCE
prs := [ <385,1>, <385,3>, <11,1>, <11,105>, <105,1>, <105,11>,
         <165,1>, <165,7>, <231,1>, <231,5> ];
for pr in prs do
    t0 := Cputime();
    Dp := pr[1]; Rl := pr[2];
    Bq := QuaternionAlgebra(Dp);
    Oq := MaximalOrder(Bq);
    OR := Rl eq 1 select Oq else Order(Oq, Rl);
    CM := Matrix(Rationals(), 4, 4, &cat[ Eltseq(x) : x in Basis(OR) ]);
    LZ := Lattice(CM);
    Bv := [ Bq! Eltseq(v) : v in Basis(LZ) ];
    den0 := Lcm([Denominator(Trace(x)) : x in Bv]);
    TrInt := Matrix(Integers(), #Bv, 1, [ Integers()!(den0*Trace(x)) : x in Bv ]);
    KintZ := KernelMatrix(TrInt);
    S0 := [ &+[ KintZ[i][j]*Bv[j] : j in [1..#Bv] ] : i in [1..Nrows(KintZ)] ];
    GG := Matrix(Rationals(), 3,3, [ Trace(S0[i]*Conjugate(S0[j])) : j in [1..3], i in [1..3] ]);
    LGr := LatticeWithGram(GG);
    reps := Representatives(Genus(LGr));
    mass := &+[ Rationals() | 1/#AutomorphismGroup(Lr) : Lr in reps ];
    printf "MASS GROSS(%o,%o): disc %o, %o cls, mass = %o   (%o s)\n",
        Dp, Rl, Determinant(GG)/2, #reps, mass, RealField(5)!Cputime(t0);
end for;
quit;
