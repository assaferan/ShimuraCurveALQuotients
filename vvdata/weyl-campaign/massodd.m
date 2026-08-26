// masses of every definite structure (D',R) with D'R | DN, for the odd-DN N=1 bases
for DN in [15, 21, 33, 35] do
    printf "== DN = %o ==\n", DN;
    for Dp in [ d : d in Divisors(DN) | d gt 1 and IsSquarefree(d) and IsOdd(#PrimeDivisors(d)) ] do
        for Rl in Divisors(DN div Dp) do
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
            printf "  GROSS(%o,%o): disc %o, %o cls, mass = %o\n",
                Dp, Rl, Determinant(GG)/2, #reps, mass;
        end for;
    end for;
end for;
quit;
