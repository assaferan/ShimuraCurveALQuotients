// Mass of genus(Gross(D',R)) for every s-law support triple, vs the weights.
// Triple at base (q,s,N): (q,N) w=-1/(s-1); (q,N*s) w=(s+1)/(2(s-1)); (q*s*N,1) w=-1/2.
// Bases: [q, s, N], q = LARGEST ramified prime of D (banked supports).
bases := [ [5,3,2], [7,3,2], [11,2,3], [5,2,7], [17,2,3], [7,5,2], [11,3,2], [11,5,2], [11,7,2] ];
bnames := [ "15_2", "21_2", "22_3", "10_7", "34_3", "35_2", "33_2", "55_2", "77_2" ];
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
    mass := &+[ 1/#AutomorphismGroup(Lr) : Lr in reps ];
    return Determinant(GG)/2, #reps, mass;
end function;
for bi->b in bases do
  try
    q := b[1]; s := b[2]; N := b[3];
    trip := [ <q, N>, <q, N*s>, <q*s*N, 1> ];
    ws := [ -1/(s-1), (s+1)/(2*(s-1)), -1/2 ];
    printf "== %o (q=%o, s=%o, N=%o)\n", bnames[bi], q, s, N;
    for t in [1..3] do
        Dp := trip[t][1]; Rl := trip[t][2];
        d0, ncls, mass := grossmass(Dp, Rl);
        printf "  (%o,%o): disc %o, cls %o, mass %o   w = %o   w*mass = %o   w/mass = %o\n",
            Dp, Rl, d0, ncls, mass, ws[t], ws[t]*mass, ws[t]/mass;
    end for;
  catch e
    printf "  ERROR at base %o: %o\n", bnames[bi], e`Object;
  end try;
end for;
quit;
