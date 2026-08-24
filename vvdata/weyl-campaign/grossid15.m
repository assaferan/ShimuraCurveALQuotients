// Match the four definite Gross lattices to the four disc-900 member genera by theta comparison.
DEPTH := 60;   // compare 60 coefficients: ample for genus identification
MM := DEPTH - 1;
grams := [
    [1,15,15, 0,0,0], [2,8,15,  0,0,2], [3,10,10, 10,0,0], [5,6,9,   6,0,0],
    [4,4,15,  0,0,2], [3,5,15,  0,0,0]
];
Ls := [ LatticeWithGram(Matrix(Rationals(), 3,3,
        [2*g[1], g[6], g[5],  g[6], 2*g[2], g[4],  g[5], g[4], 2*g[3]])) : g in grams ];
// genus averages for the four disc-900 genera:
// G6 = {[1,15,15],[4,4,15|2]}, G10 = {[2,8,15|2],[3,5,15]}, G12 = {[3,10,10|10]}, G13 = {[5,6,9|6]}
gsets := [ [1,5], [2,6], [3], [4] ];
gnames := [ "G6{1,15,15 / 4,4,15|2}", "G10{2,8,15|2 / 3,5,15}", "G12{3,10,10|10}", "G13{5,6,9|6}" ];
avgs := [* *];
for gs in gsets do
    num := [ Rationals()| 0 : m in [0..MM] ]; den := 0;
    for i in gs do
        wt := 1/#AutomorphismGroup(Ls[i]);
        T := ThetaSeries(Ls[i], 4*MM+2);
        for m in [0..MM] do num[m+1] +:= wt*Coefficient(T, 2*m); end for;
        den +:= wt;
    end for;
    Append(~avgs, [ x/den : x in num ]);
end for;

for pairDR in [ <30,1>, <2,15>, <3,10>, <5,6> ] do
    Bq := QuaternionAlgebra(pairDR[1]);
    Oq := MaximalOrder(Bq);
    OR := pairDR[2] eq 1 select Oq else Order(Oq, pairDR[2]);
    bas := Basis(OR);
    for variant in [1, 2] do
        gens := variant eq 1 select ([ Bq!1 ] cat [ 2*b : b in bas ]) else [ b : b in bas ];
        CM := Matrix(Rationals(), #gens, 4, &cat[ Eltseq(x) : x in gens ]);
        LZ := Lattice(CM);
        Bvecs := [ Bq! Eltseq(v) : v in Basis(LZ) ];
        den0 := Lcm([Denominator(Trace(x)) : x in Bvecs]);
        TrInt := Matrix(Integers(), #Bvecs, 1, [ Integers()!(den0*Trace(x)) : x in Bvecs ]);
        KintZ := KernelMatrix(TrInt);
        S0 := [ &+[ KintZ[i][j]*Bvecs[j] : j in [1..#Bvecs] ] : i in [1..Nrows(KintZ)] ];
        GG := Matrix(Rationals(), 3,3, [ Trace(S0[i]*Conjugate(S0[j])) : j in [1..3], i in [1..3] ]);
        LGr := LatticeWithGram(GG);   // norm x.GG.x = 2 Nrd
        TG := ThetaSeries(LGr, 8*MM+4);
        // b(m) := # {Nrd = m}
        bseq := [ Coefficient(TG, 2*m) : m in [0..MM] ];       // Nrd values
        b2 :=  [ Coefficient(TG, 4*m) : m in [0..MM] ];        // Nrd = 2m (form Nrd/2)
        for gi in [1..4] do
            a := avgs[gi];
            if forall{m : m in [0..MM] | a[m+1] eq bseq[m+1]} then
                printf "(D'=%o,R=%o) v%o  Nrd    == avg %o\n", pairDR[1], pairDR[2], variant, gnames[gi];
            end if;
            if forall{m : m in [0..MM] | a[m+1] eq b2[m+1]} then
                printf "(D'=%o,R=%o) v%o  Nrd/2  == avg %o\n", pairDR[1], pairDR[2], variant, gnames[gi];
            end if;
        end for;
    end for;
end for;
print "done";
quit;
