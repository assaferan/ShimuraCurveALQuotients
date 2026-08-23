// Character test: does (panel monomial)*(E-pool monomial) have trivial
// multiplier on Gamma_0(60)?  Ratio P(g tau)/((c tau + d)^2 P(tau)) for
// generators-ish g with various d mod 60.  Control: product of two panel
// monomials (weight 1) with the (s/d) character, s = square => trivial.
CC := ComplexField(50); ii := CC.1;
ds := [1,2,3,4,5,6,10,12,15,20,30,60];
rf := [-21, 8, 3, 1, 20, 8, -9, -4, -6, 0, 2, -1];
rE := [-3,6,-1,0,0,3,0,-2,0,0,0,0];
tau := CC!0.13 + CC!0.9*ii;
P := func< r, z | &*[ CC | DedekindEta(d*z)^(r[i]) : i->d in ds | r[i] ne 0 ] >;
gs := [ [1,0,60,1], [7,5,60,43], [11,2,60,11], [13,8,60,37], [59,58,60,59], [-1,0,0,-1] ];
for g in gs do
    a,b,c,d := Explode(g);
    error if a*d-b*c ne 1, "det";
    gz := (a*tau+b)/(c*tau+d);
    for pair in [ <"f*E", [ rf[i]+rE[i] : i in [1..#ds] ], 2>,
                  <"f*f", [ 2*rf[i] : i in [1..#ds] ], 1> ] do
        nm, r, wt := Explode(pair);
        rat := P(r, gz) / ( (c*tau+d)^wt * P(r, tau) );
        printf "%o gamma d=%o: ratio = %o + %o i\n", nm, d mod 60,
            RealField(8)!Re(rat), RealField(8)!Im(rat);
    end for;
end for;
quit;
