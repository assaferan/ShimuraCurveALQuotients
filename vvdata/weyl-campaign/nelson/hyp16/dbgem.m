AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/nelson3.m";
D := 2002; PREC := 30; CC := ComplexField(PREC); RR := RealField(PREC); Ic := CC.1;
Wtab := WkTable(2, 22, RR);
dat := NelsonSetupT(D, 1, 2, PREC, 22, Wtab);
A := dat`A; O := dat`O; mrep := dat`mrep;
for d in [-11, -91] do
    K := QuadraticField(d);
    R := sub< MaximalOrder(K) | 1 >;      // the order of discriminant d (fundamental here)
    R := MaximalOrder(K);
    printf "\nd = %o: quadratic order discriminant %o\n", d, Discriminant(R);
    try
        mu0, f := Embed(R, O);
        mu := A ! f(R ! (K.1));   // sqrt(d): trace 0, reduced norm |d| or |d|/4
        printf "   embedded sqrt(%o): trace %o, reduced norm %o\n", d, Trace(mu), Norm(mu);
        M := ChangeRing(mrep(mu), RR);
        a := M[1][1]; b := M[1][2]; c := M[2][1];
        printf "   matrix [%o %o ; %o %o]\n", ChangeRing(Matrix(RR,1,1,[a]),RealField(8))[1][1],
            ChangeRing(Matrix(RR,1,1,[b]),RealField(8))[1][1],
            ChangeRing(Matrix(RR,1,1,[c]),RealField(8))[1][1],
            ChangeRing(Matrix(RR,1,1,[M[2][2]]),RealField(8))[1][1];
        mm := Integers() ! Norm(mu);
        z := (a + Sqrt(RR!mm)*Ic)/c; if Imaginary(z) lt 0 then z := (a - Sqrt(RR!mm)*Ic)/c; end if;
        printf "   fixed point z0 = %o   (Im = %o)\n", ChangeRing(z, ComplexField(12)),
            ChangeRing(Matrix(RR,1,1,[Imaginary(z)]),RealField(8))[1][1];
        printf "   normalises O: %o\n", forall{ b2 : b2 in ZBasis(O) | (-mu*(A!b2)*mu/mm) in O };
    catch e
        printf "   Embed failed: %o\n", e`Object;
    end try;
end for;
printf "DBGEM DONE\n"; quit;
