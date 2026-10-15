load "/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/vw6/vw2.m";
SetColumns(0);
PREC := 60; CC := ComplexField(PREC);
A := QuaternionAlgebra(6); O := QuaternionOrder(A, 1);
G := FuchsianGroup(O); DD := UnitDisc( : Precision := PREC);
FD := FundamentalDomain(G, DD);
rho := Max([Abs(ComplexValue(v)) : v in FD]);
H := UpperHalfPlane();
p := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);
gens := Generators(G);
Q := 40; pi := Pi(RealField(PREC));
nin := 0; nout := 0; mods := [];
for m in [1..Q] do
    wc := CC ! (rho * Exp(2*pi*CC.1*m/Q));
    wpt := DD ! wc;
    red := ShimuraReduceUnit(G ! (O ! 1), gens, DD : z0 := wpt);
    delta := red[1];
    if IsScalar(Quaternion(delta)) then nin +:= 1; continue; end if;
    nout +:= 1;
    W := DiscMatrix(mrep, Quaternion(delta), CC ! p, CC);
    wp := (W[1][1]*wc + W[1][2]) / (W[2][1]*wc + W[2][2]);
    // apply delta via Magma's own action, to compare with my disc matrix
    wp2 := ComplexValue(delta * wpt);
    // is wp back in the domain?  re-reduce it
    red2 := ShimuraReduceUnit(G ! (O ! 1), gens, DD : z0 := DD ! wp);
    inD := IsScalar(Quaternion(red2[1]));
    Append(~mods, Abs(wp));
    if nout le 6 then
        printf "m=%o |w|=%o |w'|=%o  mine-vs-Magma %o  w' in D: %o\n", m, RealField(4)!Abs(wc),
           RealField(5)!Abs(wp), RealField(3)!Abs(wp - wp2), inD;
    end if;
end for;
printf "already in D: %o, reduced: %o\n", nin, nout;
printf "|w'| min %o mean %o max %o (rho = %o)\n", RealField(5)!Min(mods),
   RealField(5)!(&+mods/#mods), RealField(5)!Max(mods), RealField(5)!rho;
printf "DBGC DONE\n";
quit;
