AttachSpec("ShimuraQuotients.spec");
printf "START\n";
for b in [<34,3>,<38,5>,<38,7>] do
    D := b[1]; N := b[2];
    t := Realtime();
    Ld := ShimuraCurveLattice(D,N);
    dg := Ld`disc_grp;
    printf "DSIZE %o_%o discgrp %o det %o  (%os)\n", D, N, #dg,
           Determinant(ChangeRing(Ld`Q, Rationals())), RealField(6)!(Realtime()-t);
end for;
printf "DONE\n";
quit;
