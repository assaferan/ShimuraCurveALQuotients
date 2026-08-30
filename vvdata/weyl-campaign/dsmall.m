AttachSpec("ShimuraQuotients.spec");
printf "START\n";
for b in [<6,1>,<10,1>,<14,1>,<15,1>,<6,5>,<22,1>,<26,1>,<6,7>] do
    Ld := ShimuraCurveLattice(b[1],b[2]);
    printf "DS %o_%o d %o\n", b[1], b[2], #Ld`disc_grp;
end for;
printf "DONE\n";
quit;
