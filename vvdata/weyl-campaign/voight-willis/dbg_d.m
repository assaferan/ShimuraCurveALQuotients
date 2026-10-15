load "/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/vw6/vw2.m";
SetColumns(0);
PREC := 100; CC := ComplexField(PREC);
A := QuaternionAlgebra(6); O := QuaternionOrder(A, 1);
G := FuchsianGroup(O); DD := UnitDisc( : Precision := PREC);
FD := FundamentalDomain(G, DD);
rho := Max([Abs(ComplexValue(v)) : v in FD]);
H := UpperHalfPlane();
p := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);
samp := ReducedSamples(G, O, DD, mrep, p, rho, 300, CC);
printf "rho %o, %o samples\n", RealField(6)!rho, #samp;
for k in [4, 12] do
    printf "\n--- k = %o (expect dim %o)\n", k, k eq 4 select 1 else 3;
    for NT in [15, 25, 35, 45, 60] do
        sv := SingVals(NormaliseColumns(AutSystem(samp, k, NT, CC)));
        printf "  N=%2o rho^N=%o  sv[1..6] %o\n", NT, RealField(3)!(rho^NT),
            [RealField(3)!s : s in sv[1..Min(6,#sv)]];
    end for;
end for;
printf "DBGD DONE\n";
quit;
