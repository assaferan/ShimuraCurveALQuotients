// Which part of NelsonSetupT differs between processes: the algebra's standard generators,
// the Eichler order, or the real splitting?  Run several times and diff.
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
D := StringToInteger(DISC); N := StringToInteger(LEV);
A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
a, b := StandardForm(A);
printf "standard form a, b = %o, %o\n", a, b;
printf "A generators i^2, j^2 = %o, %o\n", (A.1)^2, (A.2)^2;
printf "ZBasis(O) = %o\n", [ Eltseq(A ! x) : x in ZBasis(O) ];
mrep := FuchsianMatrixRepresentation(A : Precision := 30);
printf "mrep(i) = %o\nmrep(j) = %o\n", mrep(A.1), mrep(A.2);
print "DETPROBE DONE";
quit;
