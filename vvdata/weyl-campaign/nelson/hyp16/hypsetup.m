// Fix the Eichler order once, for every chunk of a split evaluation.
//
// Magma's QuaternionOrder(A, N) returns a DIFFERENT Eichler order of level N from one process
// to the next (measured: three runs at D = 55, N = 7 gave two distinct Z-bases; the algebra's
// standard form and FuchsianMatrixRepresentation were identical each time).  Nelson's values
// at a point z depend on the order through the lattices R^(d) and through which point of the
// curve z names, so chunks evaluated in separate processes against separately chosen orders do
// not lie on one canonical conic (each chunk alone does; see conicdiag.m).  This writes the
// Z-basis of one order, in the coordinates 1, i, j, k of QuaternionAlgebra(D), for every chunk
// to rebuild.
//
//   magma -b DISC:=55 LEV:=7 OUT:=setup_55_7.m hypsetup.m < /dev/null
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
D := StringToInteger(DISC); N := StringToInteger(LEV);
A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
assert Discriminant(O) eq D*N;
a, b := StandardForm(A);
F := Open(OUT, "w");
fprintf F, "// Eichler order of level %o in QuaternionAlgebra(%o) (standard form %o, %o), Z-basis in 1,i,j,k\n", N, D, a, b;
fprintf F, "SA := %o; SB := %o;\nOBAS := %o;\n", a, b, [ Eltseq(A ! x) : x in ZBasis(O) ];
delete F;
printf "wrote %o: %o\n", OUT, [ Eltseq(A ! x) : x in ZBasis(O) ];
print "HYPSETUP DONE";
quit;
