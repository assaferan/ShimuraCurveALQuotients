RR := RealField(30);
M := Matrix(RR, 7, 3, [RR | Random(1,100)/7 : i in [1..21]]);
U, S, V := SingularValueDecomposition(M);
printf "M 7x3 -> U %ox%o S %ox%o V %ox%o\n", Nrows(U),Ncols(U),Nrows(S),Ncols(S),Nrows(V),Ncols(V);
printf "FULL diagonal of S (7 entries, only 3 are genuine): %o\n", [RealField(6)!S[i][i] : i in [1..7]];
printf "genuine singular values via Ncols(M)=3:            %o\n", [RealField(6)!S[i][i] : i in [1..3]];
// independent check: sqrt of eigenvalues of M^T M
ev := Eigenvalues(Transpose(M)*M);
printf "sqrt eigenvalues of M^T M:                         %o\n", Sort([RealField(6)!Sqrt(e[1]) : e in ev]);
quit;
