RR := RealField(30);
M := Matrix(RR, 7, 3, [RR | Random(1,100)/7 : i in [1..21]]);
U, S, V := SingularValueDecomposition(M);
printf "M %o x %o -> U %o x %o, S %o x %o, V %o x %o\n", Nrows(M), Ncols(M), Nrows(U), Ncols(U), Nrows(S), Ncols(S), Nrows(V), Ncols(V);
printf "check |M - U*S*Transpose(V)| = %o\n", RealField(5)!Max([Abs(x) : x in Eltseq(M - U*S*Transpose(V))]);
printf "check |M - U*S*V|            = %o\n", RealField(5)!Max([Abs(x) : x in Eltseq(M - U*S*V)]);
printf "diag S: %o\n", [RealField(5)!S[i][i] : i in [1..Min(Nrows(S),Ncols(S))]];
quit;
