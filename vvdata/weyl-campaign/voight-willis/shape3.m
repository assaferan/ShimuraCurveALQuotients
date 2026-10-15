RR := RealField(30);
M := Matrix(RR, 7, 3, [RR | Random(1,100)/7 : i in [1..21]]);
U, S, V := SingularValueDecomposition(M);
tru := Sort([RealField(8)!Sqrt(e[1]) : e in Eigenvalues(Transpose(M)*M)]);
printf "true singular values: %o\n", tru;
printf "U %ox%o  S %ox%o  V %ox%o\n", Nrows(U),Ncols(U),Nrows(S),Ncols(S),Nrows(V),Ncols(V);
function mx(X) return RealField(5)!Max([Abs(x) : x in Eltseq(X)]); end function;
printf "U orthonormal cols? |U^T U - I| = %o\n", mx(Transpose(U)*U - IdentityMatrix(RR,Ncols(U)));
printf "S orthogonal?       |S^T S - I| = %o\n", mx(Transpose(S)*S - IdentityMatrix(RR,Nrows(S)));
printf "V orthogonal?       |V^T V - I| = %o\n", mx(Transpose(V)*V - IdentityMatrix(RR,Nrows(V)));
printf "U diag entries: %o\n", [RealField(6)!U[i][i] : i in [1..3]];
// try reconstructions
try printf "|M - S*U*Transpose(V)| = %o\n", mx(M - S*U*Transpose(V)); catch e printf "S*U*V^T bad shape\n"; end try;
try printf "|M - S*U*V| = %o\n", mx(M - S*U*V); catch e printf "S*U*V bad shape\n"; end try;
quit;
