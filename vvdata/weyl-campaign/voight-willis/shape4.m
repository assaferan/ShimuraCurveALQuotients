RR := RealField(30);
M := Matrix(RR, 7, 3, [RR | Random(1,100)/7 : i in [1..21]]);
Sig, L, Rt := SingularValueDecomposition(M);
sv := [RealField(10)!Sig[i][i] : i in [1..3]];
printf "singular values (diag of 1st return, descending): %o\n", sv;
for j in [1..3] do
    vcol := Matrix(RR, 3, 1, [Rt[i][j] : i in [1..3]]);
    vrow := Matrix(RR, 3, 1, [Rt[j][i] : i in [1..3]]);
    nc := Sqrt(&+[x^2 : x in Eltseq(M*vcol)]);
    nr := Sqrt(&+[x^2 : x in Eltseq(M*vrow)]);
    printf "  j=%o sv=%o  ||M*(col j)||=%o  ||M*(row j)||=%o\n", j, RealField(8)!sv[j], RealField(8)!nc, RealField(8)!nr;
end for;
quit;
