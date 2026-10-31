// The CUSPIDAL space S_{3/2}(Gamma_0(M)) (theta^3 multiplier) to precision B, cached as
// basis_<M>_<B>.txt (the format of solve.m), a matrix over Q whose rows are the basis forms and whose column n
// holds the coefficient of q^n, n = 1..B.  Used by the sharp dual test (wtest.m), which only needs
// enough precision to pin the forms supported on the CM indices.
//   magma -b MM:=6780 BB:=1100 basisS.m
SetColumns(0);
M := StringToInteger(MM); B := StringToInteger(BB);
t0 := Realtime();
Mf := HalfIntegralWeightForms(Gamma0(M), 3/2);
printf "M = %o: dim M_{3/2} = %o, dim S_{3/2} = %o (%o s)\n", M, Dimension(Mf), Dimension(CuspidalSubspace(Mf)), Realtime() - t0;
S := CuspidalSubspace(Mf);
Bq := Basis(S, B + 2);
A0 := Matrix(Rationals(), #Bq, B, [ Coefficient(g, n) : n in [1..B], g in Bq ]);
cache := Sprintf("basis_%o_%o.txt", M, B);
Write(cache, Sprint(A0, "Magma") : Overwrite := true);
printf "basis to precision %o written to %o in %o s\n", B + 2, cache, Realtime() - t0;
quit;
