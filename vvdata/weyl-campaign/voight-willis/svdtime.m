// How long does a decomposition take?  The cost model in the README was an estimate and is
// about 700 times too fast; this is the measurement that replaces it.  Cost goes as n^2.8 for a
// real m x n matrix with m = 2n, and barely depends on the working precision.
//   magma -b svdtime.m < /dev/null
SetColumns(0);
for pr in [40, 70, 120] do
  RR := RealField(pr);
  for n in [80, 160, 320] do
    m := 2*n;
    M := Matrix(RR, m, n, [RR| (i*7 mod 101)/101 + (i mod 13)/17 : i in [1..m*n]]);
    t := Cputime();
    S := SingularValueDecomposition(M);
    printf "prec %-4o %o x %-5o : %o s\n", pr, m, n, RealField(4)!Cputime(t);
  end for;
end for;
quit;
