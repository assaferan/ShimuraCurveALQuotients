SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/nelson.m";
for DD in [6,10,15] do
A := QuaternionAlgebra(DD); O := MaximalOrder(A);
bas := [A!b : b in ZBasis(O)];
T := TraceGram(bas, A);
printf "=== D=%o  det T = %o\n", DD, Determinant(T);
for d in Divisors(DD) do
  B := RdBasis(T, DD, d);
  ls := [ &+[ B[i][j]*bas[j] : j in [1..4] ] : i in [1..4] ];
  // (i) contains R?
  Bi := B^(-1);
  containsR := &and[ IsIntegral(x) : x in Eltseq(Bi) ];
  // (ii) inside R^*?  trd(ls[i]*bas[j]) in Z
  inDual := &and[ IsIntegral(Trace(A!(x*y))) : x in ls, y in bas ];
  // (iii) index
  idx := 1/Determinant(B);
  // (iv) nrd denominators
  S := {};
  for a,b,c,e in [-2..2] do Include(~S, Denominator(Norm(a*ls[1]+b*ls[2]+c*ls[3]+e*ls[4]))); end for;
  printf "  d=%-3o contains R: %o   inside R^*: %o   index %o (want %o)   nrd denominators %o\n",
         d, containsR, inDual, idx, d^2, Sort(SetToSequence(S));
end for;
end for;
quit;
