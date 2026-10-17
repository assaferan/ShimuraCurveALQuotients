SetColumns(0);
AttachSpec("ShimuraQuotients.spec");
D := 39; N := 11;
D0 := (D*N) div 2^Valuation(D*N,2);
M := 4*D0;
printf "D=%o N=%o  D0=%o  M=%o  #Divisors(M)=%o\n", D,N,D0,M,#Divisors(M);
printf "  Divisors(M) = %o\n", Divisors(M);
printf "  Genus(Gamma0(M)) = %o\n", Genus(Gamma0(M));
// comparison: a known-good base
for pr in [<39,2>,<15,2>,<6,11>,<39,1>] do
  d := pr[1]; n := pr[2];
  d0 := (d*n) div 2^Valuation(d*n,2); m := 4*d0;
  printf "  base %o_%o : M=%o  #div=%o\n", d,n,m,#Divisors(m);
end for;
printf "DONE PART 1\n";
curves := GetHyperellipticCandidates();
printf "#curves = %o\n", #curves;
sel := [X : X in curves | X`D eq 39 and X`N eq 11];
printf "#curves at 39_11 = %o\n", #sel;
for X in sel do
  printf "  W=%o genus=%o id=%o star=%o\n", Sort(SetToSequence(X`W)), X`g, X`CurveID, IsStarCurve(X);
end for;
printf "DONE\n";
exit;
