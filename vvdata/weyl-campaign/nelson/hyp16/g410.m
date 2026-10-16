AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
for D in [410] do
  full := { Integers() | d : d in Divisors(D) | GCD(d, D div d) eq 1 };
  printf "D = %o, Hall divisors %o\n", D, Sort(SetToSequence(full));
  for W in [ {Integers()|1}, full ] do
    try g := GenusShimuraCurveQuotient(D,1,W); printf "  W = %o : genus %o\n", Sort(SetToSequence(W)), g;
    catch e printf "  W = %o : FAILS -- %o\n", Sort(SetToSequence(W)), e`Object; end try;
  end for;
end for;
printf "G410 DONE\n"; quit;
