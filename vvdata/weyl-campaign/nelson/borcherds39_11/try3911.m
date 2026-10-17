SetColumns(0);
AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 3);
curves := GetHyperellipticCandidates();
assert exists(Xs){X : X in curves | X`D eq 39 and X`N eq 11 and IsStarCurve(X)};
printf "X_0^39(11)* genus = %o\n", Xs`g;
t0 := Realtime();
try
  eqns, ws := AllEquationsAboveCovers(Xs, curves : Targets := {{1,11,39,429}});
  printf "\nRETURNED after %o s. keys = %o\n", Realtime()-t0, Keys(eqns);
  for k in Keys(eqns) do printf "  %o : %o\n", k, eqns[k]; end for;
catch e
  printf "\nFAILED after %o s:\n%o\n", Realtime()-t0, e;
end try;
printf "TRY3911 DONE\n";
exit;
