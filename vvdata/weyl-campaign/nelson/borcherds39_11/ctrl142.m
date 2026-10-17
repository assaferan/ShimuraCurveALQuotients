SetColumns(0);
AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 1);
curves := GetHyperellipticCandidates();
assert exists(Xs){X : X in curves | X`D eq 142 and X`N eq 1 and IsStarCurve(X)};
printf "X_0^142(1)* genus = %o ; covers = %o\n", Xs`g, [<curves[i]`CurveID, Sort(SetToSequence(curves[i]`W)), curves[i]`g> : i in Xs`CoveredBy];
t0 := Realtime();
try
  eqns, ws := AllEquationsAboveCovers(Xs, curves);
  printf "\nRETURNED after %o s. keys = %o\n", Realtime()-t0, Keys(eqns);
  for k in Keys(eqns) do printf "  %o : %o\n", k, eqns[k]; end for;
catch e
  printf "\nFAILED after %o s:\n%o\n", Realtime()-t0, e;
end try;
printf "CTRL142 DONE\n";
exit;
