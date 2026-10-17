SetColumns(0);
AttachSpec("ShimuraQuotients.spec");
curves := GetHyperellipticCandidates();
assert exists(T){Y : Y in curves | Y`CurveID eq 9679};   // X/[1,11,39,429], genus 3
assert exists(Xs){Y : Y in curves | Y`CurveID eq 9681};  // the star
for X in [T, Xs] do
  printf "=== W = %o (genus %o) ===\n", Sort(SetToSequence(X`W)), X`g;
  for d in [-52, -132] do
    printf "  d = %-5o degree of field of definition = %o\n", d, DegreeOfFieldOfDefinitionOfCMPoint(X, d);
  end for;
  D2 := Degree2Points(X);
  for d in [-52,-132] do
    if IsDefined(D2, d) then printf "  d = %-5o degree-2 point data: %o\n", d, D2[d]; end if;
  end for;
end for;
printf "CMF DONE\n";
exit;
