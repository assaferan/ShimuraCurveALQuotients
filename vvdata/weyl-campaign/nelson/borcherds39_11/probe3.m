SetColumns(0);
AttachSpec("ShimuraQuotients.spec");
curves := GetHyperellipticCandidates();
assert exists(Xs){X : X in curves | X`D eq 39 and X`N eq 11 and IsStarCurve(X)};
printf "Xstar 39_11: genus = %o, id = %o\n", Xs`g, Xs`CurveID;
printf "  CoveredBy = %o\n", Xs`CoveredBy;
// independent genus check via the standard formula components
printf "  NumberOfEllipticPoints q=2,3,4,6: %o\n", [NumberOfEllipticPoints(Xs,q) : q in [2,3,4,6]];
rams := RamficationPointsOfCovers(Xs, curves);
for k in Sort(SetToSequence(Keys(rams))) do
  printf "  cover id=%o W=%o g=%o  ram=%o\n", k, Sort(SetToSequence(curves[k]`W)), curves[k]`g, rams[k];
end for;
printf "=== CM points ===\n";
pts, rest := RationalandQuadraticCMPoints(Xs : bd := 2, coprime_to_level := false);
printf "  rational/quadratic pts (bd=2, no coprime filter): #=%o\n", #pts;
printf "  %o\n", pts;
ptsc, _ := RationalandQuadraticCMPoints(Xs : bd := 2, coprime_to_level := true);
printf "  coprime-filtered: #=%o  %o\n", #ptsc, ptsc;
printf "DONE\n";
exit;
