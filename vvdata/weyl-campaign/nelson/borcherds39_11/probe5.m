SetColumns(0);
AttachSpec("ShimuraQuotients.spec");
curves := GetHyperellipticCandidates();
stars := [X : X in curves | IsStarCurve(X) and X`D ne 1 and X`g eq 1 and IsSquarefree(X`N)];
cands := [];
for X in stars do
  D0 := (X`D*X`N) div 2^Valuation(X`D*X`N,2); M := 4*D0;
  if #Divisors(M) le 16 and M le 500 then
    nsub := #[i : i in X`CoveredBy | assigned curves[i]`IsSubhyp and curves[i]`IsSubhyp];
    Append(~cands, <M, #Divisors(M), X`D, X`N, nsub, #X`CoveredBy>);
  end if;
end for;
Sort(~cands);
printf "genus-1 star bases with small M (M, #div, D, N, #subhyp covers, #covers):\n";
for c in cands do printf "  %o\n", c; end for;
printf "DONE\n";
exit;
