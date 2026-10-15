AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
// genus-3 quotients of X_0^D(1) with NO genus-0 quotient one step above: for Shimura curves the
// automorphism group is the Atkin-Lehner group in this range, so such a quotient has no degree-2
// map to a rational curve and is not hyperelliptic.
for D in [1..260] do
  if not IsSquarefree(D) or D eq 1 or IsOdd(#Factorisation(D)) then continue; end if;
  ds := Set(Divisors(D));
  subs := [];
  for S in Subsets(ds) do
    if not (1 in S) or #S eq 0 then continue; end if;
    ok := true;
    for a in S do for b in S do
      if not ((a*b) div GCD(a,b)^2) in S then ok := false; end if;
    end for; end for;
    if ok then Append(~subs, S); end if;
  end for;
  for W in subs do
    if GenusShimuraCurveQuotient(D, 1, W) ne 3 then continue; end if;
    bad := false;
    for W2 in subs do
      if W subset W2 and #W2 eq 2*#W then
        if GenusShimuraCurveQuotient(D, 1, W2) eq 0 then bad := true; end if;
      end if;
    end for;
    if not bad then
      gs := [ <Sort(SetToSequence(W2)), GenusShimuraCurveQuotient(D,1,W2)>
              : W2 in subs | W subset W2 and #W2 eq 2*#W ];
      printf "D=%-5o W=%-22o genus 3, quotients above: %o\n", D, Sort(SetToSequence(W)), gs;
    end if;
  end for;
end for;
printf "SCAN4 DONE\n";
quit;
