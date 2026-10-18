// Which still-open quotients are bielliptic via an Atkin-Lehner involution, by genus; and the
// full-curve examples X_0^D(1) that are both hyperelliptic and bielliptic (all must have g <= 3,
// by Castelnuovo-Severi).
AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
cs := eval Read("data/all_curves_progress.dat");
byDN := AssociativeArray();
for X in cs do k := <X`D, X`N>; if not IsDefined(byDN,k) then byDN[k] := []; end if; Append(~byDN[k], X); end for;
open := [X : X in cs | not assigned X`TestInWhichProved or X`TestInWhichProved eq ""];
printf "open: %o\n", #open;
for X in open do
  overs := [Y : Y in byDN[<X`D,X`N>] | #Y`W eq 2*#X`W and X`W subset Y`W];
  ell := [Y`W : Y in overs | Y`g eq 1];
  printf "%o_%o W=%o g=%o  elliptic AL quotients: %o\n", X`D, X`N, X`W, X`g, #ell;
end for;
printf "\nfull curves X_0^D(1) hyperelliptic with an elliptic AL quotient:\n";
for X in cs do
  if X`N ne 1 or #X`W ne 1 or not (assigned X`IsHyp and X`IsHyp) then continue; end if;
  overs := [Y : Y in byDN[<X`D,1>] | #Y`W eq 2 and Y`g eq 1];
  if #overs gt 0 then printf "  D=%o g=%o elliptic quotients by %o\n", X`D, X`g, [Y`W : Y in overs]; end if;
end for;
