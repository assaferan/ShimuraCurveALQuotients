// Rotger records V_291/W -- the star curve of X_0^291(1) -- as a Shimura curve quotient whose
// hyperelliptic involution is EXCEPTIONAL, i.e. not Atkin-Lehner.  If its genus is 3 or 4 it is
// the positive control this detector has never had: a curve in the same structural position as
// the sixteen (no Atkin-Lehner overgroup of genus 0) that really IS hyperelliptic.
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
for D in [291] do
  halls := [ m : m in Divisors(D) | GCD(m, D div m) eq 1 ];
  printf "X_0^%o(1): Hall divisors %o\n", D, halls;
  for W in Subsets(Set(halls)) do
    if #W eq 0 or not (1 in W) then continue; end if;
    ok := true;
    for a,b in W do if not ((a*b) div GCD(a,b)^2 in W) then ok := false; end if; end for;
    if not ok then continue; end if;
    printf "   W = %-22o genus %o\n", Sort(SetToSequence(W)), GenusShimuraCurveQuotient(D,1,W);
  end for;
end for;
printf "V291 DONE\n";
quit;
