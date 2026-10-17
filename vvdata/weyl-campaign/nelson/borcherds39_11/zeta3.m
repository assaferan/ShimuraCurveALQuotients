SetColumns(0);
AttachSpec("ShimuraQuotients.spec");
P<x> := PolynomialRing(Rationals());
curves := GetHyperellipticCandidates();
assert exists(X3){Y : Y in curves | Y`CurveID eq 9679};
assert exists(Xs){Y : Y in curves | Y`CurveID eq 9681};
C3 := HyperellipticCurve(-3*x^7-22*x^6-55*x^5-69*x^4-55*x^3-22*x^2-3*x, x^4+x^3+x^2+x+1);
Cst := IntegralModel(HyperellipticCurve(-25447851/236421376*x^4 - 573417/476656*x^3 + 18083/7688*x^2 + 13*x + 13));
// counting #C(F_p^r) for r = 1..g determines the numerator of the zeta function completely.
procedure zchk(C, X, g, name, ~nchk, ~nfail, prbd)
  for p in [q : q in PrimesUpTo(prbd) | 429 mod q ne 0] do
    ok := true;
    try
      Cp := ChangeRing(C, GF(p));
      if not IsNonsingular(Cp) then ok := false; end if;
    catch e ok := false; end try;
    if not ok then printf "  %o p=%-3o model has bad reduction, skipped\n", name, p; continue; end if;
    got := []; want := [];
    for r in [1..g] do
      Append(~got,  #Points(ChangeRing(C, GF(p^r))));
      Append(~want, ComputePointsViaTrace(X, p, r));
    end for;
    nchk +:= g;
    if got eq want then
      printf "  %o p=%-3o counts over F_p..F_p^%o AGREE: %o\n", name, p, g, got;
    else
      nfail +:= g;
      printf "  %o p=%-3o MISMATCH model %o vs trace %o\n", name, p, got, want;
    end if;
  end for;
end procedure;
nchk := 0; nfail := 0;
printf "=== star X_0^39(11)*, genus 1 ===\n";
zchk(Cst, Xs, 1, "star", ~nchk, ~nfail, 100);
printf "=== X_0^39(11)/[1,11,39,429], genus 3 ===\n";
zchk(C3, X3, 3, "g3  ", ~nchk, ~nfail, 80);
printf "\nTOTAL zeta-determining counts %o, mismatches %o\n", nchk, nfail;
printf "ZETA2 DONE\n";
exit;
