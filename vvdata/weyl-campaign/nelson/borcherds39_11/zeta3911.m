SetColumns(0);
AttachSpec("ShimuraQuotients.spec");
P<x> := PolynomialRing(Rationals());
curves := GetHyperellipticCandidates();
assert exists(X){Y : Y in curves | Y`CurveID eq 9679};
assert exists(Xs){Y : Y in curves | Y`CurveID eq 9681};
C3 := HyperellipticCurve(-3*x^7-22*x^6-55*x^5-69*x^4-55*x^3-22*x^2-3*x, x^4+x^3+x^2+x+1);
Cst := HyperellipticCurve(-25447851/236421376*x^4 - 573417/476656*x^3 + 18083/7688*x^2 + 13*x + 13);
Cst := IntegralModel(Cst);

// Full L-polynomial comparison: the trace formula fixes #X(F_p^r) for r=1..g, which
// determines the whole numerator of the zeta function.  Strictly stronger than counting F_p.
procedure zetacheck(C, X, g, name, ~nchk, ~nfail)
  RT<TT> := PolynomialRing(Rationals());
  for p in [p : p in PrimesUpTo(60) | 429 mod p ne 0] do
    ok := true;
    try
      Cp := ChangeRing(C, GF(p));
      if not IsNonsingular(Cp) then ok := false; end if;
    catch e ok := false; end try;
    if not ok then printf "  %o p=%o : bad reduction of the model, skipped\n", name, p; continue; end if;
    Cp := ChangeRing(C, GF(p));
    Lmodel := RT!Numerator(ZetaFunction(Cp));
    // build L from the trace-formula counts
    cnt := [ComputePointsViaTrace(X, p, r) : r in [1..g]];
    // power sums s_r = p^r + 1 - #X(F_p^r) = sum of alpha_i^r  (alpha the 2g inverse roots)
    s := [p^r + 1 - cnt[r] : r in [1..g]];
    // Newton's identities -> coefficients e_r of prod (1 - alpha_i T) for r = 1..g
    e := [Rationals()|];
    for r in [1..g] do
      acc := s[r];
      for j in [1..r-1] do acc +:= (-1)^j * e[j] * s[r-j]; end for;
      Append(~e, (-1)^(r-1)*acc/r);
    end for;
    // L(T) = sum_{r=0}^{2g} a_r T^r with a_r = (-1)^r e_r, and functional equation a_{2g-r} = p^{g-r} a_r
    a := [Rationals()|1] cat [(-1)^r*e[r] : r in [1..g]];
    for r in [0..g-1] do Append(~a, p^(g-r)*a[r+1]); end for;
    Ltrace := &+[a[i+1]*TT^i : i in [0..2*g]];
    nchk +:= 1;
    if Ltrace ne Lmodel then
      nfail +:= 1;
      printf "  %o p=%o : MISMATCH\n    model %o\n    trace %o\n", name, p, Lmodel, Ltrace;
    else
      printf "  %o p=%o : L-polynomial AGREES (%o coefficients)\n", name, p, 2*g+1;
    end if;
  end for;
end procedure;

nchk := 0; nfail := 0;
printf "=== star X_0^39(11)* (genus 1) ===\n";
zetacheck(Cst, Xs, 1, "star", ~nchk, ~nfail);
printf "=== X_0^39(11)/[1,11,39,429] (genus 3) ===\n";
zetacheck(C3, X, 3, "g3", ~nchk, ~nfail);
printf "\nTOTAL full-zeta checks %o, failures %o\n", nchk, nfail;
printf "ZETA DONE\n";
exit;
