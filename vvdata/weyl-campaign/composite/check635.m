// The X_0^6(35) models from lovelace against the Eichler-Selberg trace formula: for every stored
// y^2 + h y = f of genus >= 1, #(X/W)(F_p) and #(X/W)(F_{p^2}) at the primes 11, 13, 17, 19 (all
// coprime to 210) must equal ComputePointsViaTrace.  The trace formula knows nothing about CM values,
// so this is the first outside test of the composite-level m = 0 multipliers (prop:composite).
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
P<x> := PolynomialRing(Rationals());
models := eval (Read(MODELS) cat "\nreturn models;");
curves := GetHyperellipticCandidates();
D := 6; N := 35;
nok := 0; nbad := 0; nskip := 0;
for k in Keys(models) do
    if #models[k] eq 0 then continue; end if;
    e := models[k][1];
    if Type(e[2]) eq MonStgElt then nskip +:= 1; continue; end if;    // CRV entries: checked through their factors
    g := e[1]; f := e[2]; h := #e ge 3 select e[3] else P!0;
    if g eq 0 then nskip +:= 1; continue; end if;                     // a conic twist is invisible to point counts
    X := rep{Y : Y in curves | Y`D eq D and Y`N eq N and Y`W eq Set(k)};
    assert X`g eq g;
    C := HyperellipticCurve(f, h);
    for p in [11, 13, 17, 19] do
        Cp := ChangeRing(C, GF(p));
        if Genus(Cp) ne g then continue; end if;
        cnt := [#Points(BaseChange(Cp, GF(p^d))) : d in [1..2]];
        exp := [ComputePointsViaTrace(X, p, d) : d in [1..2]];
        if cnt eq exp then nok +:= 1; else nbad +:= 1; printf "MISMATCH W=%o genus %o p=%o: counts %o, trace formula %o\n", k, g, p, cnt, exp; end if;
    end for;
end for;
printf "X_0^6(35): %o (prime, cover) point counts agree with the trace formula, %o disagree; %o entries skipped (genus 0 or CRV)\n", nok, nbad, nskip;
