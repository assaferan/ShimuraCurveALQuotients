// vvdata/weyl-campaign/covercount.m -- hand-run probe. Campaign branch, per CLAUDE.md's rule that
// scratch scripts live here and not in /tmp.
//
//     magma -b Dd:=10 Nn:=3 vvdata/weyl-campaign/covercount.m < /dev/null
//     PTSCOPRIME=1 magma -b Dd:=10 Nn:=3 vvdata/weyl-campaign/covercount.m < /dev/null
//
// ⇒ DOES THE PIPELINE STILL PRODUCE EVERY COVER THE COMMITTED FILE RECORDS?
// Compares, per W key, the number of covers a fresh AllEquationsAboveCovers run emits against the
// number of entries in data/models/models_D_N.m.  A key producing FEWER is a LOST COVER.
//
// WHAT IT FOUND (2026-09-23/24), which is why it is kept: 0ca6e37 dropped the coprime-to-level
// filter on the divisor-support CM pool and thereby changed the ORDER in which the (infty,P,Q)
// sweep meets anchors.  The sweep takes the FIRST workable triple, so the anchor moved, the
// hauptmodul was re-normalised, and covers went missing at three bases:
//
//     base   default        PTSCOPRIME=1    (PTSCOPRIME=1 restores the pre-0ca6e37 pool)
//     10_3   23 of 26       26 of 26
//     6_13   24 of 27       27 of 27
//     26_5   11 of 12       no shortfall
//     6_5 6_7 22_7 6_23 6_71   clean either way
//
// Fixed by sorting the pool coprime-to-N first (PR #41).  Keep this probe: it is the cheapest way
// to re-check the property if the CM pool, the anchor sweep, or the hauptmodul normalisation is
// ever touched again.
//
// ⚠ SCOPING IT CHEAPLY.  Do NOT sweep all 76 N>1 bases by hand.  Only 26 have any key with >= 2
// entries (a free grep over data/models/), and the X0_ helper's SECOND pass tests this exact
// shortfall on every committed key with model_covers on -- so every base whose X0_ test passes is
// already proven clean.  That left four bases needing a manual run.
AttachSpec("ShimuraQuotients.spec");
D := StringToInteger(Dd); N := StringToInteger(Nn);
models := eval (Read(Sprintf("data/models/models_%o_%o.m", D, N)) cat "\nreturn models;");
curves := GetHyperellipticCandidates();
assert exists(Xstar){X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)};
covers, ws := AllEquationsAboveCovers(Xstar, curves);
printf "%o_%o  PTSCOPRIME=%o\n", D, N,
       GetEnv("PTSCOPRIME") eq "" select "unset (current default)" else GetEnv("PTSCOPRIME");
tot_c := 0; tot_m := 0; short := [];
for label in Keys(covers) do
    W := Sort([Integers()| w : w in curves[label]`W]);
    if not IsDefined(models, W) then continue; end if;
    nc := #Keys(covers[label]); nm := #models[W];
    tot_c +:= nc; tot_m +:= nm;
    if nc ne nm then
        printf "   W=%-14o committed %o, produced %o   %o\n", Sprint(W), nm, nc,
               nc lt nm select "<== SHORT" else "(extra)";
    end if;
    if nc lt nm then Append(~short, W); end if;
end for;
printf "   TOTALS over keys present in the file: committed %o, produced %o%o\n",
       tot_m, tot_c, IsEmpty(short) select "   NO SHORTFALL" else Sprintf("   SHORT AT %o", short);
exit;
