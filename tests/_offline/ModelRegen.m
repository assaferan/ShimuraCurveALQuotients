// tests/_offline/ModelRegen.m
//
// DOES EACH COMMITTED MODEL STILL REGENERATE FROM CURRENT CODE?
//
// WHY THIS EXISTS. Nothing else in the suite would notice a committed model file drifting away
// from what the pipeline now produces:
//   * ModelChecks.m validates STORED models structurally (genus, Weil divisibility, trace-formula
//     point counts) -- it never runs the pipeline;
//   * GuoYangEquations.m compares STORED models to the literature -- it never runs the pipeline;
//   * the X0_D_N.m tests DO run the pipeline, but only for the ~25 bases that have such a file,
//     and each needs hand-written cover/AL data.
// So "we can no longer produce this model" was invisible. Measured 2026-09-05: THREE committed
// files -- models_22_5.m, models_22_3.m, models_15_2.m -- do not reproduce, all for the same
// reason (the unpinned-y2-scale guard, commit 1768517, postdates them). They are not WRONG --
// they pass ModelChecks and their Guo-Yang comparisons -- they are unreproducible.
//
// WHAT IT CHECKS, per base: re-run AllEquationsAboveCovers, aggregate exactly as genmodels.m
// does, and compare against the committed file:
//   [1] every populated cover key in the committed file is still produced, and
//   [2] the regenerated curve is ISOMORPHIC to the committed one.
// Isomorphism, not byte-equality: a model is only defined up to isomorphism, so a differently
// presented but isomorphic curve is not a regression. A MISSING cover is.
//
// ⚠ THIS IS A DRIFT TEST, NOT A VALIDATION. It says "current code still produces this", not
// "this is correct". Correctness comes from ModelChecks (independent trace-formula point counts)
// and GuoYangEquations (the literature). A base can pass here and be wrong, or fail here and be
// right -- 22_3/15_2/22_5 are exactly the latter.
//
// COST. This regenerates, so it is SLOW and lives in tests/_offline/ (invisible to run_tests.m's
// sweep and to the CI matrix, per the leading-underscore convention). Measured single-base times
// vary by three orders of magnitude: 14_3 37s, 51_1 240s, 22_5 680s, 10_19 ~3600s, 87_1 ~7000s.
// Run a subset with the MODELREGEN_BASES env var:
//     MODELREGEN_BASES=14_3,51_1 magma -b filename:=tests/_offline/ModelRegen.m run_tests.m
// Unset, it runs CHEAP_BASES below -- a spread that completes in ~20 min -- rather than all 81,
// which would take tens of hours.
// ⚠ It must be an ENV VAR, not a `name:=value` argument: run_tests.m executes test files via
// `eval`, and a command-line variable is NOT visible inside that context (`assigned bases` is
// false there, so the argument silently does nothing and the default list runs instead). Same
// reason INTSOL / M0PROGRESS / Y2TWIST are env-gated.
//
// NB: top-level statements, not a procedure -- run_tests.m uses `eval`, and an `eval` inside a
// procedure closing over an outer variable segfaults Magma 2.29 (same trap as ModelChecks.m).

// A cheap spread: both D parities, hyperelliptic and CRV entries, and the three known-stale bases
// so the failure they represent stays visible instead of being quietly dropped.
CHEAP_BASES := ["14_3", "51_1", "57_1", "6_11", "35_1", "38_1", "22_3", "15_2", "22_5"];

mr_sel := CHEAP_BASES;
mr_env := GetEnv("MODELREGEN_BASES");
if mr_env ne "" then
    mr_sel := Split(mr_env, ",");
end if;

printf "ModelRegen: regenerating %o base(s) and comparing to the committed models...\n", #mr_sel;

mr_ok := 0; mr_fail := 0; mr_bad := [];

for mr_b in mr_sel do
    mr_dn := Split(mr_b, "_");
    error if #mr_dn ne 2, Sprintf("ModelRegen: cannot parse base %o (want D_N)", mr_b);
    mr_D := StringToInteger(mr_dn[1]);
    mr_N := StringToInteger(mr_dn[2]);
    mr_f := Sprintf("data/models/models_%o_%o.m", mr_D, mr_N);
    mr_stored := eval (Read(mr_f) cat "\nreturn models;");

    mr_curves := GetHyperellipticCandidates();
    mr_Xstar := rep{X : X in mr_curves | X`D eq mr_D and X`N eq mr_N and IsStarCurve(X)};
    mr_covers := AllEquationsAboveCovers(mr_Xstar, mr_curves : Prec := 100);

    // aggregate by W-key exactly as genmodels.m does
    mr_fresh := AssociativeArray();
    for mr_lab in Keys(mr_covers) do
        mr_W := Sort([Integers()| w : w in mr_curves[mr_lab]`W]);
        if not IsDefined(mr_fresh, mr_W) then mr_fresh[mr_W] := [* *]; end if;
        for mr_base in Keys(mr_covers[mr_lab]) do
            Append(~mr_fresh[mr_W], mr_covers[mr_lab][mr_base]);
        end for;
    end for;

    mr_missing := []; mr_noniso := []; mr_skipped := 0; mr_cmp := 0;
    for mr_k in Keys(mr_stored) do
        if #mr_stored[mr_k] eq 0 then continue; end if;      // empty entries carry no claim
        mr_key := Sort([Integers()| w : w in mr_k]);
        mr_have, mr_list := IsDefined(mr_fresh, mr_key);
        if (not mr_have) or #mr_list eq 0 then
            Append(~mr_missing, mr_key); continue;
        end if;
        // Match as a MULTISET: each fresh cover can absorb at most one committed entry. Without
        // consuming matches, a key holding several curves passes when some were LOST, because two
        // committed entries both match the same surviving fresh one -- which is exactly how
        // 22_3's [1,66] (3 committed entries, 2 regenerated) went unreported on the first draft.
        mr_used := {};
        for mr_e in mr_stored[mr_k] do
            mr_g, mr_fo, mr_h := Explode(mr_e);
            if Type(mr_fo) eq MonStgElt then
                // CRV entry: skipped. Building the weighted-projective model needs the ambient
                // weights, which the file does not record -- see GuoYangEquations.m for how those
                // are pinned by hand.
                mr_skipped +:= 1;
                continue;
            end if;
            mr_Cst := HyperellipticCurve(mr_fo, mr_h);
            mr_found := false;
            for mr_i -> mr_c in mr_list do
                if mr_i in mr_used then continue; end if;
                if Type(mr_c) ne CrvHyp then continue; end if;
                if Genus(mr_c) ne Genus(mr_Cst) then continue; end if;
                if IsIsomorphic(mr_c, mr_Cst) then
                    mr_found := true; Include(~mr_used, mr_i); break;
                end if;
            end for;
            mr_cmp +:= 1;
            if not mr_found then Append(~mr_noniso, mr_key); end if;
        end for;
    end for;

    if IsEmpty(mr_missing) and IsEmpty(mr_noniso) then
        printf "  %o_%o: OK (%o cover(s) compared, %o CRV skipped)\n", mr_D, mr_N, mr_cmp, mr_skipped;
        mr_ok +:= 1;
    else
        printf "  %o_%o: DRIFTED -- missing %o, non-isomorphic %o (%o compared, %o CRV skipped)\n",
               mr_D, mr_N, mr_missing, mr_noniso, mr_cmp, mr_skipped;
        mr_fail +:= 1;
        Append(~mr_bad, Sprintf("%o_%o", mr_D, mr_N));
    end if;
end for;

printf "ModelRegen: %o base(s) reproduce, %o drifted%o\n", mr_ok, mr_fail,
       IsEmpty(mr_bad) select "" else Sprintf(" (%o)", mr_bad);

// KNOWN-DRIFTED, recorded rather than asserted away. Update this list -- do not delete the check.
//
//   22_3, 15_2, 22_5  unreproducible because the unpinned-y2-scale guard (1768517) POSTDATES them;
//                     the committed models are correct but the pipeline now withholds covers they
//                     contain. Y2TWIST=1 restores 22_5 and 15_2 fully, 22_3 13/14.
//   14_3              same as 39_2: produced with CMNONCOPRIME=1. Under default settings its
//                     covers are under-determined and W={1} comes out empty. Validated by the
//                     offline full-curve comparison, not by regeneration.
//   39_2              a DIFFERENT reason: it was produced with CMNONCOPRIME=1 and cannot be
//                     regenerated by default at all -- with the coprime CM filter on it sees 3
//                     points against demand 19 and dies with "Could not find enough points". Its
//                     drift here is EXPECTED and is documented in the model file header; the file
//                     is validated by GuoYangEquations.m (IsIsomorphic to Guo-Yang's published
//                     genus-7 curve) and by ModelChecks, not by regeneration.
// 14_43 added 2026-09-06: produced with INTSOL=1 (see data/models/PROVENANCE.md), so it does
// not regenerate under the default recipe. Keep this list and that table in sync.
MR_KNOWN_DRIFT := ["22_3", "15_2", "22_5", "39_2", "14_3", "14_43"];
mr_unexpected := [b : b in mr_bad | not (b in MR_KNOWN_DRIFT)];
error if not IsEmpty(mr_unexpected),
    Sprintf("ModelRegen: base(s) %o no longer reproduce and are NOT in the known-drift list",
            mr_unexpected);
