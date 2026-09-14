// Verify a freshly generated model set that is NOT yet in data/models/, AND negative-control the
// verification in the same breath.  Used 15 times on 2026-09-14; previously re-sed'd per base from
// the scratchpad, which CLAUDE.md warns against (/tmp is purged nightly and has eaten a driver).
//
//   magma -b D_s:=46 N_s:=3 MODEL:=/path/to/models_46_3.m vvdata/weyl-campaign/modelcheck.m < /dev/null
//   ... NEG:=1   also runs the negative control
//
// ⚠ WHY THE NEGATIVE CONTROL IS NOT OPTIONAL.  "A PASSING CHECK IS NOT EVIDENCE UNTIL YOU KNOW IT
// COULD HAVE FAILED."  VerifyModelSet skips entries it cannot handle, so a clean pass can mean
// "nothing was checked".  The control perturbs ONE genus>=1 entry by the non-square twist -1 --
// invariant-preserving in genus and point counts at first glance, but detectable -- and requires
// the check to go red.  Typical: 46_3 116/0 -> 116/3, 35_2 146/0 -> 146/7.
SetQuitOnError(true); SetColumns(0);
AttachSpec("ShimuraQuotients.spec");
D := StringToInteger(D_s); N := StringToInteger(N_s);
// ⚠ Read+eval, NOT `load`: Magma's load needs a LITERAL path, so a variable fails with
// "Identifier 'models' has not been declared".  This is how tests/GonzalezRotger.m does it.
models := eval (Read(MODEL) cat "\nreturn models;");
nchk, nfail := VerifyModelSet(models, D, N);
printf "\nVERIFY %o_%o : checks %o, failures %o\n", D, N, nchk, nfail;
if assigned NEG then
    done := false;
    for k in Sort(SetToSequence(Keys(models))) do
        if done then break; end if;
        for i in [1..#models[k]] do
            e := models[k][i];
            if Type(e[2]) eq MonStgElt or e[1] lt 1 then continue; end if;
            printf "PERTURBING key %o (genus %o) by the non-square twist -1\n", k, e[1];
            models[k][i] := <e[1], -1*e[2], e[3]>;
            done := true; break;
        end for;
    end for;
    error if not done, "no genus>=1 entry to perturb -- the control cannot run, so the pass above is UNVALIDATED";
    nchk2, nfail2 := VerifyModelSet(models, D, N : Verbose := false);
    printf "NEGCTL %o_%o : checks %o, failures %o  (MUST be > 0)\n", D, N, nchk2, nfail2;
    error if nfail2 eq 0,
        "NEGATIVE CONTROL FAILED: twisting a genus>=1 entry did not make VerifyModelSet fail, so "
        * "the clean pass above verifies NOTHING for this base";
end if;
quit;
