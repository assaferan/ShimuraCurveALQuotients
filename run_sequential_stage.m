// Runs a single sequential pipeline stage and writes the result to output_dat.
//
// Usage:
//   magma stage:=StageName [input_dat:=...] output_dat:=... run_sequential_stage.m
//
// input_dat is not required for FindPairs.
// The combined stage GetQuotientsAndGenera_UpdateByGenus runs both steps and
// saves under the name UpdateByGenus (matching the sequential pipeline filenames).

AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 0);
SetQuitOnError(true);

try

// Read input unless this is the first stage
if stage ne "FindPairs" then
    if not assigned input_dat then
        error "input_dat must be provided for stage " cat stage;
    end if;
    curves := eval Read(input_dat);
end if;

t0 := Realtime();

case stage:
    when "FindPairs":
        r := GetLargestPrimeIndex();
        assert r eq 7;
        curves := FindPairs(r);
        assert #curves eq 2342;

    when "UpdateGenera":
        UpdateGenera(~curves);

    when "UpdateByGenusStar":
        UpdateByGenus(~curves);

    when "HHProposition1":
        // The [HH] checks are about HH's own input, the FilterByTraceStar output, which sits next to
        // input_dat.  The input itself is FilterByTwistedTraceStar's output, and the twisted trace
        // decides further D = 1 star curves of HH Table 2 (e.g. X_0^*(396)), so it would not match.
        slash := [i : i in [1..#input_dat] | input_dat[i] eq "/"];
        dir := #slash eq 0 select "." else input_dat[1..slash[#slash]-1];
        hh := eval Read(dir cat "/curves_after_FilterByTraceStar.dat");
        VerifyHHTable2(hh);
        printf "VerifyHHTable2 passed\n";
        HHProposition1(~hh);
        VerifyHHProposition1(hh);
        printf "VerifyHHProposition1 passed\n";
        HHProposition1(~curves);

    when "SpecialFiberIsomorphismStar":
        SpecialFiberIsomorphism(~curves);

    when "GetQuotientsAndGenera_UpdateByGenus":
        // Expand star curves into all AL quotients, then classify by genus
        t1 := Realtime();
        curves := GetQuotientsAndGenera(curves);
        printf "GetQuotientsAndGenera took %o\n", Realtime() - t1;
        t0 := Realtime();
        UpdateByGenus(~curves);

    when "UpdateCurves5":
        VerifyFHTable3(curves);
        printf "VerifyFHTable3 passed\n";
        UpdateCurves(~curves);

    when "UpdateCurves1", "UpdateCurves2", "UpdateCurves3", "UpdateCurves4",
         "UpdateCurves6", "UpdateCurves7", "UpdateCurves8",
         "UpdateCurvesAfterAutomorphismGroup", "UpdateCurvesAfterTwistedTrace",
         "UpdateCurvesAfterTwistedWeilPolynomial":
        UpdateCurves(~curves);

    when "Genus3CoversGenus2":
        Genus3CoversGenus2(~curves);

    else
        error Sprintf("Unknown sequential stage: %o", stage);
end case;

printf "%o took %o s\n", stage, Realtime() - t0;

// Post-stage verifications (assertions that catch regressions)
case stage:
    when "UpdateGenera":
        VerifyHHTable1(curves);
        printf "VerifyHHTable1 passed\n";
    when "GetQuotientsAndGenera_UpdateByGenus":
        VerifyFHTheorem3(curves);
        printf "VerifyFHTheorem3 passed\n";
end case;

catch e
    WriteStderr(Sprintf("ERROR in sequential stage %o:\n", stage));
    WriteStderr(e);
    error e;  // re-raise so SetQuitOnError exits non-zero
end try;

Write(output_dat, Sprint(curves, "Magma") : Overwrite);
printf "%o: wrote %o curves to %o\n", stage, #curves, output_dat;
quit;
