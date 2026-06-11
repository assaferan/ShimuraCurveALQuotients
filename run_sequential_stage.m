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
        // Verify FilterByTraceStar output before running HHProposition1
        VerifyHHTable2(curves);
        printf "VerifyHHTable2 passed\n";
        HHProposition1(~curves);

    when "GetQuotientsAndGenera_UpdateByGenus":
        // Expand star curves into all AL quotients, then classify by genus
        t1 := Realtime();
        curves := GetQuotientsAndGenera(curves);
        printf "GetQuotientsAndGenera took %o\n", Realtime() - t1;
        t0 := Realtime();
        UpdateByGenus(~curves);

    when "UpdateCurves1", "UpdateCurves2", "UpdateCurves3", "UpdateCurves4",
         "UpdateCurves5", "UpdateCurves6", "UpdateCurves7", "UpdateCurves8":
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
    when "HHProposition1":
        VerifyHHProposition1(curves);
        printf "VerifyHHProposition1 passed\n";
    when "GetQuotientsAndGenera_UpdateByGenus":
        VerifyFHTheorem3(curves);
        printf "VerifyFHTheorem3 passed\n";
    when "UpdateCurves5":
        VerifyFHTable3();
        printf "VerifyFHTable3 passed\n";
end case;

catch e
    WriteStderr(Sprintf("ERROR in sequential stage %o:\n", stage));
    WriteStderr(e);
    error e;  // re-raise so SetQuitOnError exits non-zero
end try;

Write(output_dat, Sprint(curves, "Magma") : Overwrite);
printf "%o: wrote %o curves to %o\n", stage, #curves, output_dat;
quit;
