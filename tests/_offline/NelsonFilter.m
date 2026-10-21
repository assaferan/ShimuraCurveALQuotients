// The Noether-count stage, against the banked readings it is meant to reproduce.
//
// ⚠ THIS FILE IS OFFLINE, and lives here for the reason the other files in this directory do: it
// costs 2150 s (26_3 alone is 312 s), far past what the default sweep should carry, and files
// under tests/_offline/ are invisible to both `ls tests/*.m` in run_tests.m and the CI matrix.
// Run it deliberately with
//     magma -b filename:=tests/_offline/NelsonFilter.m run_tests.m < /dev/null
// The cheap half of the job -- that the stage is inert unless asked and marks no curve when it
// does run -- is in tests/NelsonFilterContract.m and runs in the suite.
//
// SOURCE OF EVERY EXPECTED VALUE: vvdata/weyl-campaign/nelson/stage3_verdicts.txt on the campaign
// branch, "THIRD BATCH" table, which transcribes the logs of the stage3o.m runs. The three rows
// used here are the cheapest in it -- cutoff 40, 30 digits, D*N of 78, 154 and 238 -- chosen so
// the test runs in the suite rather than overnight:
//     X_0^26(3)/{1,39}   genus 3   1 relation    HYPERELLIPTIC
//     X_0^34(7)/{1,14}   genus 3   1 relation    HYPERELLIPTIC
//     X_0^14(11)/{1,11}  genus 3   0 relations   NOT hyperelliptic
// The third is the one that matters most: a precision floor puts a nonzero number where a
// vanishing one belongs, so it is the NEGATIVE reading that a floor can fake, and this is the
// case whose verdict the code refuses to give without a confirming second precision.
//
// WHAT THIS TEST IS NOT. It does not check that these curves ARE or are NOT hyperelliptic. The
// count is numerical and no verdict here is a proof; what is checked is that the stage reproduces
// the readings that were banked, so that a change to the evaluator, the basis or the counting
// rule shows up as a disagreement with the record rather than silently.

procedure test_nelson_filter()
    printf "The Noether-count stage, against the banked readings...";

    // [1] The stage is INERT unless asked. A pipeline run that has not opted in must not pay for
    // it, and must not be marked by it.
    D := 26; N := 3;
    X := CreateShimuraQuot(D, N, {1, 39});
    X`g := GenusShimuraCurveQuotient(D, N, X`W); X`CurveID := 1;
    cs := [X];
    before := [* assigned c`IsSubhyp select c`IsSubhyp else "unset" : c in cs *];
    FilterByNelsonCount(~cs);                       // not enabled: must do nothing
    assert [* assigned c`IsSubhyp select c`IsSubhyp else "unset" : c in cs *] eq before;

    // [2] The three banked readings, reproduced at the settings they were taken at.
    expect := [* <26, 3, {1, 39}, "HYPERELLIPTIC", 1>,
                 <34, 7, {1, 14}, "HYPERELLIPTIC", 1>,
                 <14, 11, {1, 11}, "NOT hyperelliptic", 0> *];
    for e in expect do
        DD := e[1]; NN := e[2]; WW := e[3]; want := e[4]; wantrel := e[5];
        Y := CreateShimuraQuot(DD, NN, WW);
        Y`g := GenusShimuraCurveQuotient(DD, NN, WW); Y`CurveID := 1;
        assert Y`g eq 3;
        // a negative reading is refused without a second precision, so give one
        conf := want eq "NOT hyperelliptic" select 45 else 0;
        rs := NelsonFilterReadings([Y] : Cutoff := 40, Precision := 30, Points := 18,
                                         ConfirmPrecision := conf);
        assert #rs eq 1;
        got := rs[1];
        if got[6] ne want or got[9] ne wantrel then
            printf "\n  X_0^%o(%o)/%o: read %o with %o relations, banked %o with %o\n",
                   DD, NN, Sort(Setseq(WW)), got[6], got[9], want, wantrel;
        end if;
        assert got[6] eq want;
        assert got[9] eq wantrel;
        printf " %o_%o:%o", DD, NN, got[9];
    end for;

    // [3] The stage marks nothing, whatever it reads. This is the property the header argues for:
    // the pipeline's verdict fields are for proved verdicts, and these readings are not proofs.
    cs2 := [X];
    FilterByNelsonCount(~cs2 : Enabled := true, Cutoff := 40, Precision := 30, MaxLevel := 100);
    assert not assigned cs2[1]`IsSubhyp;
    printf " Done!\n";
end procedure;

test_nelson_filter();
