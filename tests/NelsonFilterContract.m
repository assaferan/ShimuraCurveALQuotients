// What the Noether-count stage must do to the pipeline, as opposed to what it computes.
//
// The readings themselves are checked against the banked table in tests/_offline/NelsonFilter.m,
// which costs 2150 s. What is checked here is the part that protects everything else and must
// therefore run every time, in under a second:
//   * the stage is INERT unless asked, so an ordinary pipeline run neither pays for it nor is
//     changed by it;
//   * it marks NO curve even when it does run and does produce a reading. Every other filter sets
//     IsSubhyp, which the rest of the pipeline treats as proved; a Noether count is a numerical
//     kernel dimension and is not. The reading here is taken at a cutoff far too small to mean
//     anything, which is the point: the assertion is about the contract, not the verdict, and it
//     would be vacuous if the stage had not actually run.

procedure test_nelson_filter_contract()
    printf "The Noether-count stage: inert unless asked, and marks nothing...";
    D := 26; N := 3;
    X := CreateShimuraQuot(D, N, {1, 39});
    X`g := GenusShimuraCurveQuotient(D, N, X`W); X`CurveID := 1;
    assert X`g eq 3;

    // [1] not enabled: no reading, no mark, and no cost
    cs := [X];
    t := Cputime();
    FilterByNelsonCount(~cs);
    assert not assigned cs[1]`IsSubhyp;
    assert Cputime(t) lt 5;

    // [2] enabled, with a cutoff too small for the reading to mean anything: the stage runs and
    // returns a tuple, and the curve is still unmarked.
    cs2 := [X];
    FilterByNelsonCount(~cs2 : Enabled := true, Cutoff := 6, Precision := 20);
    assert not assigned cs2[1]`IsSubhyp;
    rs := NelsonFilterReadings([X] : Cutoff := 6, Precision := 20);
    assert #rs eq 1;
    assert rs[1][2] eq 26 and rs[1][3] eq 3 and rs[1][5] eq 3;
    assert rs[1][6] in {"HYPERELLIPTIC", "NOT hyperelliptic", "INCONCLUSIVE", "REFUSED"};

    // [3] the cost guards do what they say: nothing of too small a genus or too large a level
    assert IsEmpty(NelsonFilterReadings([X] : MinGenus := 4));
    assert IsEmpty(NelsonFilterReadings([X] : MaxLevel := 10));
    printf " Done!\n";
end procedure;

test_nelson_filter_contract();
