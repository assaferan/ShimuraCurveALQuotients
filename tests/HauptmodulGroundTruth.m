// Regression test for HauptmodulM0Residuals -- the cheap route to the outer-m=0 multipliers.
//
// The tables below are the Hauptmodul rows as the pipeline produced them (find_signs_hauptmodul's
// own inputs), recorded verbatim so that this test is DETERMINISTIC: it exercises the solver, not
// the CM machinery, and so is immune to the randomised quaternion-order representative that makes
// the Schofer values themselves vary between runs.
//
// The decisive check is X0^15(2), whose multipliers are known independently from the vector-valued
// oracle (4 on row -1, 2 on row -2): the solver must recover BOTH, from a table in which the
// normalising discriminants are themselves non-firing -- the regime opposite to X0^14(5)'s.

procedure test_determined_bases()
    printf "  multipliers recovered on four bases...";

    // ---- X0^14(5), N = 5.  Normalisers -91 (s) and -11 (stilde) both FIRE, so the information is
    // at the non-firing -35 and -280.  Nothing is applied here (odd N), so residual = multiplier.
    ds  := [ -4, -35, -280, -11, -91, -84, -51 ];
    dg  := [  1,   1,    1,   1,   1,   1,   2 ];
    s   := [* Infinity(), 5/8,   55/128,  0, 1, 1/4, 1/8 *];
    st  := [* Infinity(), 35/8, 585/128,  1, 0, 3/4, 9/4 *];
    assert HauptmodulM0Residuals(s, st, ds, dg, 5) eq { <1, 1> };

    // ---- X0^15(2), N = 2.  Normalisers -120 (s) and -40 (stilde) are NON-firing, so the
    // information is at the firing -15 and -7 and the exponents flip sign.  The pipeline applies
    // m0_multiplier = 4 to both rows here, and the true multipliers are 4 and 2 -- hence residuals
    // 0 and -2.  Three further points (-52, -88, -132) carry no information and must simply hold.
    ds  := [ -12, -40, -120, -52, -88, -132, -15, -7 ];
    dg  := [   1,   1,    1,   1,   1,    1,   1,  1 ];
    s   := [* Infinity(), 0,   8,   4,   16,   4, 5,   1 *];
    st  := [* Infinity(), 2/3, 0, 1/3,  2/3,   1, 1, 7/3 *];
    assert HauptmodulM0Residuals(s, st, ds, dg, 2) eq { <0, -2> };

    // ---- X0^6(5) and X0^10(3), N = 5 and 3.  Both come out 0, matching the oracle.  X0^6(5) is
    // the base whose signs are (-,+): a solver that folds the signs into the unknowns reports the
    // impossible "-1" here instead of 0.
    ds  := [ -4, -24, -40, -120, -19, -115, -51, -84 ];
    dg  := [  1,   1,   1,    1,   1,    1,   1,   1 ];
    s   := [* Infinity(), 0, 1, 3/128, 1/8, 23/648, 9/8,  3/2 *];
    st  := [* Infinity(), 1, 9, 19/16,   0, 104/81,   8,   13 *];
    assert HauptmodulM0Residuals(s, st, ds, dg, 5) eq { <0, 0> };

    ds  := [ -3, -8, -20, -120, -35, -123, -68, -155 ];
    dg  := [  1,  1,   1,    1,   1,    1,   2,    2 ];
    s   := [* 1, 0, 1/5,  8/35, Infinity(),  49/76, 4/19, 7/200 *];
    st  := [* 4, 1,   0,   1/7, Infinity(), 169/76, 16/19, 5/14 *];
    assert HauptmodulM0Residuals(s, st, ds, dg, 3) eq { <0, 0> };

    printf " ok\n";
end procedure;

procedure test_signs_are_per_discriminant()
    printf "  the per-discriminant signs are honoured...";
    // X0^6(5) again, restricted to its three non-firing points.  The relation holds there with
    // signs (-,+): -8 + 9 = 1, -3/16 + 19/16 = 1, -23/81 + 104/81 = 1.  Reading those as (+,+)
    // forces u = -1, which is not a power of N, and the multiplier would come out "no solution".
    for pair in [ <8, 9>, <3/16, 19/16>, <23/81, 104/81> ] do
        assert -pair[1] + pair[2] eq 1;
        assert pair[1] + pair[2] ne 1;
    end for;
    printf " ok\n";
end procedure;

procedure test_not_measurable_is_reported_not_guessed()
    printf "  an uninformative table returns no solution rather than a wrong one...";
    // Every discriminant shares the normalisers' firing status (N = 7 divides none of these), so the
    // correction cancels everywhere and nothing can be measured.  The answer must be "empty", not a
    // spurious pair.
    ds := [ -4, -8, -3, -11, -19 ];
    dg := [  1,  1,  1,   1,   1 ];
    s  := [* Infinity(), 0, 1, 1/2, 1/4 *];
    st := [* Infinity(), 1, 0, 1/2, 3/4 *];
    assert HauptmodulM0Residuals(s, st, ds, dg, 7) eq {};
    printf " ok\n";
end procedure;

procedure test_HauptmodulGroundTruth()
    printf "Testing the Hauptmodul-consistency route to the m=0 multipliers...\n";
    test_determined_bases();
    test_signs_are_per_discriminant();
    test_not_measurable_is_reported_not_guessed();
    printf "Done!\n";
end procedure;

test_HauptmodulGroundTruth();
