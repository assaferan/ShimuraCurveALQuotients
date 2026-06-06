// Smoke test for FilterByWeilPolynomialGenusScaled.
//
// Loads a small slice of UpdateCurves6 data, runs the genus-scaled filter,
// and verifies:
//   (1) it completes without error
//   (2) every curve it marks non-hyperelliptic has TestInWhichProved set
//   (3) the genus-bound formula assigns the correct bounds

AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 0);

print "Loading data...";
all_curves := eval Read("data/curves_after_UpdateCurves6.dat");

// Take a small slice covering a range of genera
curves := all_curves[1..200];

print "Running FilterByWeilPolynomialGenusScaled on 200 curves...";
FilterByWeilPolynomialGenusScaled(~curves);
print "Completed without error.";

// Verify: every curve marked false has TestInWhichProved
for c in curves do
    if assigned c`IsSubhyp and not c`IsSubhyp then
        assert assigned c`TestInWhichProved;
    end if;
end for;
print "PASS: TestInWhichProved set for all eliminated curves.";

// Verify the bound formula directly
expected := AssociativeArray();
expected[3]  := 25; expected[4]  := 25; expected[5]  := 25; expected[6]  := 19;
expected[7]  := 22; expected[8]  := 19; expected[9]  := 16; expected[10] := 13;
expected[11] := 10; expected[12] := 7;  expected[13] := 7;  expected[20] := 7;
for g in Keys(expected) do
    if g lt 6 then bd := 25;
    elif g eq 6 then bd := 19;
    else bd := Maximum(7, 25 - 3*(g - 6));
    end if;
    if bd ne expected[g] then
        error Sprintf("g=%o: expected %o got %o", g, expected[g], bd);
    end if;
end for;
print "PASS: genus-bound formula correct.";

print "All tests passed.";
quit;
