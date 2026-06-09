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

// Verify the per-curve database-tight prime bound: WeilClassNumberPrimeBound(Qmax, g)
// is the largest b with 4*Qmax*b^g <= 2^40 (the range covered by the class-group tables).
DB := ClassNumberDataMaxAbsDisc();
for g in [3..10] do
    for Qmax in [1, 2, 37, 1000, 14430, 2*10^6] do
        b := WeilClassNumberPrimeBound(Qmax, g);
        assert 4*Qmax*b^g le DB;          // bound stays within the tables
        assert 4*Qmax*(b+1)^g gt DB;      // and is maximal
    end for;
end for;
print "PASS: WeilClassNumberPrimeBound is database-tight.";

print "All tests passed.";
quit;
