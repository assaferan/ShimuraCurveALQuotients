// tests/DataFileNoWrap.m
//
// The pipeline writes its curve lists with Sprint(curves, "Magma").  With Magma's default line
// width of 80 columns a long verdict label is split over several lines with "\n" inserted, and the
// label read back from the file is not the one written.  The writers now set the line width to 0
// first.
// This test writes a curve with a label longer than 80 characters the way the writers do, reads
// it back, and checks the label is unchanged; it also checks that the default width does damage
// it, so that the test tells the two apart.

printf "Testing that data files are written without line wrapping...";
X := CreateShimuraQuot(6, 23, {Integers()|1, 23});
X`g := 3;
X`IsSubhyp := false; X`IsHyp := false;
X`TestInWhichProved := "SpecialFiberIsomorphism, isomorphic over F_2 to curve 169 (source: WeilPolynomial with p = 2, a label long enough to be wrapped at eighty columns)";
assert #X`TestInWhichProved gt 80;
fn := GetTempDir() cat "/nowrap_test.dat";

SetColumns(80);
Write(fn, Sprint([X], "Magma") : Overwrite);
Y := eval Read(fn);
assert Y[1]`TestInWhichProved ne X`TestInWhichProved and Position(Y[1]`TestInWhichProved, "\n") gt 0;

SetColumns(0);
Write(fn, Sprint([X], "Magma") : Overwrite);
Y := eval Read(fn);
assert Y[1]`TestInWhichProved eq X`TestInWhichProved;
assert Position(Read(fn), "\\n") eq 0;

// the three writers and the in-process runner set the width before printing
for f in ["run_sequential_stage.m", "parallel_filter_worker.m", "parallel_merge.m", "workingcode.m"] do
    t := Read(f);
    assert Position(t, "SetColumns(0)") gt 0 and Position(t, "SetColumns(0)") lt Position(t, "Sprint(");
end for;
printf " ok\n";
