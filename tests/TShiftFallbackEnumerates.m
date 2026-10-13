// tests/TShiftFallbackEnumerates.m
//
// A spanning set the solver script derives from a lower rung and the shift set can fall short of
// the rank the ring of forms needs.  The search must then discard the derived file and enumerate
// the rung with Normaliz instead of climbing: with a cache holding only the enumerated (60, 3, 0)
// (one point), the level-60 shift set (one quotient, pole 8) and the two products these give for
// (60, 11, 0) -- a derived file of two points where the ring of {oo}-weakly holomorphic forms of
// level 60 at the first rung needs rank 6 -- the ring routine for X_0^6(5) must finish with the
// cache file (60, 11, 0) replaced by the lattice points of the polytope (the committed file, 162
// points, read as a set) and return the same basis as it does from the committed cache.
//
// The ring is computed in a separate Magma process started in a temporary directory, since the
// cache path is relative to the working directory.  Needs Normaliz (NORMALIZ_BIN, or normaliz on
// the path): the enumeration of (60, 11, 0) takes about a second.

// evaluated at top level: `eval` inside a procedure segfaults under the harness
committed := Sort(eval Read("polymake/polymake_solution_60_11_0"));
here := GetCurrentDirectory();
tmp := Pipe("mktemp -d", "");
tmp := tmp[1..#tmp-1];
assert #tmp gt 0 and tmp[1] eq "/";
ok := System(Sprintf("mkdir -p %o/polymake", tmp));
assert ok eq 0;
// the shift set: the weight-0 quotients of level 60 with poles at oo of order at most 8
rc := System(Sprintf("cd %o && python3 %o/nmzsolve.py 60 8 0 polymake/tshift_w0_60.txt 0 1 0 2> gen.err", tmp, here));
assert rc eq 0;
rc := System(Sprintf("cd %o && python3 %o/nmzsolve.py 60 3 0 polymake/polymake_solution_60_3_0 2>> gen.err", tmp, here));
assert rc eq 0;
// the derived (60, 11, 0): the products of the rung (60, 3, 0) with the shifts, two points
rc := System(Sprintf("cd %o && NORMALIZ_BIN=/nonexistent python3 %o/nmzsolve.py 60 11 0 polymake/polymake_solution_60_11_0 2>> gen.err", tmp, here));
assert rc eq 0;
first := Pipe(Sprintf("head -1 %o/polymake/polymake_solution_60_11_0", tmp), "");
assert first[1..10] eq "// derived";
short := eval Read(Sprintf("%o/polymake/polymake_solution_60_11_0", tmp));
assert #short lt 6;

// the ring routine, from the temporary cache and from the committed one, printing the exponents
// of the eta quotients it returns
script := Sprintf("AttachSpec(\"%o/ShimuraQuotients.spec\");\nSetVerbose(\"ShimuraQuotients\", 2);\nSetColumns(0);\n" cat
                  "E, n, n0, t, etas := WeaklyHolomorphicBasis(6, 5);\n" cat
                  "printf \"\\nBASIS %%o %%o %%o\\n\", n, n0, [Exponents(e) : e in etas];\nquit;\n", here);
Write(Sprintf("%o/whb.m", tmp), script : Overwrite := true);
// the search finds the solver script through NMZSOLVE when it is not in the working directory
rc := System(Sprintf("cd %o && NMZSOLVE=%o/nmzsolve.py magma -b whb.m < /dev/null > whb.log 2>&1", tmp, here));
assert rc eq 0;
log := Read(Sprintf("%o/whb.log", tmp));
rc := System(Sprintf("cd %o && magma -b %o/whb.m < /dev/null > ref.log 2>&1", here, tmp));
assert rc eq 0;
ref := Read(Sprintf("%o/ref.log", here));
_ := System(Sprintf("rm -f %o/ref.log", here));
after_first := Pipe(Sprintf("head -1 %o/polymake/polymake_solution_60_11_0", tmp), "");
after := Sort(eval Read(Sprintf("%o/polymake/polymake_solution_60_11_0", tmp)));
_ := System(Sprintf("rm -rf %o", tmp));

procedure test_TShiftFallbackEnumerates(log, ref, after_first, after, committed)
    printf "Testing that the search enumerates a rung whose derived set falls short ((60, 11, 0) for X_0^6(5))...";
    assert "enumerating the rung" in log;                        // the derived set was rejected
    assert after_first[1..10] ne "// derived";                    // and replaced by an enumeration
    assert after eq committed;                                    // the lattice points of the polytope
    basis_line := func<s | [l : l in Split(s, "\n") | #l ge 5 and l[1..5] eq "BASIS"]>;
    bl := basis_line(log); br := basis_line(ref);
    assert #bl eq 1 and #br eq 1 and bl[1] eq br[1];              // the same ring as from the committed cache
    printf "Done! (%o points enumerated)\n", #after;
end procedure;

test_TShiftFallbackEnumerates(log, ref, after_first, after, committed);
