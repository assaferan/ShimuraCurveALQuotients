// tests/TShiftFallbackEnumerates.m
//
// A spanning set the solver script derives from a lower rung and the shift set can fall short of
// the rank the ring of forms needs.  The search must then discard the derived file and enumerate
// the rung with Normaliz instead of climbing: with a cache holding only the enumerated (60, 3, 0)
// (one point), the level-60 shift set (one quotient, pole 8) and the two products these give for
// (60, 11, 0) -- a derived file of two points where the ring of {oo}-weakly holomorphic forms of
// level 60 at the first rung needs rank 9 -- the ring routine for X_0^6(5) must finish with the
// cache file (60, 11, 0) replaced by the lattice points of the polytope (the committed file, 162
// points, read as a set) and return the same eta quotients as it does from the committed cache.
//
// The cache path is relative to the working directory, so the ring is computed once here and once
// after changing into a temporary directory holding the small cache and a copy of the solver
// script; the working directory is restored afterwards.  Needs Normaliz (NORMALIZ_BIN, or normaliz
// on the path): the enumeration of (60, 11, 0) takes about a second.

// evaluated at top level: `eval` inside a procedure segfaults under the harness
committed := Sort(eval Read("polymake/polymake_solution_60_11_0"));
here := GetCurrentDirectory();
tmp := Pipe("mktemp -d", "");
tmp := tmp[1..#tmp-1];
assert #tmp gt 0 and tmp[1] eq "/";
ok := System(Sprintf("mkdir -p %o/polymake && cp %o/nmzsolve.py %o/", tmp, here, tmp));
assert ok eq 0;
// the shift set: the weight-0 quotients of level 60 with poles at oo of order at most 8
rc := System(Sprintf("cd %o && python3 nmzsolve.py 60 8 0 polymake/tshift_w0_60.txt 0 1 0 2> gen.err", tmp));
assert rc eq 0;
rc := System(Sprintf("cd %o && python3 nmzsolve.py 60 3 0 polymake/polymake_solution_60_3_0 2>> gen.err", tmp));
assert rc eq 0;
// the derived (60, 11, 0): the products of the rung (60, 3, 0) with the shifts, two points
rc := System(Sprintf("cd %o && python3 nmzsolve.py 60 11 0 polymake/polymake_solution_60_11_0 2>> gen.err", tmp));
assert rc eq 0;
first := Pipe(Sprintf("head -1 %o/polymake/polymake_solution_60_11_0", tmp), "");
assert first[1..10] eq "// derived";
short := eval Read(Sprintf("%o/polymake/polymake_solution_60_11_0", tmp));
assert #short lt 9;

// the ring from the committed cache, then from the temporary one
_, n_ref, n0_ref, _, etas_ref := WeaklyHolomorphicBasis(6, 5);
ChangeDirectory(tmp);
_, n_tmp, n0_tmp, _, etas_tmp := WeaklyHolomorphicBasis(6, 5);
ChangeDirectory(here);
after_first := Pipe(Sprintf("head -1 %o/polymake/polymake_solution_60_11_0", tmp), "");
after := Sort(eval Read(Sprintf("%o/polymake/polymake_solution_60_11_0", tmp)));
_ := System(Sprintf("rm -rf %o", tmp));

procedure test_TShiftFallbackEnumerates(after_first, after, committed, n_ref, n0_ref, etas_ref, n_tmp, n0_tmp, etas_tmp)
    printf "Testing that the search enumerates a rung whose derived set falls short ((60, 11, 0) for X_0^6(5))...";
    assert after_first[1..10] ne "// derived";                    // the derived set was replaced by an enumeration
    assert after eq committed;                                    // the lattice points of the polytope
    assert n_tmp eq n_ref and n0_tmp eq n0_ref;
    assert [Exponents(e) : e in etas_tmp] eq [Exponents(e) : e in etas_ref];   // the same ring as from the committed cache
    printf "Done! (%o points enumerated, %o eta quotients)\n", #after, #etas_ref;
end procedure;

test_TShiftFallbackEnumerates(after_first, after, committed, n_ref, n0_ref, etas_ref, n_tmp, n0_tmp, etas_tmp);
