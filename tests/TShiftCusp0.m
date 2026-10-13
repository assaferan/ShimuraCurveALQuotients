// tests/TShiftCusp0.m
//
// The solver script's shift fallback at a cusp-0 order m > 0: with only the enumerated rung
// (204, 20, 32) and the level-204 weight-0 shift set in a cache, a request for (204, 45, 32) is
// answered by shifted points, every one of which is a lattice point of the enumerated (204, 45, 32),
// and the shifted set spans the same space of forms: the rank of the q-expansions at oo equals that
// of the full rung.  The output file carries the first-line mark of a derived spanning set.
//
// PROVENANCE.  The two rungs are the committed cache files (Normaliz enumerations).  The rank is
// computed here twice with plain power series, independently of the EtaQuot code, and the two
// numbers are compared with each other, not with a stored constant; the value, 63, matches the one
// Sachi computed in the review of PR #77 (63 = 63) and the campaign-branch check
// vvdata/weyl-campaign/cusp0ladder/fbrank.m.

// The cache files are evaluated at top level: `eval` inside a procedure segfaults under the harness.
M := 204; nlo := 20; nhi := 45; m := 32;
full := Sort(eval Read(Sprintf("polymake/polymake_solution_%o_%o_%o", M, nhi, m)));
tmp := Pipe("mktemp -d", "");
tmp := tmp[1..#tmp-1];                                       // strip the newline
assert #tmp gt 0 and tmp[1] eq "/";
ok := System(Sprintf("mkdir -p %o/polymake && cp polymake/polymake_solution_%o_%o_%o polymake/tshift_w0_%o.txt %o/polymake/",
                     tmp, M, nlo, m, M, tmp));
assert ok eq 0;
here := GetCurrentDirectory();
rc := System(Sprintf("cd %o && NORMALIZ_BIN=/nonexistent python3 %o/nmzsolve.py %o %o %o polymake/polymake_solution_%o_%o_%o 2> nmzsolve.err",
                     tmp, here, M, nhi, m, M, nhi, m));
assert rc eq 0;                                              // served from the cache, no Normaliz
out := Sprintf("%o/polymake/polymake_solution_%o_%o_%o", tmp, M, nhi, m);
first := Pipe(Sprintf("head -1 %o", out), "");
shifted := Sort(eval Read(out));
_ := System(Sprintf("rm -rf %o", tmp));

procedure test_TShiftCusp0(M, nhi, m, full, shifted, first)
    printf "Testing the shift fallback at a cusp-0 order, (204, 20, 32) -> (204, 45, 32)...";
    ds := Divisors(M);
    assert first[1..10] eq "// derived";                         // marked as a spanning set
    assert #shifted gt 0 and #shifted lt #full;
    assert forall{r : r in shifted | r in full};                  // every shifted point is a lattice point

    // rank of the q-expansions at oo over GF(p): prod_d eta(d tau)^{r_d} = q^{(sum d r_d)/24} prod (1 - q^{dn})^{r_d}
    index := Integers()!(M * &*[1 + 1/p : p in PrimeDivisors(M)]);
    T := index div 24 + nhi + m + 2;
    function rank_of(pts)
        F := GF(1000003);
        R<q> := PowerSeriesRing(F, T + nhi + 1);
        E := [&*[R | 1 - q^(d*n) : n in [1..(T + nhi) div d]] : d in ds];
        rows := [];
        for r in pts do
            v24 := &+[ds[i]*r[i] : i in [1..#ds]];
            assert v24 mod 24 eq 0 and v24 div 24 ge -nhi;
            f := &*[R | E[i]^r[i] : i in [1..#ds] | r[i] ne 0];
            f := q^(v24 div 24 + nhi) * f;
            Append(~rows, [Coefficient(f, j) : j in [0..T + nhi]]);
        end for;
        return Rank(Matrix(F, rows));
    end function;
    rk_full := rank_of(full);
    rk_shift := rank_of(shifted);
    assert rk_full eq rk_shift;
    printf "Done! (%o shifted points of %o, rank %o = %o)\n", #shifted, #full, rk_shift, rk_full;
end procedure;

test_TShiftCusp0(M, nhi, m, full, shifted, first);
