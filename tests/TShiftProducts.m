// tests/TShiftProducts.m
//
// The solver script's shift fallback for a tall rung with poles at infinity only: with the
// enumerated rung (380, 92, 0) and the level-380 weight-0 shift set in a cache, a request for
// (380, 158, 0) is answered by the products of the rung with the shifts (no Normaliz), every one of
// which is a lattice point of the enumerated (380, 158, 0), and the products span the same space of
// forms: the rank of the q-expansions at infinity equals that of the full rung.  The file carries
// the first-line mark of a derived spanning set.
//
// The two rungs are the committed cache files (Normaliz enumerations).  The rank is computed here
// twice with plain power series, independently of the eta-quotient code, and the two numbers are
// compared with each other, not with a stored constant.

// evaluated at top level: `eval` inside a procedure segfaults under the harness
M := 380; nlo := 92; nhi := 158;
full := Sort(eval Read(Sprintf("polymake/polymake_solution_%o_%o_0", M, nhi)));
tmp := Pipe("mktemp -d", "");
tmp := tmp[1..#tmp-1];
assert #tmp gt 0 and tmp[1] eq "/";
ok := System(Sprintf("mkdir -p %o/polymake && cp polymake/polymake_solution_%o_%o_0 polymake/tshift_w0_%o.txt %o/polymake/",
                     tmp, M, nlo, M, tmp));
assert ok eq 0;
here := GetCurrentDirectory();
rc := System(Sprintf("cd %o && NORMALIZ_BIN=/nonexistent python3 %o/nmzsolve.py %o %o 0 polymake/polymake_solution_%o_%o_0 2> nmzsolve.err",
                     tmp, here, M, nhi, M, nhi));
assert rc eq 0;                                              // served from the cache, no Normaliz
out := Sprintf("%o/polymake/polymake_solution_%o_%o_0", tmp, M, nhi);
first := Pipe(Sprintf("head -1 %o", out), "");
shifted := Sort(eval Read(out));
_ := System(Sprintf("rm -rf %o", tmp));

procedure test_TShiftProducts(M, nhi, full, shifted, first)
    printf "Testing the shift fallback, (380, 92, 0) x shifts -> (380, 158, 0)...";
    ds := Divisors(M);
    assert first[1..10] eq "// derived";                         // marked as a spanning set
    assert #shifted gt 0 and #shifted lt #full;
    assert forall{r : r in shifted | r in full};                  // every product is a lattice point
    // rank of the q-expansions at infinity over GF(p): prod_d eta(d tau)^{r_d} = q^{(sum d r_d)/24} prod (1 - q^{dn})^{r_d};
    // a weight-1/2 form holomorphic away from infinity with a pole of order at most nhi there has
    // at most index/24 + nhi zeros at infinity, so this many coefficients determine it
    index := Integers()!(M * &*[1 + 1/p : p in PrimeDivisors(M)]);
    T := index div 24 + nhi + 2;
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
    printf "Done! (%o products of %o points, rank %o = %o)\n", #shifted, #full, rk_shift, rk_full;
end procedure;

test_TShiftProducts(M, nhi, full, shifted, first);
